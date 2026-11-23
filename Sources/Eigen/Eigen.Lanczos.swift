internal import Threefry

extension Eigen {
    /// Shift-invert Lanczos on `(A − σB)⁻¹ B`, self-adjoint in the `B` inner product, sliced over shifts and checked
    /// by the inertia of `A − σB`.
    @usableFromInline
    internal struct Lanczos<Scalar: BinaryFloatingPoint> {
        @usableFromInline let a: Sparse.Symmetric<Scalar>
        @usableFromInline let b: Sparse.Symmetric<Scalar>
        @usableFromInline let options: Eigen.Pencil<Scalar>.Options
        @usableFromInline let order: Int
        @usableFromInline let elimination: Sparse.Elimination
        @usableFromInline let normA: Scalar
        @usableFromInline let normB: Scalar
        /// Every eigenpair found so far, `B`-normalized, in the order found.
        @usableFromInline var values: [Scalar] = []
        @usableFromInline var vectors: [[Scalar]] = []
        @usableFromInline var stream: UInt64 = 0

        @inlinable
        init(pencil: Eigen.Pencil<Scalar>, options: Eigen.Pencil<Scalar>.Options) throws(Eigen.Error) {
            self.a = pencil.a
            self.b = pencil.b
            self.options = options
            self.order = pencil.a.order
            let union: Sparse.Symmetric<Scalar>
            do throws(Sparse.Error) {
                union = try pencil.a.adding(pencil.b, scaledBy: 0)
            } catch {
                throw .order
            }
            self.elimination = Sparse.Elimination(union.pattern)
            self.normA = pencil.a.norm
            self.normB = pencil.b.norm
        }
    }
}

extension Eigen.Lanczos {
    @usableFromInline
    typealias Factorization = Sparse.Symmetric<Scalar>.Factorization

    /// The factorization of `A − σB`, nudging `σ` off an exact eigenvalue of a leading block when a pivot is zero.
    @inlinable
    func factor(_ shift: Scalar) throws(Eigen.Error) -> (shift: Scalar, factorization: Factorization) {
        var sigma = shift
        let scale = Swift.max(normA / Swift.max(normB, Scalar.leastNormalMagnitude), Scalar.leastNormalMagnitude)
        for attempt in 0..<8 {
            do throws(Sparse.Error) {
                let matrix = try a.adding(b, scaledBy: -sigma)
                return (sigma, try Factorization(matrix, elimination: elimination))
            } catch {
                sigma += Scalar(attempt + 1) * Scalar.ulpOfOne.squareRoot() * Swift.max(Swift.abs(sigma), scale) * 1e-3
            }
        }
        throw .shift
    }

    /// A shift with no eigenvalue below it, and its factorization.
    @inlinable
    func start() throws(Eigen.Error) -> (shift: Scalar, factorization: Factorization) {
        let diagonalA = a.diagonal
        let diagonalB = b.diagonal
        var least = Scalar.infinity
        for index in 0..<order {
            guard diagonalB[index] > 0 else { throw .shift }
            least = Swift.min(least, diagonalA[index] / diagonalB[index])
        }
        let scale = normA / normB
        var shift: Scalar
        if let given = options.shift {
            shift = given
        } else if least > 0, let zero = try? factor(0), zero.shift == 0, zero.factorization.inertia.negative == 0 {
            // A positive definite `A`: zero lies below the spectrum and as close to its low end as any shift can.
            return zero
        } else if least < 0 {
            shift = least * 2
        } else {
            // A singular or nearly singular `A`, as with rigid-body modes: just below zero, clear of its rounding.
            shift = -scale * Scalar.ulpOfOne.squareRoot()
        }
        for _ in 0..<64 {
            let found = try factor(shift)
            if found.factorization.inertia.negative == 0 { return found }
            shift = found.shift - Swift.max(Swift.abs(found.shift), scale * Scalar.ulpOfOne) * 9
        }
        throw .shift
    }

    @inlinable
    mutating func lowest(_ count: Int) throws(Eigen.Error) -> [Eigen.Pair<Scalar>] {
        var (sigma, factorization) = try start()
        var bound = sigma
        var verified = 0
        var reach: Scalar = 0
        var extra = 1
        var missing = 0
        var stalled = 0
        var check: (shift: Scalar, negative: Int)? = nil
        let separation = options.tolerance.squareRoot()
        while verified < count {
            // A run deflates the eigenvectors found near its shift; the rest were purged from its start vector and
            // return only through rounding, which acceptance checks for.
            let nearby = values.indices.filter { Swift.abs(values[$0] - sigma) <= reach }
            let pending = values.reduce(0) { $1 > bound ? $0 + 1 : $0 }
            let wanted = Swift.max(count - verified - pending, missing, 1)
            let size = Swift.min(
                options.basis.map { Swift.max($0, 1) } ?? Swift.min(Swift.max(2 * wanted + 20, 40), 300) * extra,
                order - values.count
            )
            guard size > 0 else { throw .incomplete }
            let run = self.run(
                shift: sigma,
                factorization: factorization,
                nearby: nearby,
                bound: bound,
                wanted: wanted,
                size: size
            )
            var accepted = 0
            for pair in run.accepted where accept(pair.value, pair.vector, nearby: nearby, separation: separation) {
                accepted += 1
                reach = Swift.max(reach, Swift.abs(pair.value - sigma))
            }
            if accepted == 0 { extra = Swift.min(extra * 2, 64) }
            // The check: midway between the highest value found below the first unconverged Ritz value and the value
            // after it, never inside a cluster closer than rounding can tell apart. Its inertia counts the values
            // below it.
            let unconverged = run.unconverged.min() ?? Scalar.infinity
            let above = values.filter { $0 > bound }.sorted() + [unconverged]
            var settled = -1
            for index in 0..<(above.count - 1) where above[index] < unconverged {
                let gap = above[index + 1] - above[index]
                let scale = Swift.max(Swift.abs(above[index]), Swift.abs(above[index + 1]))
                if !above[index + 1].isFinite || gap > separation * scale { settled = index }
            }
            var next: (shift: Scalar, factorization: Factorization)? = nil
            if settled >= 0 {
                let candidate: Scalar
                if above[settled + 1].isFinite {
                    candidate = above[settled] + (above[settled + 1] - above[settled]) / 2
                } else {
                    let step = Swift.max(above[settled] - bound, Swift.abs(above[settled]) * separation)
                    candidate = above[settled] + step / 2
                }
                if check?.shift != candidate {
                    let shifted = try factor(candidate)
                    next = shifted
                    check = (shifted.shift, shifted.factorization.inertia.negative)
                }
            }
            guard let (shift, negative) = check else {
                stalled += 1
                guard stalled <= 24 else { throw .incomplete }
                continue
            }
            let below = values.reduce(0) { $1 <= shift ? $0 + 1 : $0 }
            if negative == below {
                // Everything below the check is found: it becomes the bound, and the next shift.
                bound = shift
                verified = negative
                missing = 0
                stalled = 0
                if let next {
                    sigma = next.shift
                    factorization = next.factorization
                } else {
                    (sigma, factorization) = try factor(shift)
                }
                reach = 0
                check = nil
            } else {
                // Copies of repeated values are still missing below the check. A run at the same shift from a start
                // vector purged of every value found has them as its nearest; after two such runs the shift moves
                // up to the check anyway. More found than counted comes only from a value within rounding of it.
                missing = Swift.max(negative - below, 0)
                stalled += 1
                guard stalled <= 24 else { throw .incomplete }
                if stalled % 3 == 0 {
                    if let next {
                        sigma = next.shift
                        factorization = next.factorization
                    } else {
                        (sigma, factorization) = try factor(shift)
                    }
                    reach = 0
                }
            }
        }
        let ascending = values.indices.sorted { values[$0] < values[$1] }.prefix(count)
        return ascending.map { Eigen.Pair(value: values[$0], vector: vectors[$0]) }
    }

    /// Adds a converged pair to those found unless it repeats one: a pair whose eigenvalue lies within rounding of a
    /// found one that was not deflated is `B`-orthogonalized against it, dropped if little is left, and its residual
    /// tested again.
    @inlinable
    mutating func accept(_ value: Scalar, _ vector: [Scalar], nearby: [Int], separation: Scalar) -> Bool {
        var x = vector
        let close = values.indices.filter { index in
            Swift.abs(values[index] - value) <= separation * Swift.max(Swift.abs(value), Swift.abs(values[index]))
        }
        if !close.isEmpty && close.contains(where: { !nearby.contains($0) }) {
            var bx = b.applied(to: x)
            for index in close {
                let overlap = Self.dot(vectors[index], bx)
                Self.add(&x, -overlap, vectors[index])
            }
            b.apply(to: x, into: &bx)
            let length = Self.dot(x, bx).squareRoot()
            guard length > 1 / 2 else { return false }
            Self.scale(&x, by: 1 / length)
            guard residual(value, x) else { return false }
        }
        values.append(value)
        vectors.append(x)
        return true
    }

    /// Whether `‖A x − λ B x‖₂ ≤ tolerance · (‖A‖ + |λ| ‖B‖) · ‖x‖₂`.
    @inlinable
    func residual(_ value: Scalar, _ x: [Scalar]) -> Bool {
        let ax = a.applied(to: x)
        let bx = b.applied(to: x)
        var sum: Scalar = 0
        for row in 0..<order {
            let r = ax[row] - value * bx[row]
            sum += r * r
        }
        return sum.squareRoot() <= options.tolerance * (normA + Swift.abs(value) * normB) * Self.dot(x, x).squareRoot()
    }
}

extension Eigen.Lanczos {
    @usableFromInline
    struct Run {
        @usableFromInline var accepted: [Eigen.Pair<Scalar>]
        @usableFromInline var unconverged: [Scalar]

        @inlinable
        init(accepted: [Eigen.Pair<Scalar>], unconverged: [Scalar]) {
            self.accepted = accepted
            self.unconverged = unconverged
        }
    }

    /// One Lanczos run at `shift`: the eigenpairs above `bound` whose residual meets the tolerance, and the Ritz values
    /// above `bound` that do not yet.
    @inlinable
    mutating func run(
        shift: Scalar,
        factorization: Factorization,
        nearby: [Int],
        bound: Scalar,
        wanted: Int,
        size: Int
    ) -> Run {
        let n = order
        let tolerance = options.tolerance
        stream += 1
        var w = Self.random(count: n, seed: options.seed, stream: stream)
        var basis = [Scalar](repeating: 0, count: n * size)
        var alphas: [Scalar] = []
        var betas: [Scalar] = []
        var q = [Scalar](repeating: 0, count: n)
        b.apply(to: w, into: &q)
        orthogonalize(&w, &q, against: Array(values.indices), basis: basis, columns: 0)
        var beta = Self.dot(w, q).squareRoot()
        guard beta > 0 else { return Run(accepted: [], unconverged: []) }
        Self.scale(&w, by: 1 / beta)
        Self.scale(&q, by: 1 / beta)
        var p = q
        var steps = 0
        while steps < size {
            basis.replaceSubrange((steps * n)..<((steps + 1) * n), with: w)
            factorization.solve(p, into: &w)
            let alpha = Self.dot(w, p)
            alphas.append(alpha)
            b.apply(to: w, into: &q)
            orthogonalize(&w, &q, against: nearby, basis: basis, columns: steps + 1)
            beta = Self.dot(w, q).squareRoot()
            steps += 1
            if !(beta > Scalar.ulpOfOne * Swift.abs(alpha)) {
                beta = 0
                break
            }
            betas.append(beta)
            if steps == size { break }
            Self.scale(&w, by: 1 / beta)
            Self.scale(&q, by: 1 / beta)
            p = q
            if steps >= wanted && steps % 10 == 0 {
                let estimate = Eigen.Tridiagonal.decomposed(
                    diagonal: alphas,
                    offDiagonal: Array(betas.prefix(steps - 1)),
                    rows: 1
                )
                var converged = 0
                for index in 0..<steps {
                    let theta = estimate.values[index]
                    guard theta != 0, shift + 1 / theta > bound else { continue }
                    if Swift.abs(beta * estimate.vectors[index]) <= tolerance * Swift.abs(theta) { converged += 1 }
                }
                if converged >= wanted { break }
            }
        }
        let ritz = Eigen.Tridiagonal.decomposed(
            diagonal: alphas,
            offDiagonal: Array(betas.prefix(steps - 1)),
            rows: steps
        )
        var accepted: [Eigen.Pair<Scalar>] = []
        var unconverged: [Scalar] = []
        var bx = [Scalar](repeating: 0, count: n)
        for index in 0..<steps {
            let theta = ritz.values[index]
            guard theta != 0 else { continue }
            let value = shift + 1 / theta
            guard value > bound else { continue }
            let last = ritz.vectors[index * steps + steps - 1]
            guard Swift.abs(beta * last) <= 10 * tolerance * Swift.abs(theta) else {
                unconverged.append(value)
                continue
            }
            var x = [Scalar](repeating: 0, count: n)
            Self.accumulate(&x, basis, Array(ritz.vectors[(index * steps)..<((index + 1) * steps)]))
            b.apply(to: x, into: &bx)
            Self.scale(&x, by: 1 / Self.dot(x, bx).squareRoot())
            if residual(value, x) {
                accepted.append(Eigen.Pair(value: value, vector: x))
            } else {
                unconverged.append(value)
            }
        }
        return Run(accepted: accepted, unconverged: unconverged)
    }

    /// Two passes of classical Gram–Schmidt in the `B` inner product against the found eigenvectors at `indices` and
    /// the first `columns` columns of the basis; `q` is `B w` on entry and on exit.
    @inlinable
    func orthogonalize(_ w: inout [Scalar], _ q: inout [Scalar], against indices: [Int], basis: [Scalar], columns: Int) {
        for _ in 0..<2 {
            let found = indices.map { Self.dot(vectors[$0], q) }
            let coefficients = Self.products(basis, columns: columns, q)
            for (position, index) in indices.enumerated() { Self.add(&w, -found[position], vectors[index]) }
            Self.accumulate(&w, basis, coefficients.map { -$0 })
            b.apply(to: w, into: &q)
        }
    }

    /// `count` draws from [−1, 1): the start vector of the run on `stream`.
    @usableFromInline
    static func random(count: Int, seed: UInt64, stream: UInt64) -> [Scalar] {
        let key = Threefry.Key(seed: seed)
        return (0..<count).map { index in
            Threefry.unit(key: key, stream: stream, index: UInt64(index), as: Scalar.self) * 2 - 1
        }
    }
}
