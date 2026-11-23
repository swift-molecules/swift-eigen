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
        } else if least > 0 {
            shift = -least / 1000
        } else if least < 0 {
            shift = least * 2
        } else {
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
        var slices = 0
        let budget = 64 + 4 * count
        while verified < count {
            slices += 1
            guard slices <= budget else { throw .incomplete }
            let pending = values.reduce(0) { $1 > bound ? $0 + 1 : $0 }
            let wanted = Swift.max(count - verified - pending, 1)
            let deflation = values.indices.filter { values[$0] > bound || values[$0] >= sigma - reach }
            let size = Swift.min(
                options.basis.map { Swift.max($0, 1) } ?? Swift.min(Swift.max(2 * wanted + 20, 40), 300) * extra,
                order - deflation.count
            )
            guard size > 0 else { throw .incomplete }
            let run = self.run(
                shift: sigma,
                factorization: factorization,
                deflation: deflation,
                bound: bound,
                wanted: wanted,
                size: size
            )
            if run.accepted.isEmpty {
                extra *= 2
                continue
            }
            for pair in run.accepted {
                values.append(pair.value)
                vectors.append(pair.vector)
                reach = Swift.max(reach, 2 * Swift.abs(pair.value - sigma))
            }
            // The next shift: midway between the highest eigenvalue found below the first unconverged Ritz value and
            // the next value above it.
            let unconverged = run.unconverged.min() ?? Scalar.infinity
            let above = values.filter { $0 > bound }.sorted()
            let settled = above.last { $0 < unconverged } ?? above[0]
            let next = Swift.min(above.first { $0 > settled } ?? Scalar.infinity, unconverged)
            let candidate: Scalar
            if next.isFinite {
                candidate = settled + (next - settled) / 2
            } else {
                let step = Swift.max(settled - bound, Swift.abs(settled) * Scalar.ulpOfOne.squareRoot())
                candidate = settled + step / 2
            }
            let shifted = try factor(candidate)
            let below = values.reduce(0) { $1 <= shifted.shift ? $0 + 1 : $0 }
            let negative = shifted.factorization.inertia.negative
            if negative == below {
                bound = shifted.shift
                verified = negative
            } else if negative < below {
                throw .incomplete
            }
            sigma = shifted.shift
            factorization = shifted.factorization
        }
        let ascending = values.indices.sorted { values[$0] < values[$1] }.prefix(count)
        return ascending.map { Eigen.Pair(value: values[$0], vector: vectors[$0]) }
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
        deflation: [Int],
        bound: Scalar,
        wanted: Int,
        size: Int
    ) -> Run {
        let n = order
        let tolerance = options.tolerance
        stream += 1
        var w = Self.random(count: n, seed: options.seed, stream: stream)
        var basis: [[Scalar]] = []
        basis.reserveCapacity(size)
        var alphas: [Scalar] = []
        var betas: [Scalar] = []
        var q = [Scalar](repeating: 0, count: n)
        b.apply(to: w, into: &q)
        orthogonalize(&w, &q, deflation: deflation, basis: basis)
        var beta = Self.dot(w, q).squareRoot()
        guard beta > 0 else { return Run(accepted: [], unconverged: []) }
        Self.scale(&w, by: 1 / beta)
        Self.scale(&q, by: 1 / beta)
        var p = q
        var steps = 0
        while steps < size {
            let v = w
            basis.append(v)
            factorization.solve(p, into: &w)
            let alpha = Self.dot(w, p)
            alphas.append(alpha)
            Self.add(&w, -alpha, v)
            if steps > 0 { Self.add(&w, -betas[steps - 1], basis[steps - 1]) }
            b.apply(to: w, into: &q)
            orthogonalize(&w, &q, deflation: deflation, basis: basis)
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
        var ax = [Scalar](repeating: 0, count: n)
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
            for column in 0..<steps {
                Self.add(&x, ritz.vectors[index * steps + column], basis[column])
            }
            b.apply(to: x, into: &bx)
            let length = Self.dot(x, bx).squareRoot()
            Self.scale(&x, by: 1 / length)
            Self.scale(&bx, by: 1 / length)
            a.apply(to: x, into: &ax)
            var residual: Scalar = 0
            for row in 0..<n {
                let r = ax[row] - value * bx[row]
                residual += r * r
            }
            let bound = tolerance * (normA + Swift.abs(value) * normB) * Self.dot(x, x).squareRoot()
            if residual.squareRoot() <= bound {
                accepted.append(Eigen.Pair(value: value, vector: x))
            } else {
                unconverged.append(value)
            }
        }
        return Run(accepted: accepted, unconverged: unconverged)
    }

    /// Two passes of classical Gram–Schmidt in the `B` inner product against the deflated eigenvectors and the basis;
    /// `q` is `B w` on entry and on exit.
    @inlinable
    func orthogonalize(_ w: inout [Scalar], _ q: inout [Scalar], deflation: [Int], basis: [[Scalar]]) {
        for _ in 0..<2 {
            var coefficients: [Scalar] = []
            coefficients.reserveCapacity(deflation.count + basis.count)
            for index in deflation { coefficients.append(Self.dot(vectors[index], q)) }
            for vector in basis { coefficients.append(Self.dot(vector, q)) }
            var position = 0
            for index in deflation {
                Self.add(&w, -coefficients[position], vectors[index])
                position += 1
            }
            for vector in basis {
                Self.add(&w, -coefficients[position], vector)
                position += 1
            }
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

    @inlinable
    static func dot(_ x: [Scalar], _ y: [Scalar]) -> Scalar {
        var s0: Scalar = 0
        var s1: Scalar = 0
        var s2: Scalar = 0
        var s3: Scalar = 0
        let n = x.count
        var index = 0
        while index + 4 <= n {
            s0 += x[index] * y[index]
            s1 += x[index + 1] * y[index + 1]
            s2 += x[index + 2] * y[index + 2]
            s3 += x[index + 3] * y[index + 3]
            index += 4
        }
        while index < n {
            s0 += x[index] * y[index]
            index += 1
        }
        return (s0 + s1) + (s2 + s3)
    }

    /// `y += factor · x`.
    @inlinable
    static func add(_ y: inout [Scalar], _ factor: Scalar, _ x: [Scalar]) {
        if factor == 0 { return }
        for index in y.indices { y[index] += factor * x[index] }
    }

    @inlinable
    static func scale(_ x: inout [Scalar], by factor: Scalar) {
        for index in x.indices { x[index] *= factor }
    }
}
