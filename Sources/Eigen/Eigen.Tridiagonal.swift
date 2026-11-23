extension Eigen {
    /// A symmetric tridiagonal matrix: its diagonal and the entries just beside it.
    public struct Tridiagonal<Scalar: BinaryFloatingPoint> {
        /// The diagonal, `order` entries.
        public let diagonal: [Scalar]
        /// The entries at `(i, i + 1)` and `(i + 1, i)`, `order − 1` of them.
        public let subdiagonal: [Scalar]

        /// Fails with `.order` unless there is one fewer off-diagonal entry than diagonal ones (none for an empty matrix).
        @inlinable
        public init(diagonal: [Scalar], subdiagonal: [Scalar]) throws(Eigen.Error) {
            guard subdiagonal.count == diagonal.dropFirst().count else { throw .order }
            self.diagonal = diagonal
            self.subdiagonal = subdiagonal
        }
    }
}

extension Eigen.Tridiagonal: Equatable {}
extension Eigen.Tridiagonal: Hashable where Scalar: Hashable {}
extension Eigen.Tridiagonal: Sendable where Scalar: Sendable {}

extension Eigen.Tridiagonal {
    /// The eigenvalues, ascending.
    @inlinable
    public var values: [Scalar] {
        Self.decomposed(diagonal: diagonal, subdiagonal: subdiagonal, rows: 0).values
    }

    /// The eigenpairs, ascending by value, each vector of unit length; together the vectors are orthonormal.
    @inlinable
    public var pairs: [Eigen.Pair<Scalar>] {
        let order = diagonal.count
        let result = Self.decomposed(diagonal: diagonal, subdiagonal: subdiagonal, rows: order)
        return (0..<order).map { index in
            Eigen.Pair(value: result.values[index], vector: Array(result.vectors[(index * order)..<((index + 1) * order)]))
        }
    }

    /// Implicit QL with Wilkinson shifts (the EISPACK `tql2` iteration), ascending.
    ///
    /// `rows` is the number of trailing rows of the eigenvector matrix to carry: the order for every vector, 1 for the
    /// last component of each, which is what a Lanczos convergence test needs, or 0 for none. `vectors` holds them
    /// eigenvector by eigenvector.
    @inlinable
    internal static func decomposed(
        diagonal: [Scalar],
        subdiagonal: [Scalar],
        rows: Int
    ) -> (values: [Scalar], vectors: [Scalar]) {
        let order = diagonal.count
        var d = diagonal
        var e = subdiagonal + [0]
        if order == 0 { return ([], []) }
        var z = [Scalar](repeating: 0, count: order * rows)
        var row = 0
        while row < rows {
            // WHY: row `row` of the carried block is row `order − rows + row` of the identity.
            z[(order - rows + row) * rows + row] = 1
            row += 1
        }
        let epsilon = Scalar.ulpOfOne
        var l = 0
        while l < order {
            var iterations = 0
            while true {
                var m = l
                while m < order - 1 {
                    let size = Swift.abs(d[m]) + Swift.abs(d[m + 1])
                    if Swift.abs(e[m]) <= epsilon * size { break }
                    m += 1
                }
                if m == l || iterations == 60 { break }
                iterations += 1
                var g = (d[l + 1] - d[l]) / (2 * e[l])
                var r = Self.hypotenuse(g, 1)
                g = d[m] - d[l] + e[l] / (g + (g >= 0 ? r : -r))
                var s: Scalar = 1
                var c: Scalar = 1
                var p: Scalar = 0
                var i = m - 1
                var underflow = false
                while i >= l {
                    var f = s * e[i]
                    let b = c * e[i]
                    r = Self.hypotenuse(f, g)
                    e[i + 1] = r
                    if r == 0 {
                        d[i + 1] -= p
                        e[m] = 0
                        underflow = true
                        break
                    }
                    s = f / r
                    c = g / r
                    g = d[i + 1] - p
                    r = (d[i] - g) * s + 2 * c * b
                    p = s * r
                    d[i + 1] = g + p
                    g = c * r - b
                    var row = 0
                    while row < rows {
                        f = z[(i + 1) * rows + row]
                        z[(i + 1) * rows + row] = s * z[i * rows + row] + c * f
                        z[i * rows + row] = c * z[i * rows + row] - s * f
                        row += 1
                    }
                    i -= 1
                }
                if underflow { continue }
                d[l] -= p
                e[l] = g
                e[m] = 0
            }
            l += 1
        }
        let sorted = (0..<order).sorted { d[$0] < d[$1] }
        var vectors = [Scalar](repeating: 0, count: order * rows)
        for (index, original) in sorted.enumerated() {
            var row = 0
            while row < rows {
                vectors[index * rows + row] = z[original * rows + row]
                row += 1
            }
        }
        return (sorted.map { d[$0] }, vectors)
    }

    @inlinable
    internal static func hypotenuse(_ a: Scalar, _ b: Scalar) -> Scalar {
        let x = Swift.abs(a)
        let y = Swift.abs(b)
        let large = Swift.max(x, y)
        if large == 0 { return 0 }
        let ratio = Swift.min(x, y) / large
        return large * (1 + ratio * ratio).squareRoot()
    }
}
