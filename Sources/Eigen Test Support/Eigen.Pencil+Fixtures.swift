public import Eigen
import Sparse_Test_Support

extension Eigen.Pencil {
    /// A string of unit length, tension and density, fixed at both ends, in `points` interior nodes of linear elements
    /// with the consistent mass matrix: `A = (1/h) tridiag(−1, 2, −1)`, `B = (h/6) tridiag(1, 4, 1)`, `h = 1/(points+1)`.
    ///
    /// Its eigenvalues are `(6/h²)(1 − c)/(2 + c)` with `c = cos(pπh)` for `p` in `1...points`; the continuum's are `(pπ)²`.
    public static func string(points: Int) -> Self {
        let h = 1 / Scalar(points + 1)
        return try! Self(
            .tridiagonal(order: points, diagonal: 2 / h, offDiagonal: -1 / h),
            .tridiagonal(order: points, diagonal: 4 * h / 6, offDiagonal: h / 6)
        )
    }

    /// The same string free at both ends, in `points` nodes including the ends, `h = 1/(points − 1)`: its lowest
    /// eigenvalue is zero, with the rigid translation as its mode, so `A` is singular.
    ///
    /// Its eigenvalues are `(6/h²)(1 − c)/(2 + c)` with `c = cos(pπ/(points − 1))` for `p` in `0..<points`.
    public static func freeString(points: Int) -> Self {
        let h = 1 / Scalar(points - 1)
        var stiffness: [Sparse.Symmetric<Scalar>.Entry] = []
        var mass: [Sparse.Symmetric<Scalar>.Entry] = []
        for element in 0..<(points - 1) {
            let i = element
            let j = element + 1
            stiffness += [.init(i, i, 1 / h), .init(j, j, 1 / h), .init(i, j, -1 / h)]
            mass += [.init(i, i, 2 * h / 6), .init(j, j, 2 * h / 6), .init(i, j, h / 6)]
        }
        return try! Self(
            Sparse.Symmetric(order: points, entries: stiffness),
            Sparse.Symmetric(order: points, entries: mass)
        )
    }

    /// A simply supported rectangular plate of unit bending stiffness and mass per area, `width` by `height`, on a
    /// `rows × columns` grid of interior points: `A` is the squared five-point Laplacian over `h⁴` and `B` the identity,
    /// with `h = width/(columns+1) = height/(rows+1)`.
    ///
    /// Its eigenvalues are `((4/h²)(sin²(pπ/(2(columns+1))) + sin²(qπ/(2(rows+1)))))²`; the continuum's are
    /// `π⁴((p/width)² + (q/height)²)²`.
    public static func plate(rows: Int, columns: Int, spacing h: Scalar) -> Self {
        let order = rows * columns
        return try! Self(
            Sparse.Symmetric.biharmonic(rows: rows, columns: columns).scaled(by: 1 / (h * h * h * h)),
            .tridiagonal(order: order, diagonal: 1, offDiagonal: 0)
        )
    }

    /// The largest relative residual of `pairs`: `‖A x − λ B x‖ / ((‖A‖ + |λ| ‖B‖) ‖x‖)` in the 2-norm.
    @inlinable
    public func residual(of pairs: [Eigen.Pair<Scalar>]) -> Scalar {
        let normA = a.norm
        let normB = b.norm
        var worst: Scalar = 0
        for pair in pairs {
            let ax = a.applied(to: pair.vector)
            let bx = b.applied(to: pair.vector)
            var residual: Scalar = 0
            var length: Scalar = 0
            for row in ax.indices {
                let r = ax[row] - pair.value * bx[row]
                residual += r * r
                length += pair.vector[row] * pair.vector[row]
            }
            worst = Swift.max(worst, residual.squareRoot() / ((normA + Swift.abs(pair.value) * normB) * length.squareRoot()))
        }
        return worst
    }

    /// The largest entry of `Xᵀ B X − I` for the vectors of `pairs`: how far they are from `B`-orthonormal.
    @inlinable
    public func orthogonality(of pairs: [Eigen.Pair<Scalar>]) -> Scalar {
        let products = pairs.map { b.applied(to: $0.vector) }
        var worst: Scalar = 0
        for i in pairs.indices {
            for j in i..<pairs.count {
                var sum: Scalar = 0
                for row in products[j].indices { sum += pairs[i].vector[row] * products[j][row] }
                worst = Swift.max(worst, Swift.abs(sum - (i == j ? 1 : 0)))
            }
        }
        return worst
    }
}
