public import Eigen
import Sparse_Test_Support

extension Eigen.Pencil {
    /// A string of unit length, tension and density in `points` nodes of linear elements with the consistent mass
    /// matrix, `A = (1/h) tridiag(−1, 2, −1)` and `B = (h/6) tridiag(1, 4, 1)` away from the ends.
    ///
    /// Fixed at both ends, the nodes are interior, `h = 1/(points + 1)`, and the eigenvalues are
    /// `(6/h²)(1 − c)/(2 + c)` with `c = cos(pπh)` for `p` in `1...points`; the continuum's are `(pπ)²`.
    /// Free at both ends, the nodes include the ends, `h = 1/(points − 1)`, `A` is singular with the rigid translation
    /// as the mode of eigenvalue zero, and the eigenvalues are the same expression with `c = cos(pπh)` for `p` in
    /// `0..<points`.
    public static func string(points: Cardinal, fixed: Bool = true) throws(Eigen.Error) -> Self {
        let points = Int(bitPattern: points)
        let h = fixed ? 1 / Scalar(points + 1) : 1 / Scalar(points - 1)
        var stiffness: [Sparse.Symmetric<Scalar>.Entry] = []
        var mass: [Sparse.Symmetric<Scalar>.Entry] = []
        let elements = fixed ? points + 1 : points - 1
        var element = 0
        while element < elements {
            let i = fixed ? element - 1 : element
            let j = i + 1
            if i >= 0 {
                stiffness.append(.init(i, i, 1 / h))
                mass.append(.init(i, i, 2 * h / 6))
            }
            if j < points {
                stiffness.append(.init(j, j, 1 / h))
                mass.append(.init(j, j, 2 * h / 6))
            }
            if i >= 0 && j < points {
                stiffness.append(.init(i, j, -1 / h))
                mass.append(.init(i, j, h / 6))
            }
            element += 1
        }
        let a: Sparse.Symmetric<Scalar>
        let b: Sparse.Symmetric<Scalar>
        do throws(Sparse.Error) {
            a = try Sparse.Symmetric(order: points, entries: stiffness)
            b = try Sparse.Symmetric(order: points, entries: mass)
        } catch {
            throw .order
        }
        return try Self(a, b)
    }

    /// A simply supported rectangular plate of unit bending stiffness and mass per area on a `rows × columns` grid of
    /// interior points `h` apart: `A` is the squared five-point Laplacian over `h⁴` and `B` the identity.
    ///
    /// Its eigenvalues are `((4/h²)(sin²(pπ/(2(columns+1))) + sin²(qπ/(2(rows+1)))))²`; the continuum's, for a plate
    /// `(columns + 1)h` by `(rows + 1)h`, are `π⁴((p/width)² + (q/height)²)²`.
    public static func plate(rows: Cardinal, columns: Cardinal, spacing h: Scalar) throws(Eigen.Error) -> Self {
        let a: Sparse.Symmetric<Scalar>
        let b: Sparse.Symmetric<Scalar>
        do throws(Sparse.Error) {
            let biharmonic = try Sparse.Symmetric<Scalar>.biharmonic(rows: rows, columns: columns)
            a = biharmonic.scaled(by: 1 / (h * h * h * h))
            b = try .tridiagonal(order: Cardinal(UInt(biharmonic.order)), diagonal: 1, offDiagonal: 0)
        } catch {
            throw .order
        }
        return try Self(a, b)
    }

    /// The largest relative residual of `pairs`: `‖A x − λ B x‖ / ((‖A‖ + |λ| ‖B‖) ‖x‖)` in the 2-norm.
    @inlinable
    public func residual(of pairs: [Eigen.Pair<Scalar>]) -> Scalar {
        let normA = a.norm
        let normB = b.norm
        return pairs.reduce(0) { worst, pair in
            let ax = a.applied(to: pair.vector)
            let bx = b.applied(to: pair.vector)
            var residual: Scalar = 0
            var length: Scalar = 0
            ax.indices.forEach { row in
                let r = ax[row] - pair.value * bx[row]
                residual += r * r
                length += pair.vector[row] * pair.vector[row]
            }
            return Swift.max(worst, residual.squareRoot() / ((normA + Swift.abs(pair.value) * normB) * length.squareRoot()))
        }
    }

    /// The largest entry of `Xᵀ B X − I` for the vectors of `pairs`: how far they are from `B`-orthonormal.
    @inlinable
    public func orthogonality(of pairs: [Eigen.Pair<Scalar>]) -> Scalar {
        let products = pairs.map { b.applied(to: $0.vector) }
        var worst: Scalar = 0
        pairs.indices.forEach { i in
            pairs.indices.filter { $0 >= i }.forEach { j in
                var sum: Scalar = 0
                products[j].indices.forEach { row in sum += pairs[i].vector[row] * products[j][row] }
                worst = Swift.max(worst, Swift.abs(sum - (i == j ? 1 : 0)))
            }
        }
        return worst
    }
}
