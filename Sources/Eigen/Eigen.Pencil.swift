extension Eigen {
    /// The generalized symmetric eigenproblem `A x = λ B x`, with `A` and `B` sparse and symmetric and `B` positive
    /// definite: a stiffness and a mass matrix, whose eigenpairs are a structure's squared angular frequencies and modes.
    public struct Pencil<Scalar: BinaryFloatingPoint> {
        /// `A`, the stiffness.
        public let a: Sparse.Symmetric<Scalar>
        /// `B`, the mass: positive definite.
        public let b: Sparse.Symmetric<Scalar>

        /// Fails with `.order` when the two matrices have different orders.
        @inlinable
        public init(_ a: Sparse.Symmetric<Scalar>, _ b: Sparse.Symmetric<Scalar>) throws(Eigen.Error) {
            guard a.order == b.order else { throw .order }
            self.a = a
            self.b = b
        }
    }
}

extension Eigen.Pencil: Equatable {}
extension Eigen.Pencil: Sendable where Scalar: Sendable {}

extension Eigen.Pencil {
    /// How the lowest eigenpairs are searched for.
    public struct Options {
        /// The first shift, which must lie below the lowest eigenvalue; nil picks one from the diagonals and lowers it
        /// until the factorization of `A − σB` has no negative pivot.
        public var shift: Scalar?
        /// The bound on each returned pair's relative residual,
        /// `‖A x − λ B x‖ ≤ tolerance · (‖A‖ + |λ| ‖B‖) · ‖x‖` in the 2-norm, with the matrix norms the largest row sums.
        public var tolerance: Scalar
        /// The most Lanczos vectors one shift keeps; nil sizes it from the number still wanted, at most 300.
        public var basis: Int?
        /// The seed of the start vectors.
        public var seed: UInt64

        @inlinable
        public init(
            shift: Scalar? = nil,
            tolerance: Scalar = Scalar.ulpOfOne.squareRoot() * Scalar.ulpOfOne.squareRoot().squareRoot(),
            basis: Int? = nil,
            seed: UInt64 = 0
        ) {
            self.shift = shift
            self.tolerance = tolerance
            self.basis = basis
            self.seed = seed
        }
    }
}

extension Eigen.Pencil.Options: Equatable {}
extension Eigen.Pencil.Options: Sendable where Scalar: Sendable {}

extension Eigen.Pencil {
    /// The `count` lowest eigenpairs, ascending, each eigenvalue counted as often as its multiplicity.
    ///
    /// Shift-invert Lanczos with full reorthogonalization in the `B` inner product, sliced over shifts: each shift's
    /// run finds the eigenvalues nearest it, the next shift is placed just above the ones found, and the number of
    /// negative pivots of `A − σB` (Sylvester's law of inertia) checks that none below it was missed. The vectors are
    /// `B`-orthonormal: `xᵢᵀ B xⱼ` is one for `i == j` and zero otherwise, up to rounding.
    ///
    /// Fails with `.count` for a count below zero or above the order, with `.shift` when no shift below the spectrum is
    /// found, and with `.incomplete` when the slices do not verify `count` eigenvalues.
    @inlinable
    public func lowest(_ count: Int, options: Options = Options()) throws(Eigen.Error) -> [Eigen.Pair<Scalar>] {
        guard count >= 0, count <= a.order else { throw .count }
        if count == 0 { return [] }
        var solver = try Eigen.Lanczos(pencil: self, options: options)
        return try solver.lowest(count)
    }
}
