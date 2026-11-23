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
        public var basis: Cardinal?
        /// The seed of the start vectors.
        public var seed: UInt64

        @inlinable
        public init(
            shift: Scalar? = nil,
            tolerance: Scalar = Scalar.ulpOfOne.squareRoot() * Scalar.ulpOfOne.squareRoot().squareRoot(),
            basis: Cardinal? = nil,
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
