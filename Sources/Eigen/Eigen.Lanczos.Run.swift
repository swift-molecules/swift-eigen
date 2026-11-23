extension Eigen.Lanczos {
    /// What one run at one shift found.
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
}
