extension Eigen {
    /// Why an eigenproblem could not be solved.
    public enum Error: Swift.Error, Hashable, Sendable {
        /// The two matrices of a pencil have different orders.
        case order
        /// More eigenpairs were asked for than the problem has, or fewer than none.
        case count
        /// No shift below the spectrum was found, or the second matrix of the pencil is not positive definite.
        case shift
        /// The eigenpairs found could not be shown to be the lowest ones within the solver's budget of shifts.
        case incomplete
    }
}
