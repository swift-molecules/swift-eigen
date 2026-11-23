extension Eigen {
    /// An eigenvalue and an eigenvector that belongs to it.
    public struct Pair<Scalar> {
        /// The eigenvalue.
        public let value: Scalar
        /// The eigenvector.
        public let vector: [Scalar]

        @inlinable
        public init(value: Scalar, vector: [Scalar]) {
            self.value = value
            self.vector = vector
        }
    }
}

extension Eigen.Pair: Equatable where Scalar: Equatable {}
extension Eigen.Pair: Hashable where Scalar: Hashable {}
extension Eigen.Pair: Sendable where Scalar: Sendable {}
