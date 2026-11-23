public import Eigen

extension Eigen.Pencil {
    /// How a string fixture is held at its two ends.
    public enum Ends: Hashable, Sendable {
        /// Held still: the ends are not nodes.
        case fixed
        /// Free to move: the ends are nodes, and the string has a rigid mode.
        case free
    }
}
