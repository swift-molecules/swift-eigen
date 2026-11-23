import Eigen
import Eigen_Test_Support
import Foundation
import Sparse_Test_Support
import Testing

enum `Eigen Tests` {}

extension `Eigen Tests` {
    static func relative(_ a: [Double], _ b: [Double]) -> Double {
        precondition(a.count == b.count)
        return zip(a, b).reduce(0) { Swift.max($0, Swift.abs($1.0 - $1.1) / Swift.max(Swift.abs($1.1), 1e-300)) }
    }

    static func string(points: Int) -> [Double] {
        let h = 1 / Double(points + 1)
        return (1...points).map { p in
            let c = cos(Double(p) * Double.pi * h)
            return 6 / (h * h) * (1 - c) / (2 + c)
        }.sorted()
    }

    static func plate(rows: Int, columns: Int, spacing h: Double) -> [Double] {
        func mu(_ p: Int, _ n: Int) -> Double {
            let s = sin(Double(p) * Double.pi / Double(2 * (n + 1)))
            return 4 * s * s / (h * h)
        }
        var values: [Double] = []
        for p in 1...columns {
            for q in 1...rows {
                let sum = mu(p, columns) + mu(q, rows)
                values.append(sum * sum)
            }
        }
        return values.sorted()
    }

    /// The generalized eigenvalues of dense `a`, `b` by Cholesky reduction and cyclic Jacobi rotations.
    static func dense(_ a: [[Double]], _ b: [[Double]]) -> [Double] {
        let n = a.count
        var l = [[Double]](repeating: [Double](repeating: 0, count: n), count: n)
        for i in 0..<n {
            for j in 0...i {
                var sum = b[i][j]
                for k in 0..<j { sum -= l[i][k] * l[j][k] }
                l[i][j] = i == j ? sum.squareRoot() : sum / l[j][j]
            }
        }
        // c = L⁻¹ A L⁻ᵀ
        var y = a
        for column in 0..<n {
            for i in 0..<n {
                var sum = y[i][column]
                for k in 0..<i { sum -= l[i][k] * y[k][column] }
                y[i][column] = sum / l[i][i]
            }
        }
        var c = y
        for row in 0..<n {
            for i in 0..<n {
                var sum = y[row][i]
                for k in 0..<i { sum -= l[i][k] * c[row][k] }
                c[row][i] = sum / l[i][i]
            }
        }
        for _ in 0..<100 {
            var off = 0.0
            for i in 0..<n { for j in 0..<n where i != j { off += c[i][j] * c[i][j] } }
            if off < 1e-28 { break }
            for p in 0..<n {
                for q in (p + 1)..<n where c[p][q] != 0 {
                    let theta = (c[q][q] - c[p][p]) / (2 * c[p][q])
                    let t = (theta >= 0 ? 1 : -1) / (Swift.abs(theta) + (theta * theta + 1).squareRoot())
                    let cs = 1 / (t * t + 1).squareRoot()
                    let sn = t * cs
                    for k in 0..<n {
                        let ckp = c[k][p]
                        let ckq = c[k][q]
                        c[k][p] = cs * ckp - sn * ckq
                        c[k][q] = sn * ckp + cs * ckq
                    }
                    for k in 0..<n {
                        let cpk = c[p][k]
                        let cqk = c[q][k]
                        c[p][k] = cs * cpk - sn * cqk
                        c[q][k] = sn * cpk + cs * cqk
                    }
                }
            }
        }
        return (0..<n).map { c[$0][$0] }.sorted()
    }
}

extension `Eigen Tests` {
    @Suite
    struct Tridiagonal {
        @Test
        func `The eigenvalues of tridiag(−1, 2, −1) are 2 − 2cos(pπ/(n+1))`() throws {
            let n = 60
            let matrix = try Eigen.Tridiagonal(
                diagonal: [Double](repeating: 2, count: n),
                offDiagonal: [Double](repeating: -1, count: n - 1)
            )
            let expected = (1...n).map { 2 - 2 * cos(Double($0) * Double.pi / Double(n + 1)) }.sorted()
            #expect(zip(matrix.values, expected).allSatisfy { Swift.abs($0 - $1) < 1e-13 })
        }

        @Test
        func `The pairs are orthonormal eigenpairs`() throws {
            let diagonal = (0..<25).map { Double(($0 * 37) % 11) - 5 }
            let offDiagonal = (0..<24).map { Double(($0 * 13) % 7) - 3 }
            let matrix = try Eigen.Tridiagonal(diagonal: diagonal, offDiagonal: offDiagonal)
            let pairs = matrix.pairs
            for pair in pairs {
                for row in 0..<25 {
                    var product = diagonal[row] * pair.vector[row]
                    if row > 0 { product += offDiagonal[row - 1] * pair.vector[row - 1] }
                    if row < 24 { product += offDiagonal[row] * pair.vector[row + 1] }
                    #expect(Swift.abs(product - pair.value * pair.vector[row]) < 1e-12)
                }
            }
            for i in 0..<25 {
                for j in 0..<25 {
                    let dot = zip(pairs[i].vector, pairs[j].vector).reduce(0) { $0 + $1.0 * $1.1 }
                    #expect(Swift.abs(dot - (i == j ? 1 : 0)) < 1e-12)
                }
            }
            #expect(pairs.map(\.value) == matrix.values)
        }

        @Test
        func `A diagonal of the wrong length is refused`() {
            #expect(throws: Eigen.Error.order) { try Eigen.Tridiagonal<Double>(diagonal: [1, 2], offDiagonal: []) }
        }
    }

    @Suite
    struct String {
        @Test
        func `A fixed string's modes match the discrete closed form`() throws {
            let pencil = Eigen.Pencil<Double>.string(points: 400)
            let pairs = try pencil.lowest(30)
            #expect(`Eigen Tests`.relative(pairs.map(\.value), Array(`Eigen Tests`.string(points: 400).prefix(30))) < 1e-10)
            #expect(pencil.residual(of: pairs) < 1e-10)
            #expect(pencil.orthogonality(of: pairs) < 1e-10)
        }

        @Test
        func `A fixed string's lowest modes approach (pπ)²`() throws {
            let pairs = try Eigen.Pencil<Double>.string(points: 400).lowest(5)
            for (index, pair) in pairs.enumerated() {
                let p = Double(index + 1)
                // Linear elements with consistent mass overestimate by about (pπh)²/12.
                let h = 1.0 / 401
                #expect(Swift.abs(pair.value / (p * p * Double.pi * Double.pi) - 1) < pow(p * Double.pi * h, 2) / 10)
            }
        }

        @Test
        func `A free string has a rigid mode at zero and the rest in closed form`() throws {
            let points = 300
            let pencil = Eigen.Pencil<Double>.freeString(points: points)
            let pairs = try pencil.lowest(12)
            let h = 1 / Double(points - 1)
            let expected = (0..<12).map { p in
                let c = cos(Double(p) * Double.pi / Double(points - 1))
                return 6 / (h * h) * (1 - c) / (2 + c)
            }
            #expect(Swift.abs(pairs[0].value) < 1e-8)
            #expect(`Eigen Tests`.relative(pairs.dropFirst().map(\.value), Array(expected.dropFirst())) < 1e-10)
            #expect(pencil.residual(of: pairs) < 1e-10)
            #expect(pencil.orthogonality(of: pairs) < 1e-10)
        }

        @Test
        func `Many modes are found over many shifts and none is skipped`() throws {
            let pencil = Eigen.Pencil<Double>.string(points: 1500)
            let pairs = try pencil.lowest(400, options: .init(basis: 120))
            #expect(`Eigen Tests`.relative(pairs.map(\.value), Array(`Eigen Tests`.string(points: 1500).prefix(400))) < 1e-9)
            #expect(pencil.residual(of: pairs) < 1e-10)
            #expect(pencil.orthogonality(of: pairs) < 1e-9)
        }
    }

    @Suite
    struct Plate {
        @Test
        func `A square plate's modes, repeated ones included, match the discrete closed form`() throws {
            let side = 30
            let h = 1 / Double(side + 1)
            let pencil = Eigen.Pencil<Double>.plate(rows: side, columns: side, spacing: h)
            let pairs = try pencil.lowest(60)
            let expected = Array(`Eigen Tests`.plate(rows: side, columns: side, spacing: h).prefix(60))
            #expect(`Eigen Tests`.relative(pairs.map(\.value), expected) < 1e-9)
            #expect(pencil.residual(of: pairs) < 1e-10)
            #expect(pencil.orthogonality(of: pairs) < 1e-10)
        }

        @Test
        func `A rectangular plate's lowest modes approach π⁴((p/a)² + (q/b)²)²`() throws {
            // 1.5 by 1, h = 1/40.
            let h = 1.0 / 40
            let pencil = Eigen.Pencil<Double>.plate(rows: 39, columns: 59, spacing: h)
            let pairs = try pencil.lowest(4)
            let pi4 = pow(Double.pi, 4)
            let continuum = [(1, 1), (2, 1), (3, 1), (1, 2)].map { p, q in
                let s = pow(Double(p) / 1.5, 2) + pow(Double(q), 2)
                return pi4 * s * s
            }.sorted()
            #expect(`Eigen Tests`.relative(pairs.map(\.value), continuum) < 5e-3)
        }
    }

    @Suite
    struct Pencils {
        @Test(arguments: [1, 2, 3])
        func `Every eigenvalue of a scattered pencil matches a dense solve`(seed: UInt64) throws {
            let n = 40
            let a = try Sparse.Symmetric<Double>.scattered(order: n, perRow: 3, seed: seed)
                .adding(.tridiagonal(order: n, diagonal: 1, offDiagonal: 0), scaledBy: -3)
            let b = Sparse.Symmetric<Double>.scattered(order: n, perRow: 2, seed: seed + 10)
            let pencil = try Eigen.Pencil(a, b)
            let pairs = try pencil.lowest(n)
            let expected = `Eigen Tests`.dense(a.dense, b.dense)
            #expect(expected[0] < 0)
            #expect(zip(pairs.map(\.value), expected).allSatisfy { Swift.abs($0 - $1) < 1e-9 * Swift.max(1, Swift.abs($1)) })
            #expect(pencil.residual(of: pairs) < 1e-10)
            #expect(pencil.orthogonality(of: pairs) < 1e-9)
        }

        @Test
        func `Bad questions are refused`() throws {
            let pencil = Eigen.Pencil<Double>.string(points: 5)
            #expect(throws: Eigen.Error.count) { try pencil.lowest(6) }
            #expect(throws: Eigen.Error.count) { try pencil.lowest(-1) }
            #expect(try pencil.lowest(0).isEmpty)
            #expect(throws: Eigen.Error.order) {
                try Eigen.Pencil(Sparse.Symmetric<Double>.laplacian(rows: 2, columns: 2), .laplacian(rows: 1, columns: 2))
            }
            #expect(throws: Eigen.Error.shift) {
                try Eigen.Pencil(Sparse.Symmetric<Double>.laplacian(rows: 2, columns: 2), .laplacian(rows: 2, columns: 2).scaled(by: -1))
                    .lowest(1)
            }
        }
    }
}
