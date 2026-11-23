extension Eigen.Lanczos {
    /// `xᵀ y` over four running sums.
    @inlinable
    static func dot(_ x: [Scalar], _ y: [Scalar]) -> Scalar {
        precondition(x.count == y.count)
        return x.withUnsafeBufferPointer { x in
            y.withUnsafeBufferPointer { y in
                unsafe dot(x.baseAddress!, y.baseAddress!, x.count)
            }
        }
    }

    @inlinable
    static func dot(_ x: UnsafePointer<Scalar>, _ y: UnsafePointer<Scalar>, _ n: Int) -> Scalar {
        var s0: Scalar = 0
        var s1: Scalar = 0
        var s2: Scalar = 0
        var s3: Scalar = 0
        var index = 0
        while index + 4 <= n {
            unsafe s0 += x[index] * y[index]
            unsafe s1 += x[index + 1] * y[index + 1]
            unsafe s2 += x[index + 2] * y[index + 2]
            unsafe s3 += x[index + 3] * y[index + 3]
            index += 4
        }
        while index < n {
            unsafe s0 += x[index] * y[index]
            index += 1
        }
        return (s0 + s1) + (s2 + s3)
    }

    /// `y += factor · x`.
    @inlinable
    static func add(_ y: inout [Scalar], _ factor: Scalar, _ x: [Scalar]) {
        precondition(x.count == y.count)
        if factor == 0 { return }
        y.withUnsafeMutableBufferPointer { y in
            x.withUnsafeBufferPointer { x in
                for index in 0..<y.count { unsafe y[index] += factor * x[index] }
            }
        }
    }

    @inlinable
    static func scale(_ x: inout [Scalar], by factor: Scalar) {
        x.withUnsafeMutableBufferPointer { x in
            for index in 0..<x.count { unsafe x[index] *= factor }
        }
    }

    /// `Vᵀ q` for the first `columns` columns of the column-major `basis`, whose columns have `q.count` rows.
    @inlinable
    static func products(_ basis: [Scalar], columns: Int, _ q: [Scalar]) -> [Scalar] {
        let n = q.count
        var result = [Scalar](repeating: 0, count: columns)
        basis.withUnsafeBufferPointer { basis in
            q.withUnsafeBufferPointer { q in
                let v = basis.baseAddress!
                let q = q.baseAddress!
                var column = 0
                while column + 4 <= columns {
                    let a = unsafe v + column * n
                    let b = unsafe a + n
                    let c = unsafe b + n
                    let d = unsafe c + n
                    var sa: Scalar = 0
                    var sb: Scalar = 0
                    var sc: Scalar = 0
                    var sd: Scalar = 0
                    for row in 0..<n {
                        let x = unsafe q[row]
                        unsafe sa += a[row] * x
                        unsafe sb += b[row] * x
                        unsafe sc += c[row] * x
                        unsafe sd += d[row] * x
                    }
                    result[column] = sa
                    result[column + 1] = sb
                    result[column + 2] = sc
                    result[column + 3] = sd
                    column += 4
                }
                while column < columns {
                    result[column] = unsafe dot(v + column * n, q, n)
                    column += 1
                }
            }
        }
        return result
    }

    /// `w += V c` over the first `c.count` columns of the column-major `basis`, whose columns have `w.count` rows.
    @inlinable
    static func accumulate(_ w: inout [Scalar], _ basis: [Scalar], _ c: [Scalar]) {
        let n = w.count
        let columns = c.count
        basis.withUnsafeBufferPointer { basis in
            w.withUnsafeMutableBufferPointer { w in
                let v = basis.baseAddress!
                let w = w.baseAddress!
                var column = 0
                while column + 4 <= columns {
                    let a = unsafe v + column * n
                    let b = unsafe a + n
                    let cc = unsafe b + n
                    let d = unsafe cc + n
                    let ca = c[column]
                    let cb = c[column + 1]
                    let ccc = c[column + 2]
                    let cd = c[column + 3]
                    for row in 0..<n {
                        unsafe w[row] += ca * a[row] + cb * b[row] + ccc * cc[row] + cd * d[row]
                    }
                    column += 4
                }
                while column < columns {
                    let a = unsafe v + column * n
                    let ca = c[column]
                    for row in 0..<n { unsafe w[row] += ca * a[row] }
                    column += 1
                }
            }
        }
    }
}
