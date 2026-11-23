# swift-eigen

Use this package for the lowest eigenpairs of a large sparse symmetric pencil `A x = λ B x` with `B` positive definite: the squared angular frequencies and mode shapes of a structure from its stiffness and mass matrices. Start from `Eigen.Pencil<Scalar>`, built from two `Sparse.Symmetric` matrices of swift-sparse, and ask for `lowest(_:options:)`. `Eigen.Tridiagonal` holds a symmetric tridiagonal matrix and its eigenpairs. Counts are `Cardinal`s. `Eigen Test Support` holds a string, fixed or free, and a simply supported plate, each with a closed-form spectrum, and the residual and `B`-orthogonality laws as functions.

## Laws

- `lowest(k)` returns the `k` lowest eigenvalues, ascending, each counted as often as its multiplicity, with an eigenvector each. That none was skipped is checked, not assumed: the number of negative pivots of `A − σB` is the number of eigenvalues below `σ` (Sylvester's law of inertia), and the solver stops only when that count, at a shift above the `k`th value, equals the number of values it found below the shift. Otherwise it fails with `.incomplete`.
- Every returned pair meets `‖A x − λ B x‖₂ ≤ tolerance · (‖A‖ + |λ| ‖B‖) · ‖x‖₂`, with the matrix norms the largest absolute row sums. Each pair is tested against this before it is accepted. The default tolerance is ε^¾ of the scalar type, about 1.8 × 10⁻¹² for `Double`.
- The vectors are `B`-orthonormal, `xᵢᵀ B xⱼ = δᵢⱼ`, up to rounding (below 10⁻⁹ in every test here). Vectors of a repeated eigenvalue are an orthonormal basis of its eigenspace, not a particular one.
- A singular `A` is allowed, as with a free structure's rigid-body modes: the first shift lies below every eigenvalue, so `A − σB` stays nonsingular.
- The same pencil, options and seed give the same pairs on every machine: the start vectors are Threefry draws.
- `Eigen.Tridiagonal.values` and `pairs` are ascending, and the `pairs` vectors are orthonormal.

## Method

Lanczos runs on `(A − σB)⁻¹ B`, which is self-adjoint in the `B` inner product and whose largest eigenvalues `θ = 1/(λ − σ)` belong to the `λ` nearest `σ`. Each step costs one solve with the `L D Lᵀ` factor of `A − σB` and three products with `B`. It is reorthogonalized twice by classical Gram–Schmidt against the run's basis and against the eigenvectors already found near `σ`. Every start vector is purged of all eigenvectors found so far, and a pair that would repeat one is caught when it is accepted. A run stops when enough Ritz values have converged, or after `basis` steps.

After each run a check point goes midway between the highest value found and the next Ritz value, never inside a cluster that rounding cannot separate. When its inertia matches the values found below it, the check becomes the verified bound and the next shift. When the inertia counts more, copies of a repeated eigenvalue are missing, since one Lanczos vector sees only one direction of an eigenspace. The next run then starts from a fresh vector at the same shift, reusing its factorization, and the missing copies are now the nearest values it has not been purged of. Every shift reuses one symbolic elimination, so each costs one numeric factorization.

The test laws come from closed forms: a string of linear elements with consistent mass, fixed and free; a simply supported square plate with repeated modes and a rectangular one, as the squared five-point Laplacian; and random pencils against a dense solve.

## Speed

Measured once on an Apple-silicon Mac (macOS 27, Swift 6.4.2, release build, one thread), with the default options, in seconds. The plate is the squared Laplacian on a square grid with `B = I`; the string is free, with consistent mass. Fill counts the factor's off-diagonal entries; the solve is one forward and back substitution.

| Pencil | n | Fill | Factor | Solve | k = 10 | k = 100 | k = 1000 |
|---|---:|---:|---:|---:|---:|---:|---:|
| String, free | 100,000 | 199,936 | 0.003 | 0.002 | 1.1 (k = 20) | | |
| Plate 100 × 100 | 10,000 | 723,007 | 0.17 | 0.003 | 0.87 | 2.1 | 73 |
| Plate 200 × 200 | 40,000 | 3,589,248 | 1.2 | 0.028 | 3.9 | 24 | 619 |
| Plate 316 × 316 | 99,856 | 10,085,536 | 5.1 | 0.034 | 29 | 46 | 2,903 |

Every run met the residual law with the largest scaled residual below 1.1 × 10⁻¹², and its sampled vectors were `B`-orthonormal to 1.3 × 10⁻¹¹. Thousands of modes of a large two-dimensional problem take long because each Lanczos step is one memory-bound solve and the reorthogonalization grows with the basis; a supernodal factor and block Lanczos are the known next steps.
