# Major corrections and changes to the DAFT book

This log accompanies the modernization of *The Discrete Algebra of the Fourier
Transform*, Gabriel Peyré's adapted English translation of *L'algèbre discrète
de la transformée de Fourier* (Ellipses, 2004, ISBN 9782729818678).

Work is in progress. Entries describe changes already made to the revised
sources; the final PDF, figure comparison, and execution checks are recorded
separately when they have been completed. This is a record of substantial
corrections, not an exhaustive list of punctuation and spelling edits.

## Preservation and attribution

- Preserved all 856 files of the previous edition, including its 296-page PDF,
  LaTeX sources, figures, MATLAB programs, and Maple worksheets, in
  `legacy/pre-python-2026-09-08/`. Verified the copy against SHA-256 hashes;
  `legacy/manifest.json` records the original paths, sizes, and hashes.
- Added explicit attribution to the French original on the revised PDF cover,
  including publisher, year, ISBN, and a link to the publisher's page.

## Editorial and production changes

- Established shared typography matching *Mathematical Foundations of Data
  Sciences*: Palatino-style text and mathematics, Source Sans Pro headings,
  restrained teal accents, improved margins, and compact Python listings.
  Final compilation and visual checks are still pending.

## Chapter 1 — Fourier analysis on finite groups

- Rewrote the exposition and all eleven exercises in standard mathematical
  English, preserving every existing cross-reference label. Added a short
  chapter overview and clarified what the reader gains from the finite setting.
- **Dirac basis normalization:** with the normalized inner product, the Dirac
  vectors have squared norm `1/|G|`; they are orthogonal, not orthonormal.
  Corrected the proof and stated the required `sqrt(|G|)` rescaling.
- **Cyclic orthogonality:** replaced the erroneous exponent `pq` by `p−q` in
  the geometric sum proving orthogonality of two characters.
- **Character extension:** rebuilt the proof. If `n` is the order of `xH`, the
  extension requires `ζ^n = χ(x^n)`, not `ζ^n = 1`. Proved uniqueness of the
  coset representation, checked multiplication across a wrap of the exponent,
  and justified induction on the subgroup index. Added the instructive
  `Z/4Z` counterexample to the previous root choice.
- **Duality:** distinguished canonical pullback and bidual maps from a chosen
  isomorphism with the dual. Required invariant factors to be at least two
  for uniqueness, and corrected the misleading discussion of infinite-group
  biduality. Corrected the quotient in the abelianization diagram to `G/D(G)`.
- **Noncommutative groups:** distinguished linear characters from characters
  of higher-dimensional representations; stated the restriction `n ≥ 2` for
  the sign character of `S_n`; included fixed points in cycle counting; and
  repaired the orientation-dependent statement about products of transpositions.
- **Convolution:** proved associativity directly on the Dirac basis, specified
  the identity, and required algebra homomorphisms to preserve the unit.
  Made the additive versus multiplicative probability notation explicit.
- **Fourier conventions:** explicitly related the positive-character convention
  in the group-theoretic chapters to the negative-exponential numerical DFT.
  Gave the normalized Plancherel identity and the unitary scaling.
- **Exercises:** required complex linearity in the Heisenberg commutant
  statement and allowed the scalar to be zero; defined an orthonormal orbit
  as an indexed family; justified positive square roots via spectral
  projectors; and included translation invariance in the indicator exercise.
- **Random walks:** reconciled conditional probabilities with column-vector
  evolution (`P[j,i] = p[i,j]`), corrected the asymmetric increment law, and
  stated the spectral condition for convergence to the uniform distribution.
- **Uncertainty principle:** required an abelian group, replaced the supremum
  of complex values by their maximum modulus, and restored modulation to the
  description of all equality cases.
- **Figures 1.1–1.4:** regenerated all four figures with Python, vector output,
  boxed axes, consistent typography, and explicit mathematical labels.
  Corrected the centering of the bell function and its Fourier transform.
  Corrected the circular self-convolution: for the displayed box supported on
  `{5,…,10}` in `Z/16Z`, the maximum is at `15`, not `7`; the new plot shows
  the required wrap-around. Added the uniform law to the random-walk panels.
- **Checks completed:** five independent computational checks covering
  character extension for cyclic groups through order 24, Fourier signs and
  normalization, convolution, the Heisenberg intertwining identity, rotation
  traces, indicator bounds, and uncertainty equality examples.
  All four figures were generated and visually checked. A provisional chapter
  build (including its solutions) succeeds without overfull or underfull boxes.
  Cross-references to later chapters remain pending until the global build.

### Chapter 1 solutions and programs

- Rewrote all eleven solutions. Restored the missing `1/|G|` in the counting
  argument, corrected the circulant entry to `f(i−j)`, and matched the first-row
  convention of the displayed matrix.
- Replaced the false claim that translations are group automorphisms by the
  correct unit-modulus Fourier factor.
- Replaced the unsupported appeal to Schur's lemma by a direct commutant proof
  using modulations and translations. Corrected the claim that the Heisenberg
  group is finite and identified the Fourier transform itself as the canonical
  intertwining map for the specified actions.
- Gave a rigorous periodization proof for orthogonal integer translates, with
  the `1/(2π)` factor and the bound `A^(−1/2)`. Restored missing factors in the
  finite-group projector calculation and corrected the normalized energy to
  `1/|G|`, rather than `1`.
- Restored the square in the pointwise probability bound, the frequency index
  in the random-walk eigenvalues, and the distinction between the group size
  and iteration number in the parity argument.
- Corrected the uncertainty proof's support set and normalization; gave a
  cyclic-gap argument that also handles nondivisible support sizes.
- Converted Programs `listing-calcul-det-circ`, `listing-calcul-trace-beta`,
  and `listing-calcul-proba-pk` to NumPy/SymPy. All three were executed and
  compared with independent matrix or symbolic references. Eight tests now
  pass in total (five formula checks and three program checks).
- Fixed the mistyped bibliography key for Matusiak–Özaydin–Przebinda, which
  previously left the uncertainty-principle citation unresolved.

## Chapter 2 — Gauss sums (completed section review)

- Corrected the opening quadratic-congruence example: division requires the
  denominator to be nonzero in the field, and modular ratios cannot be
  interpreted as unrestricted integer equalities for `x/y`.
- Corrected the domains, codomains, and index ranges of additive and
  multiplicative characters, including the missing `a−b` in the proof that
  additive characters are distinct. Required a prime power for `F_q`.
- Repaired the Gauss-sum norm proof: zero extension gives squared normalized
  norm `(q−1)/q`, so Plancherel gives `q−1`, not `q`, for the squared norm of
  the transformed function. Separating its zero coefficient yields the correct
  magnitude `sqrt(q)`. Required a nonzero multiplier in the change of variable.
- Reworked the determinant argument for the sign of the quadratic Gauss sum,
  correcting the two-by-two block orientation, the summation endpoint, missing
  differences in the Vandermonde product, and the count of phase factors.
- **Extension-field sign:** the phase depends on the characteristic `p`, not
  merely on `q mod 4`. The corrected Hasse–Davenport expression is
  `(-1)^(s−1) i^(s(p−1)^2/4) sqrt(q)` for `q=p^s`. In particular, the
  canonical quadratic sum over `F_9` is `+3`, whereas the old formula gave `−3`.
- Refined the quadratic-reciprocity proof by specifying the coefficient ring,
  justifying reduction back to integer congruences, and explaining why the
  Gauss sum can be canceled modulo the other prime. Corrected confused
  additive/multiplicative character symbols. Clarified that this proof needs
  the square of the Gauss sum, not its separately evaluated sign.
- Three additional independent checks pass: character and Gauss identities in
  nine prime and extension fields, determinant phases for seven odd primes,
  and quadratic reciprocity for all distinct pairs in a set of fourteen primes.



### Walsh transforms, Poisson summation, and Chapter 2 exercises

- Rewrote all remaining theory, all ten exercises, and all ten corresponding
  solutions; retained all existing labels. Clarified the distinction between
  Hamming degree on the Boolean cube and sequency along a sampled signal.
- **Walsh convolution:** replaced the false assertion that Walsh transforms
  have no convolution property with the XOR-convolution theorem, distinguishing
  it from cyclic convolution. Corrected the inverse sum to contain all `2^k`
  coefficients and the explanation of the fast algorithm's complexity.
- **Poisson summation:** corrected the spaces containing linear annihilators,
  the quotient-space basis in the trace proof, and the convolution eigenvalue's
  inverse character. Gave the continuity, interchange, and Fourier uniqueness
  arguments needed for the continuous formula.
- **MacWilliams identities:** restored the missing subtraction in the weight
  enumerator's exponents and corrected the binary substitution to `(X+Y,X−Y)`.
  For `F_q`, the substitution is `(X+(q−1)Y,X−Y)`, and the code must be
  `F_q`-linear. An additive subgroup need not have that property. Added the
  `F_2` subgroup of `F_4` as a counterexample and proved the annihilator
  identification using scalar closure and the nondegenerate trace pairing.
  Removed the false claim that an extension-field additive character is injective.
- **Geometric reciprocity proof:** replaced an invalid reflection across an
  oblique diagonal, and a reference to an undefined point, with an exact pairing
  of columns `a` and `p−a`. Specified which boundary points are excluded and
  counted the half-rectangle rigorously. The new diagram uses `r` consistently.
- **Fermat equation:** applied the additive counting estimate to `−H_k`,
  justified the number of preimages of a power, and supplied an explicit
  positive polynomial difference proving the sufficient size condition.
  Treated exponents not dividing `q−1` through `gcd(k,q−1)`.
- **Factorial contrasts:** corrected malformed sign indices, distinguished
  an averaged main effect from fixing other factors, and computed all seven
  effects and the mean. Specified the binary encoding and the different
  normalizations of the mean and nonconstant contrasts.
- **Haar indexing and normalization:** corrected the wavelet range to
  `j≥0, 0≤k<2^j` and the index to `n=2^j+k`. Restored the scaling term in
  projections, repaired the detail-space and refinement indices, and changed
  the detail sign to `(even−odd)/sqrt(2)` to match the defined mother wavelet.
- **Haar analysis:** corrected the initial cell-value factor to `2^(−J/2)`,
  the intermediate vector lengths to `2^(J−i)`, and the final mean coefficient
  to `sqrt(N)` times the input mean. The elementary block is a reflection,
  not a rotation. Corrected the packed vector, removed the repeated coarse
  energy term, and identified the analysis basis as the rows of the transform.
- **Haar convergence:** replaced an unfinished argument by a proof for every
  partial sum using averages on an adaptive dyadic partition. Specified endpoint
  representatives and derived completeness in `L^2` from density.
- **Walsh approximation:** made every retained count exact, corrected off-by-one
  indexing, and specified a complete two-dimensional ordering including ties.
  Distinguished orthogonal approximation from finite-bit encoding. Added a
  tested image-approximation routine alongside the separable transform.
- **Paley construction:** corrected `Q[i,j]` from `η(ji)` to `η(j−i)`, the
  diagonal correlation to a sum of squares, a repeated block count, and block
  transpose dimensions. Replaced the misleading determinant-minimization
  discussion with the correct maximization and a proof that a vertex maximizer
  exists in the real entrywise cube.
- **Tensor transforms:** corrected the tensor matrix size to `st×st`, the
  inverse formula, and the operation count. Stated orthogonal normalization
  and distinguished bitwise complementation from bit-order reversal.
- **Sampling:** repaired the Dirac-comb scaling and the missing `s` in its
  proof. Rebuilt the sampling argument from a convergent Fourier series,
  corrected the sinc shift to `t−nT`, its inverse-transform factor, and its
  squared norm `T`. Proved the sample energy identity, `L^2` and uniform
  convergence, and the correct projection factor `1/T`. Statements about
  arbitrary `L^2` Fourier transforms now hold almost everywhere, not at
  undefined band endpoints.
- **Theta evaluation:** separated a rigorous truncation-error bound from
  floating-point rounding; removed a claim of accuracy beyond the arithmetic
  precision being used.
- **Python conversion:** ten of the original 88 numbered programs now read
  directly from tested Python functions. They include circulant determinants,
  symbolic rotation traces, random walks, Walsh and Haar analysis, two-dimensional
  Walsh transforms, sequency ordering, Paley matrices, and tensor powers.
- **Checks completed:** 22 tests now cover the corrected formulas and programs.
  New checks include exact MacWilliams polynomial identities over prime and
  extension fields, lattice pairing, all factorial effects, image projections
  against explicit basis sums, theta identities, and the truncation bound.

- **Chapter 2 figures:** regenerated all twelve with Python and TikZ, including
  the Boolean cube, reciprocity geometry, Haar cascade and scale-separated
  coefficients, Walsh orderings, and exact coefficient-count approximations.
  The provisional two-chapter build has no box warnings. The 17-page working
  comparator now contains sixteen figure pairs and compiles without warnings.

## Chapter 3 — the discrete Fourier transform

- Rewrote the complete chapter, all eleven exercises, and all eleven solutions,
  preserving the labels and connecting every algorithm to an explicit storage,
  sign, and normalization convention.
- **DFT formulas:** corrected the sample spacing to `(b−a)/N`, the Plancherel
  endpoints to `N−1`, and the inverse-transform index to `N−n`. Reconciled
  the negative-exponential DFT with the character indexing of Chapter 1.
- **FFT recurrences:** restored the missing subtraction in the second output
  half, corrected the DIF phase to depend on the input summation index, and
  made inverse-sign changes apply to both subtransforms and twiddle factors.
  Repaired the cost recurrence and the input-matrix orientation in Cooley–Tukey.
- **Bit reversal:** replaced a confused change-of-input argument by the actual
  recursion-tree permutation. Distinguished cyclic reversal, bitwise
  complementation, and bit-order reversal. Corrected the recursive copying
  cost and removed a false claim that array-based Python uses no temporary memory.
- **Convolution:** restored `n−k`, `i−k`, and `j−l` wherever differences had
  become products. Gave the precise padding length and negative-index placement.
  Distinguished support interval length from sparse support cardinality and
  explained the different boundary models without claiming zero padding removes
  all boundary artifacts.
- **Zero padding:** clarified that padding input data changes the sampled
  frequencies; it is not the same length-N DFT. Padding for linear convolution
  is exact when the no-wrap condition holds.
- **Matrices and images:** corrected interpolation degree to at most `N−1`,
  distinguished symmetry from Hermitian symmetry, retained the permutation in
  the FFT matrix factorization, and gave correctly scaled frequencies for
  rectangular images. Distinguished a spectrum's magnitude from its phase.
- **Spectral decomposition:** corrected `F^4` to `N^2 I`, restored the missing
  factor `1/(2 sqrt(N))` in the sum of the four spectral components, and noted
  that a zero component is not an eigenvector. Gave explicit orthogonal projectors.
- **Fractional transforms:** fixed the missing normalization in the claimed
  square root and used unitary powers for energy-preserving comparisons. Fixed
  logarithm branches consistently across all real exponents to justify the
  composition law. Explained why changing a basis within an eigenspace alone
  does not change a fixed scalar branch of a matrix power.
- **DFT eigenbasis:** the original proof falsely asserted that the commuting
  matrix S always has distinct eigenvalues. At N=4, eigenvalue −4 has two
  independent eigenvectors. Rebuilt the construction by diagonalizing separate
  symmetric and antisymmetric blocks. Each block is irreducible tridiagonal;
  the relevant nonzero entries are on the subdiagonal, not the diagonal.
  Gave a complete simple-spectrum proof for each block and deterministic
  ordering/sign conventions, including ties between blocks.
- **Good–Thomas and split radix:** corrected the Chinese-remainder output
  indexing and removed a false bound before reduction modulo N. Repaired
  split-radix signs and summation endpoints; supplied a complete DIT recursion
  and explicit arithmetic-count recurrences instead of an unsupported universal
  optimality claim.
- **Interpolation and differentiation:** preserved complex inputs in spectral
  interpolation, specified its scaling and odd-length assumption, restored the
  factor two in the Chebyshev recurrence proof, and distinguished the book's
  cosine-transform normalization from SciPy's. Corrected the derivative
  multiplier to `(i ξ)^m`, specified a fractional branch and zero mode, and
  included the sample period and a consistent treatment of real data.
- **Cyclic orthogonalization:** restored the missing modulation/translation
  sign and the unitary scaling. For the book's normalized inner product, the
  required Fourier magnitude is `sqrt(n)`, not one. Repaired correlation
  conjugates, the `1/n` coefficient factor, and reconstruction.
- **Python checks:** 27 of the 88 numbered programs are now converted and
  executed. Twelve additional tests compare FFT variants against explicit
  Fourier matrices, overlap-add against direct convolution, cosine transforms
  against their defining sums and SciPy, and interpolation against polynomial
  and trigonometric references. They also check fractional composition laws,
  integer derivatives, all parity eigenbases for sizes 3–40, and normalized
  translation Gram matrices. A total of 34 tests pass.
- **Figures and layout:** regenerated all eighteen Chapter 3 figures and
  integrated them into a successful three-chapter, 72-page provisional build.
  Selected pages and every generated figure were visually reviewed. The
  comparator now contains 34 pairs and compiles without warnings.
- **Figure content:** replaced the absolute real part mislabeled as a matrix
  magnitude by the actual complex magnitude; showed real and imaginary parts
  of consistently normalized fractional transforms; used Chebyshev roots,
  matching the exercise, instead of the old generator's Lobatto nodes.
  Corrected the eight-point butterfly wiring and labeled bit-reversal swaps.
  Replaced a linear surface convolution mislabeled as circular by a true
  circular computation. The regenerated figure code checks the relevant
  normalization, energy, interpolation, and inverse identities.

## Chapter 4 — applications of the DFT

- Rewrote the full chapter, all twelve exercises, and all twelve solutions.
- **Continuous Fourier analysis:** corrected the Plancherel norm factor to
  sqrt(2 pi), fixed the proof that real-line characters are exponentials,
  restored the inverse-transform sign, and excluded the zero function from
  the compact-support impossibility statement. The DFT of delta_0 is the
  constant vector, not an exponential sampled with a nonzero frequency.
- **Sampling and quadrature:** corrected frequency storage order, the window
  phase, and the distinction between zero padding, longer observation, and
  denser sampling. Removed the false claim that periodic FFT quadrature is
  necessarily first order or intrinsically inaccurate; supplied the exact
  aliasing identity explaining its rapid convergence for smooth periodic data.
- **Filtering:** distinguished impulse response from transfer function, stated
  convergence and stability conditions, corrected the frequency-response sign,
  and made Gaussian normalization, truncation, and boundary handling explicit.
- **Polygon geometry:** corrected the iteration exponent, midpoint direction,
  centroid and regular-polygon radius factors 1/N, the area sign, and every
  Plancherel factor. Conditions for divergence or nonconvergence now apply to
  active frequencies. Clarified asymptotic periodicity, vertex versus area
  centroids, and the fact that magnitude descriptors define a pseudometric.
- **Heat equation:** separated continuous and discontinuous initial data,
  distinguished a periodized Gaussian from an ordinary Gaussian, and removed
  the false equivalence between arbitrary symmetric low-pass filters and heat
  evolution. Fixed the example that evaluated time zero at every requested time.
- **Heat schemes:** corrected the exact grid-dependent stability bounds, the
  unconditional theta threshold (1/2, not 1/4), the reversed stability inequality,
  and the two-dimensional explicit bound (1/4, not 1/2). Distinguished stable
  Crank–Nicolson oscillations from divergence. Corrected the expansion of CFL
  to Courant–Friedrichs–Lewy.
- **Poisson solver:** restored the missing -h^2/4 factor in the Fourier division,
  defined odd extensions consistently, and enforced boundary compatibility.
  Corrected the three-dimensional stencil center from +6 to -6 and restored
  h^-2. Equipotential surfaces are level sets of u, not f.
- **Sine transforms:** V is scaled orthogonal, with V^T V=(N/2)I, not orthogonal.
  Restored h^-2 in its eigenvalues and corrected the -2i and -4 factors relating
  one- and two-dimensional sine transforms to DFTs of odd extensions. Complex
  inputs are preserved rather than reduced to an imaginary or real part.
- **Correlation:** repaired missing coordinate differences, the squared-distance
  expansion and its -2 correlation term, and the direction of optimization.
  Restricted comparisons to complete template windows; restored powers in
  sliding-sum recurrences and used the Frobenius norm. Correlation one means
  positive proportionality after centering, not equality. Constant patches and
  templates are handled as undefined, with numerical cancellation accounted for.
- **FFT rotation:** repaired the rotation matrices and shear signs, specified
  noninteger interpolation and centered signed frequencies, and removed the
  false finite-grid rotation composition law. Odd grids avoid silently deleting
  the Nyquist component; complex values are preserved. Invertibility of the
  shear product is verified independently of the false group-law claim.
- **Products and interpolation:** distinct evaluation nodes are now assumed
  before identifying the kernel ideal; both Lagrange products exclude the
  interpolation node itself. Repaired a missing minus in the example polynomial,
  polynomial degree endpoints, and the base expansion. Exact digit conversion
  avoids logarithms and fixed-width overflow; carry propagation handles its tail.
  The floating-point FFT example verifies the product with exact Python integers.
- **Candy sharing:** supplied a finite-consensus proof using minimum sets and
  distinguished odd cycles from even-cycle oscillations for two-neighbor sharing.
- **Splines:** fixed shifted arguments, sample notation, support endpoints, and
  all extraneous 2 pi factors. Supplied a full positivity proof for both parities
  of spline degree, retaining alternating signs, and a dominated-convergence
  proof of the L2 limit to sinc. Natural endpoint conditions concern second,
  not first, derivatives. Periodic prefiltering is distinguished from natural,
  clamped, and not-a-knot finite-interval interpolation.
- **Programs:** converted and executed thirty additional numbered programs.
  Sixteen further checks cover stencils, exact Fourier modes, manufactured
  Poisson solutions and second-order convergence, sine-transform scaling,
  template matching against direct patch calculations, Fourier resampling,
  exact polynomial and integer products, polygon identities, spline symbols
  and convergence, periodic quadrature, and candy consensus. Fifty tests pass.
- **Inventory:** the original source has 67 figure environments but 68 numbered
  figures: Figures 4.10 and 4.11 shared one environment. Registered both separately
  so the final comparator will not silently omit Figure 4.11.
- **Figures:** regenerated all nineteen Chapter 4 figures from corrected calculations,
  including the two separately numbered Poisson solutions. The first four chapters
  and their solutions compile to a 100-page provisional PDF without box warnings;
  selected pages were visually checked. The 53-pair comparator compiles without warnings.

## Chapter 5 — Transforms related to the DFT

- Rewrote the chapter and all eleven exercises and solutions, preserving the
  cross-reference labels and distinguishing exact transforms from quadrature.
- **Hartley signs and scaling:** corrected H = Re(F) - Im(F), supplied the
  complex-linear version for complex inputs, and repaired the orthogonality
  proof without mixing normalized and unnormalized inner products. Restored
  the missing halves in the sine/cosine decomposition and the inversion scale.
- **Fast Hartley transform:** repaired half-length indices and paired-frequency
  twiddles. Removed unsupported universal speed, memory, and rounding-error
  claims relative to a real-input FFT. The original convolution routine
  overwrote its second input with zeros; the replacement transforms both inputs.
- **Correlation:** separated self-convolution from autocorrelation and supplied
  the correct Hartley formulas for each.
- **Hartley powers and generalization:** restored the orthogonal projector
  factors, separated unitary powers from unnormalized powers, and corrected
  the generalized-kernel inverse by a missing factor 1/N. Corrected its summation
  endpoint and identified the singular phases. Interpolation now uses signed
  frequencies, the correct dimension, and the post-padding normalization.
- **Z-transform convergence:** included absolute convergence and the region of
  convergence in the definition, covered empty interiors and boundary-only
  convergence, and repaired the convolution indices. Rational expressions alone
  do not identify causal versus anticausal sequences. Added coefficient inversion.
- **Recursive filters:** stated the zero-state convention, corrected causal and
  anticausal stability ranges, and proved the pole criterion for the reduced
  transfer function. Repeated poles in partial fractions are treated correctly.
- **Notch and RLC filters:** corrected the notch denominator sign and feedback
  coefficient indices, and the analog-to-digital substitution and all delay
  powers in the RLC transfer function. Laplace derivative formulas now account
  for initial conditions. An inverse DFT of a sampled response is a periodization,
  not an exact truncation of the impulse response.
- **Chirp algorithm:** excluded z=0, supplied the exact Vandermonde invertibility
  condition, fixed the chirp-kernel sign, and required a consistent logarithm.
  Derived the rectangular Toeplitz embedding with length at least N+M-1.
  Removed the claim that evenness of the kernel causes this length reduction.
- **Two meanings of fractional Fourier transform:** renamed the frequency-grid
  construction to distinguish it from unitary DFT powers. Its zero parameter
  returns the sum of the samples, not the identity. Corrected the sign in the
  Python transform so alpha=1 gives the forward DFT. Prime-length transforms
  are covered by Bluestein and Rader algorithms rather than declared infeasible.
- **Spectral leakage and adjustment:** distinguished the spectrum of a finite
  record from a continuous infinite sinusoid. A spectral zoom does not eliminate
  sampling or noise error. Replaced the erroneous zero-padding instruction by
  record shortening, derived its actual sidelobe formula and a valid local bound,
  and removed unsupported global bounds. The figure distinguishes the generating
  frequency from the estimated peak and compares normalized magnitudes.
- **Integration filters:** restored the sampling interval, corrected Simpson's
  two-step recurrence and its factor 1/3, and stated the initialization and its
  additional pole at -1.
- **Cubic spline recursions:** fixed the anticausal sign and derived exact mirror
  initialization by summing the periodic extension. Repaired both parallel and
  cascade implementations, including endpoint equations and the one-sample case.
  Both methods reproduce constants and agree with an independent spline filter.
- **Rader and Fourier quadrature:** corrected g^(a-b), input/output permutations,
  and modular inverses. Restored physical grid spacing, modulation and outer phase
  in Fourier-integral quadrature. Simpson's grid includes both endpoints and
  differs from the rectangle grid.
- **References:** corrected the Bailey–Swarztrauber journal citation to SIAM Review
  33(3), 389–404 (1991), and checked the author's errata. Recorded the spline
  reference's anticausal-initialization erratum.
- **Figures:** regenerated thirteen additional figures, including the RLC diagram
  in TikZ. Corrected the notch zeros (the original contour script placed them at
  radius 1.1 although the filter has unit-radius zeros), frequency indexing, and
  complex-spectrum displays. All were visually checked and their generators run
  with warnings treated as errors.
- **Programs:** converted fourteen further numbered programs. Thirteen additional
  tests cover direct transform matrices, convolution, signed interpolation,
  unitary powers, finite fields, rectangular Toeplitz products, prime FFTs,
  physical quadrature, exact mirror equations, and spectral shortening.
  The complete test suite now has 63 passing tests.

## Chapter 6 — Finite fields, rings, and error-correcting codes

- Rewrote the full chapter, all thirteen exercises, and their solutions.
- **Roots and extensions:** required transform length coprime to the characteristic;
  enlarging the extension cannot supply roots whose order is divisible by it.
  Distinguished a generator of the extension field from a generator of its
  multiplicative group. Corrected the constant term of Phi_6 and proved the
  cyclotomic-factor degree by Frobenius orbits.
- **Exact integer recovery:** bounded convolution coefficients explicitly. A larger
  extension cardinality alone does not increase the range of integer residues;
  the characteristic or the CRT product must exceed the appropriate bound.
- **Ring transforms:** clarified regularity, proved the root sums and regularity
  of the transform length, and stated the additional unit condition for inversion.
  Restored the missing +1 in Fermat moduli and numbers. Dyadic transforms use odd
  moduli of the form 2^a+1, not powers of two. The root-lifting proof now actually
  verifies the nth-power identity, and principal square roots are distinguished
  from arbitrary square roots in a ring with zero divisors.
- **Code definitions and Hamming examples:** corrected Hamming distance to w(x-y),
  distinguished length, dimension, and cardinality, and used a consistent column
  convention for generators. Recomputed the sixteen cyclic Hamming words and
  matched the parity-check matrix to their coordinate order. Corrected the
  systematic check matrix to [-A | I] over general fields.
- **Checksums:** corrected the ISBN-10 numerical example and distinguished its
  two-position transposition guarantee from the ISBN-13 rule used since 2007.
- **Cyclic and BCH codes:** defined the generator through a union of cyclotomic
  cosets, corrected the length-31 factor table and the designed-distance-five
  generator, and supplied a Vandermonde proof of the designed-distance bound.
  General finite-field magnitudes, Frobenius closure, and even designed distances
  are treated explicitly.
- **BCH decoding:** replaced the incorrect future-index recurrence and the claim
  that its coefficient matrix is triangular. The Hankel rank determines the actual
  error count; the smaller locator system is proved invertible by a Vandermonde
  factorization. This covers zero errors and fewer errors than the maximum radius.
  Spectrum completion and root/magnitude recovery are independently implemented
  and compared, including nonbinary errors.
- **Duality and enumerators:** restored n-m, X^(n-i), and the X-Y substitution
  throughout. Corrected mass-versus-average normalization and the rational
  coefficients of the length-eleven Hadamard example. Dual distance coefficients
  are nonnegative, not necessarily strictly positive or integral. Their total
  2^n/|C| does not require the code to contain zero.
- **Boolean functions and learning:** included the zero Walsh coefficient in the
  nonlinearity formula, corrected the count 2^(2^n), the squared-error expression,
  the direction of the sign-error bound, and the training-sample confidence.
  The simultaneous estimate uses the exact number of retained coefficients and
  a union bound; Walsh degree is distinguished from algebraic normal-form degree.
- **Coding exercises:** repaired the BCH realization of Hamming codes, the simplex
  and Hamming weight enumerators, the extended code's self-duality proof, and the
  repetition dual's minimum distance and characteristic-dependent self-duality.
  Walsh Hadamard codes are punctured linear/affine truth tables with different
  cardinalities. Restored the MDS moment identities and their triangular solution.
  The linear program now allows distance d and excludes only 1,...,d-1.
- **Python and validation:** converted fourteen more numbered programs, including
  the exact-arithmetic appendix. Ten new checks pass: exact transforms, composite
  moduli, large integer products, code enumerations, and exhaustive binary
  two-error patterns, together with nonbinary BCH examples.
- **Figure:** rebuilt the Hamming cube in TikZ with matching coordinate labels and
  a highlighted shortest path; inspected the rendered result.

## Chapter 7 — Linear representations of finite groups

- Rewrote the full chapter, all eleven exercises, and their solutions.
- **Group actions and representations:** distinguished the 24 proper cube rotations
  from all 48 symmetries; corrected inverse substitutions on polynomial and
  function spaces. Stated the unital algebra convention and the exact hypothesis
  for Maschke's averaging argument. Irreducible representations are nonzero.
- **Schur and Reynolds:** corrected the type of intertwining maps, the averaging
  projector and its invariant complement. Explained division algebras without
  calling them noncommutative fields. Restored complex conjugations in character
  inner products and the direction of the regular action.
- **Noncommutative Fourier analysis:** fixed the transform convention that respects
  convolution order and explained its relation to the earlier negative phase.
  Proved the center/class-function correspondence in both directions; distinguished
  canonical isotypic components from noncanonical individual copies.
- **Twists:** tensoring by a nontrivial scalar character need not produce a new
  irreducible representation; the standard representation of S3 is a counterexample.
- **Invariant polynomials:** repaired the multinomial coefficients, constant terms,
  and characteristic-zero divisions in Noether's bound. Supplied the missing
  quartic generator for the quarter-turn group and its algebraic relation.
  Restored the factor 1/|G| in Molien's formula and its inverse-matrix convention.
  Matching a rational denominator alone is not a proof of generation.
- **Counting and algebraic integers:** restored the missing half-turn in the
  necklace calculation, giving thirteen binary necklaces. Character values are
  sums of roots of unity, not necessarily roots of unity. Distinguished finite
  modules from finitely generated algebras in the degree-divisibility proof.
- **Group determinants:** restored the normalization from the Dirac basis and
  supplied the primitivity argument needed for irreducibility of a generic determinant.
- **Finite-field wavelets:** corrected affine-group multiplication, inverse action,
  modulation signs, conjugation and the frame constant. Separated reconstruction
  of the mean-zero subspace from the full tight-frame condition. The redundant
  coefficient map has a left inverse on its range, not a two-sided inverse on
  arbitrary arrays. Replaced an incorrect complexity claim with the FFT algorithm.
- **Validation:** five new tests check complex wavelets against direct sums, frame
  identities, affine composition, necklace counts, Hadamard codes, and linear
  programming bounds. Regenerated and visually checked the wavelet figure.

## Chapter 8 — Applications of representations

- Rewrote the complete chapter, all nine exercises, and their solutions.
- **Character tables:** restored complex conjugation in column orthogonality and
  placed the class-size matrix on the correct side. Corrected dihedral rotation
  angles, reflection indices, equivalent representations, and the even-order cases.
- **Cube actions and normal subgroups:** recomputed the face, vertex and edge
  characters and their decompositions. Distinguished the spatial representation
  from the standard S4 representation. Replaced the erroneous cyclic description
  of the Klein four-group by its four actual elements and the quotient S3.
  The character criterion for simplicity excludes the identity element.
- **Fourier coefficients:** restored the group-order factor and conjugation between
  normalized matrix-coefficient coordinates and matrix-valued Fourier blocks.
  Distinguished basis-dependent entries from invariant block energies and traces.
- **Ranking data:** specified that the three coordinates are the ranks assigned to
  city, suburb and countryside. Corrected the two three-cycle matrices and the
  resulting exact Fourier expansion. The dominant contrast concerns the last-ranked
  option; the totals are 987, 40 and 412, summing to 1439. The earlier index and
  permutation ambiguity changed the interpretation of the data.
- **Cauchy's theorem:** replaced rotation on all tuples by rotation on tuples whose
  product is the identity. Only that constraint makes fixed tuples give elements
  of order p. Completed the degree-two/simple-group argument, including faithfulness,
  determinant, diagonalizability and the center obstruction.
- **Quaternions:** separated the four-dimensional quaternion algebra from the
  eight-dimensional real group algebra, and identified the quotient ideal.
  Unitarity holds for unit quaternions, not arbitrary quaternions. Corrected the
  explicit Fourier basis, its conjugation and normalization, and product-group
  cardinality 8^n.
- **Cube invariants:** repaired the polynomial factors and their geometric meaning.
  Proved that the quartic itself is invariant and explained the additional
  degree-nine invariant for proper rotations.
- **Self-dual codes:** even word weights require diag(1,-1), not merely -I.
  The resulting invariant group has order sixteen. Proved generation by invariants
  of degrees two and eight, including algebraic independence, and recovered the
  extended Hamming enumerator P^4-4Q. The Python Reynolds routine returns a
  generating set that may be redundant; it does not claim automatic minimality.
- **Validation:** five new exact checks cover complex orthogonality, S4 decompositions,
  all S3 products and the ranking expansion, quaternion matrices, all 24 cube
  rotations, Reynolds averages, and both Molien series.

## Computational appendices and final editorial pass

- Replaced both language-specific appendices by numerical and exact-arithmetic
  Python appendices. Corrected claims about memory use and arbitrary-precision
  integers; documented module imports, indexing, phases, polynomial ordering,
  boundary conditions and the independent checks.
- All 88 displayed programs now come directly from executable Python functions.
  The test runner records which displayed functions actually execute; all are covered
  by the mathematical checks. No MATLAB or Maple execution is needed.
- Replaced the expensive symbolic cyclotomic search by exact truncated power-series
  arithmetic and justified the skipped indices. Independently verified all 27
  first-occurrence values, including the nonmonotone entries for 22 and 23.
- Corrected bibliography author names, titles, publishers and dates, including
  Boneh, Babai (2001), Bailey (1988), Cohen (1993), Cox–Little–O'Shea (1997),
  Elkies (2000), Kushilevitz–Mansour (1993), Pasalic–Johansson (1999), Vich,
  and Unser–Aldroubi–Eden. Kept only the 67 cited records in the active bibliography;
  the complete previous database remains in the legacy snapshot.
- Added separate navigation headings for the solutions to all eight chapters.
  Anchored programs near the text that introduces them and aligned figure-only
  pages at the top instead of spreading their contents vertically.

## Validation status

- Original edition archived and SHA-256 verified: 856 files.
- Mathematical and writing review: all eight chapters, 88 exercises and solutions,
  preface, and both computational appendices reviewed.
- Python verification: 84 tests pass; all functions printed in all 88 programs execute.
- Figures: all 68 regenerated and their generators run with warnings as errors.
- Complete PDF: 186 pages, with no LaTeX warnings, unresolved references or box warnings.
- Figure comparator: 69 pages, with all 68 figure numbers, captions and book pages
  matched against the compiled edition. All comparison pages visually reviewed.
- HTML presentation: revised cover, chapter links, Python download, figure comparator,
  correction log, and French-edition attribution; display and local links checked.
- Reproduction instructions, pinned Python dependencies, standalone Python download,
  and machine-readable verification reports included with the sources.
