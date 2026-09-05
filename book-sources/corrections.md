# Corrections to Mathematical Foundations of Data Sciences

Revision dates: 4–5 September 2026.

This revision covers all 19 chapters actively included by `FundationsDataScience.tex`, the six additional section files they load, and the preface. The record below identifies mathematical corrections, repaired arguments, qualified claims, figure changes, and grouped sentence-level edits. Section names and source labels provide locations that remain useful when pagination changes. Repeated punctuation, article, agreement, and spelling corrections are grouped rather than listed character by character.

Some original statements needed additional hypotheses or a different conclusion. Those changes are recorded explicitly below, including the revised ridge bound, compression theorem, Barron formulation, and slightly stronger sufficient condition in the constructive RIP certificate proof.

## Second editorial pass — 5 September 2026

Following the request for another pass, reread all 19 active chapters, their six section inputs, and the preface. This pass focuses on sentence structure, definitions, proof transitions, captions, terminology, and the precision of mathematical explanations. The two short optimization chapter wrappers required no changes; their included sections were revised throughout. The chapter records below group the sentence-level changes and identify every substantive qualification found during this pass. Selected before/after examples make the writing changes directly reviewable.

The displayed equations and theorem numbering are preserved. Inline corrections include the two-dimensional wavelet array dimensions and one figure's scale index; other mathematical changes clarify definitions, hypotheses, or explanatory models. The original correction record follows this second-pass record. The validation section describes the current PDFs.

### Preface, sampling, source coding, and compression

Read the preface and every active paragraph, proof, and caption in `shannon.tex`, `shannon-coding.tex`, and `compression.tex` again. Reworked complete sentences and paragraphs, retaining the established notation and displayed identities.

#### Preface (`abstract.tex`)

- Tightened the overview of the mathematical tools and their role in algorithm design.
- Clarified the relationship between the book and the Numerical Tours companion exercises.

#### Chapter 1: Shannon Sampling Theory (`shannon.tex`)

- Explained the sampler's blur/aliasing tradeoff in terms of the impulse response and sampling scale.
- Split the Plancherel extension argument into the norm identity, continuous extension, and meaning of the truncated integral limit.
- Simplified the Fourier inversion and coefficient definitions; clarified the distinctions among pointwise, normal, and uniform convergence.
- Connected the Sobolev discussion to approximation and denoising later in the book.
- Specified absolute almost-everywhere convergence of an integrable function's periodization and explicitly used **unnormalized Lebesgue measure** for its displayed `L¹` norm identity. This resolves ambiguity with the normalized inner product used earlier for Fourier series.
- Rewrote the Poisson proof to explain that decay of `f` makes its transform differentiable and compact spectral support makes periodization locally finite. Removed the misleading reference to a regularity implication in the opposite direction. Explained the sum–integral interchange directly.
- Clarified the reconstruction statement, rescaling argument, and spectral-view caption for the sampling theorem.
- Replaced the unsupported historical claim that Nyquist reproved Whittaker's theorem with a qualified account of Shannon's acknowledgment of earlier interpolation and transmission-rate work. Checked the original paper's discussion and footnotes in [Shannon, *Communication in the Presence of Noise* (1949), §II](https://course.khoury.northeastern.edu/csg250/ShannonNoise.pdf).
- Split the B-spline discussion into definition, support, regularity, interpolation equations, and the distinction between samples and spline coefficients. Explicitly identified the reconstruction kernel as `φ = φ_k`.
- Explained that finitely many samples do not determine a general bandlimited signal; exact finite-data reconstruction requires an additional finite-dimensional model or comparable assumption.
- Clarified midpoint dequantization and its worst-case error bound.

#### Chapter 3: Shannon Coding Theory (`shannon-coding.tex`)

- Simplified the alphabet and fixed-length-code examples and used “codeword” consistently.
- Made the prefix property and the benefit of shorter words for frequent symbols more direct.
- Defined entropy for a probability distribution rather than only for an empirical histogram.
- Separated differential entropy from the concavity calculation and removed the unnecessary suggestion that the density itself must be continuous.
- Introduced Kraft's inequality as a characterization of possible codeword lengths. Reworked both directions of the proof: disjoint descendant sets in the forward argument and aligned dyadic leaf blocks in the converse. Replaced the erroneous strict “less than” description of the packing capacity by the correct non-strict statement.
- Improved the Kraft and Huffman figure captions and renamed the probability-mismatch paragraph.
- Separated run-length modeling from the Markov entropy-rate discussion and made the performance qualification explicit.

#### Chapter 6: Compression (`compression.tex`)

- Rewrote the transform-coding overview and the comparison of thresholding, quantization, and coefficient reconstruction.
- Defined the comparison approximation `f_M` as retaining exactly the unquantized coefficients in `I_T`, with `M = #I_T`. Replaced the unsupported claim that the two errors are generally comparable by the existing precise additive error bound; the proof now identifies its two terms explicitly.
- Connected support coding to the bit budget rather than suggesting that a storage method is needed to measure quantization error itself.
- Added the bounded-domain hypothesis to the wavelet finite-access example: a bounded signal on the unit cube with periodic boundary conditions and periodized compactly supported wavelets. Explained the scale-wise coefficient count underlying `N(T) = O(T^-2)`; boundedness alone on all of `R^d` does not control the number of translations.
- Clarified how discrete transform coefficients relate to continuous inner products under the stated compatibility condition.
- Standardized “Entropy Coding” and rewrote the statistical-model and Huffman introductions. Allowed unequal probabilities to have equal codeword lengths, explained the integer-length overhead of symbol-by-symbol coding, and stated the arithmetic-coding comparison as an **expected** rate under the model.
- Refined the JPEG/JPEG-2000 comparison and presented dyadic thresholds as a schematic description of bit-plane refinement.
- Replaced claims of optimal packing at almost every rate with qualified coding-pass allocation and quality layers. Corrected square-only codeblocks, arbitrary-prefix decoding, stripe orientation, and the context description. These details were checked against [ITU-T T.800 (11/2015), §3.21, Annex D, and §J.13](https://www.itu.int/rec/dologin_pub.asp?id=T-REC-T.800-201511-S%21%21PDF-E&lang=e&type=items).
- Explained that contexts use significance/sign states and refinement history already known to the decoder, with rules specific to each bit type. Replaced the inaccurate generic previous-plane-bit formula by this description.

#### Representative sentence revisions

| Before this pass | After this pass |
|---|---|
| “The entropy of such a histogram is” | “The entropy of the probability distribution `p` is” |
| “Kraft's inequality, which describes the set of prefix codes using an inequality” | “Kraft's inequality characterizes the possible codeword lengths of a prefix code and is the key to the proof.” |
| “To measure the size of the additional error term … one needs to choose a method to store the quantized coefficients” | “To express the error bound … in terms of a bit budget, we must relate `M` and `T` to the cost of encoding the coefficients.” |
| “in practice these two errors have comparable magnitudes” | “The following bound separates the approximation error from the additional error caused by rounding.” |
| “Huffman coding algorithm in action.” | “Successive merges in Huffman's coding algorithm.” |

Displayed equations, source labels, figure assets, and citation keys were preserved. Inline mathematical changes are the explicit definitions and hypothesis clarifications recorded above; the removed context-bit formula was an inaccurate explanatory model, not a coding identity. Build and PDF validation are recorded in the consolidated report.

### Fourier, wavelets, approximation, and denoising

#### Fourier

- **Hilbert spaces and transforms, lines 10–104:** tightened the motivation, energy-preservation explanation, Gram–Schmidt construction, Shannon-basis example, and the transition between continuous and periodic transforms. Rewrote diagram descriptions and captions so that the operation or comparison is explicit.
- **Convolution and translation invariance, lines 128–291:** clarified the operator interpretation of Young's inequality, regularization by convolution, independence in the probability example, the distinction between continuous-output and distribution-kernel settings, and the two directions of the characterization proofs. Replaced vague antecedents with the named functional, kernel, or multiplier. Explained “translation equivariant” as input translations producing matching output translations.
- **Discrete transforms and interpolation, lines 323–597:** clarified sample and pixel counts, transform cost, the recursive FFT vectors and twiddle factors, matrix diagonalization, the impulse response, polynomial multiplication, and the distinction between approximating a continuous transform and exactly evaluating a trigonometric interpolant.
- **Several dimensions and differential operators, lines 623–832:** improved the tensor-product explanation, frequency–space pairing, row/column algorithm, operation count, and Fourier-convolution transition. Corrected the geometric descriptions listed below and made the finite-difference caption plural (“spectra”).
- **Groups, representations, and general domains, lines 856–1098:** refined the construction of the dual group, product factorization, finite/infinite cases, representation versus trace character, invariant complements, irreducible representations, change of basis, and the transition to Laplacian eigenfunctions on surfaces and graphs. Improved the source attribution to Diaconis and captions for the torus and surface Laplacian.

##### Mathematical terminology clarified

- **Lines 623–627:** the product in the displayed formula is a product of complex exponentials, hence a plane wave. Its level sets are orthogonal to the wave vector. The former description of a sine wave “moving orthogonally” confused the geometry with motion and described the wrong direction.
- **Line 826:** the displayed first difference is forward, but the displayed second difference is centered. The prose now names them separately; neither formula changed.
- **Line 917:** the finite-group isomorphism with the dual depends on generators and corresponding roots of unity, whereas the double-dual identification is canonically defined by evaluation. The revision makes this distinction explicit without changing either map.
- **Line 1089:** the spherical-harmonic index corresponds to frequency magnitude, not the “amplitude of Fourier frequencies.”

#### Wavelets

- **Multiresolution spaces and orthogonalization, lines 9–251:** refined the nested-space and projection explanations, sampling versus projection distinction, orthogonalization criterion, detail-space construction, convergence, and Haar/Shannon/spline comparisons. Caption wording now says what the projections and frequency partitions show.
- **Periodic and discrete transforms, lines 264–564:** made the period convention explicit, distinguished exact acquisition inner products from the approximation by point samples, stated when intermediate approximation coefficients can be discarded, and clarified reflected filters, periodic convolution, recursion, adjoints, and the inverse transform.
- **Two-dimensional transforms, lines 589–799:** clarified separable row/column computation, tensor-product bases, the three detail arrays, sample counts per direction, coefficient storage, and the input to the algorithm.
- **Filter design and cancellation, lines 881–1205:** distinguished a filter's periodic Fourier series from the scaling function's continuous Fourier transform, stated which direction of the construction theorem is proved, improved transitions in the proof, and clarified filter complementarity, vanishing moments, localization, smoothness, and the interpretation of the Daubechies examples.

##### Mathematical terminology and dimensional corrections

- **Lines 698, 792, and 799:** this section defines `N = 2^{-J}` samples in each coordinate direction. Its coefficient array therefore has `N²` pixels. Corrected the storage sentence from `$N$` to `$N^2$` and the 2-D algorithm input from `$f \in \CC^N$` to `$f \in \CC^{N\times N}$`. The final visual review also corrected the scale index in the first panel label of Figure 4.22, as recorded below.
- **Line 1075:** vanishing moments quantify cancellation of polynomials, rather than literally counting oscillations.
- **Line 1126:** admissibility guarantees **at least** one vanishing moment. Additional moments correspond to higher-order vanishing of the filter's Fourier transform at the specified frequency; “flatter” was replaced with this precise description.

#### Approximation

- **Basis approximation and image models, lines 6–315:** tightened the chapter motivation, linear/nonlinear index selection, ties in coefficient thresholding, coefficient-decay interpretation, Sobolev seminorm explanation, piecewise smooth models, total variation, cartoon geometry, and blur model.
- **Efficiency and applications, lines 338–606:** clarified the roles of the approximation exponent and signal-dependent constant, distinguished comparison across images from comparison across bases, improved the compression/denoising/inverse-problem transitions, and confined claims about basis performance to the displayed examples. Explained finite-sample saturation and the relation between smooth and Sobolev models.
- **Wavelet rates, lines 633–854:** clarified the role of vanishing moments, regular versus singular index sets, coefficient counting, cutoff scales, comparison approximations, the difference between 1-D and 2-D singularities, and the limits of isotropic supports for cartoon images.
- **Triangulations and curvelets, lines 873–1013:** refined the construction and implementation discussion, anisotropic dimensions, translation and rotation, orientation-dependent sampling, tight-frame thresholding, cancellation near and away from edges, redundancy, and the transition to denoising. Rewrote figure captions and standardized nonlinear terminology.

##### Hypotheses, interpretation, and reference correction

- **Lines 588 and 608:** the previous prose attached the jump-dominated Fourier error bounds `O(M^{-1})` in 1-D and `O(M^{-1/2})` in 2-D to the whole piecewise `C^alpha` model. That model allows regular pieces whose smoothness is too low for these rates. The prose now gives sufficient assumptions: piecewise Lipschitz regularity, together with finite-length edges in 2-D. The displayed bounds are unchanged. This deliberately states a sufficient regime rather than asserting a sharp smoothness threshold.
- **Line 357:** efficient sparse approximation alone does not guarantee inverse-problem recovery. The revision connects recovery explicitly to the measurement operator retaining enough information about sparse combinations, an existing qualification that previously appeared only after a broader claim.
- **Line 941:** the sentence about spatial and frequency localization referred to `fig-curvelet-discretization`. Corrected it to `fig-curvelet`, the figure that actually displays that localization. Both figure definitions and all assets are unchanged.
- **Line 975:** approximately `5N` curvelet atoms describes one cited discrete construction, not an intrinsic redundancy shared by every discretization.
- **Line 984:** the frame-thresholded approximation need not be the best M-term approximation because the frame is not orthogonal; the revision keeps this qualification directly beside the construction.

#### Denoising

- **Models and risk, lines 8–101:** distinguished random variables from observed realizations and the deterministic estimator from its random output. Made the acquisition model, Gaussian independence, theoretical risk, experimental SNR, and need for a clean reference in benchmarking explicit.
- **Filtering, oracle, and Wiener estimators, lines 119–352:** clarified the bias–detail tradeoff in smoothing, tuning within a filter family, the use of an oracle to optimize SNR, diagonal estimation, the distinction between known clean coefficients and a random-signal prior, stationarity, and the transfer of approximation rates to sampled data.
- **Thresholding, lines 379–537:** improved the sparse-coefficient motivation, hard/soft rules, coarse-coefficient treatment, sparsity explanation, Gaussian maximum discussion, and universal-threshold interpretation.
- **Translation invariance and shrinkage variants, lines 564–831:** clarified periodic shifts, cycle spinning, frame redundancy, repeated atoms, semi-soft and Stein shrinkage, block dependence, and experiment captions. Tied parameter choices and rankings to the examples actually shown.
- **Signal-dependent noise, lines 861–1034:** clarified Poisson counts, intensity-dependent variance, variance stabilization as an approximation, random multiplicative factors, independent SAR looks, Gamma distributions, and the scope of displayed denoising comparisons.

##### Statistical interpretation and scope clarified

- **Lines 93 and 96–101:** the denoising map is deterministic; its output is random. The observation has mean equal to the clean signal in the additive model. Risk averages over noise, whereas the plotted SNR evaluates an observed realization against known clean data.
- **Lines 207, 213, and 220:** the SNR-optimal filter width is an experimental oracle choice that uses clean data. The 1-D figure is correctly captioned “Noisy signal,” replacing “Noisy image.”
- **Line 392:** after a real orthogonal transform, white-noise coefficients remain an independent family of centered Gaussians, rather than describing an individual scalar coefficient as “white noise.”
- **Line 452:** soft thresholding reduces retained coefficient **magnitudes**, which describes both positive and negative coefficients correctly.
- **Lines 467, 613, 771, 821–822, and 831:** limited claims about optimal thresholds, hard versus soft thresholding, Stein rules, block sizes, and relative reconstruction quality to the displayed experiments. The numerical figures remain unchanged; they are evidence for those examples rather than universal rankings.
- **Lines 524 and 537:** distinguished asymptotic concentration of the Gaussian maximum from numerical illustration and retained the logarithmic-factor interpretation of the universal-threshold risk bound.
- **Lines 763 and 789:** preserved the distinction between deterministic shrinkage vanishing at large amplitudes and statistical unbiasedness, and explained block attenuation in terms of grouped coefficients rather than isolated scalar decisions.
- **Lines 894 and 896:** the relative Poisson standard deviation applies to positive intensities. Replaced the description of clean-image “quantization” with different maximum mean intensities; no unavailable numerical acquisition values were inferred.
- **Lines 917–923 and 946:** variance stabilization gives an approximate noise model whose quality depends on intensity; the observed improvement is scoped to the moderate count levels in the example.
- **Line 985:** averaging K independent looks divides the additive variance by K. Approximating the result by multiplicative noise alone additionally requires the remaining additive component to be negligible. Increasing K by itself does not justify discarding that component.

#### Five representative before/after examples

1. **Fourier, line 826**
   - Before: “one typically considers forward finite differences (first and second order)”
   - After: “consider the forward first difference and the centered second difference”

2. **Wavelets, line 308**
   - Before: “During the algorithm, the previously computed vector a_j can be discarded, and only the d_j are kept.”
   - After: “After an approximation vector a_j has been used to compute the next scale, it can be discarded; the detail vectors d_j are retained.”

3. **Approximation, line 357**
   - Before: “A basis that is efficient for approximating the high resolution signal is needed to recover missing information efficiently.”
   - After: “Efficient approximation of the unknown signal helps recover missing information, provided the measurement operator retains enough information about sparse combinations of basis atoms.”

4. **Denoising, line 93**
   - Before: “so that D is a deterministic function. It is thus also a random vector that depends on the noise w.”
   - After: “The map D is deterministic, but its output is random because it depends on the noise w.”

5. **Denoising, line 467**
   - Before: “These results also show that numerically, for thresholding in orthogonal bases, soft thresholding is slightly better than hard thresholding on natural signals and images.”
   - After: “The best soft-thresholded estimate has a slightly higher SNR.” The preceding sentence now explicitly limits the comparison to this experiment.

Mathematical notation is simplified in these short prose examples only; the source preservation checks use the exact TeX.

#### Gaussian plot label correction

Final visual QA identified overlap between the x-axis label and the empirical curve in Figure 7.16 (PDF page 92 before this repair). Added `figures/denoising/max-gaussian-std-corrected.pdf`, retained the original, and changed only its include at `chapters/denoising.tex:528`. The optional pypdf script `scripts/correct_max_gaussian_std.py` preserves the embedded empirical image and all original plotting commands, translates the separate outlined label below the plot, and enlarges the lower PDF bound. No sample values were inferred or regenerated. Image bytes and ICC profile are unchanged; 2× rendering comparison shows zero changes outside the old label region. Neither asset uses font resources, so the outlined label creates no font-embedding issue. The script regenerates the checked-in PDF deterministically and adds no dependency to the normal LaTeX build. Exact geometry, method, and visual evidence are recorded in `build/pass2-reviews/qa-signal.md`.

### Inverse problems and sparsity

#### Variational priors

Source: `chapters/variational-priors.tex`; 57 revised paragraphs or lines.

Reread all active prose, captions, definitions, and derivations. Reworked the motivation for image priors, the passage from continuous domains to pixel arrays, the definitions of discrete operators, the description of diffusion flows, and the transition from flows to regularized estimators. Shortened repetitive explanations of parameter selection and filtering.

Clarifications affecting mathematical meaning:

- Describe functions on continuous domains without implying that membership in `L²` makes a function continuous. Describe a discrete image as a pixel array in a finite-dimensional space.
- Say that the prior energy is real valued **when finite**, consistently with the later extension by `+∞`.
- Describe the smoothness classes as motivating Sobolev priors, avoiding an unsupported exact identification of classes.
- Identify the displayed Sobolev energy as **one half** of the squared `H¹` seminorm, matching its existing normalization.
- Identify divergence as the negative adjoint of the gradient. Clarify the grid-spacing scaling needed when comparing the discrete Laplacian with its continuous counterpart.
- Distinguish the ODE for a finite-dimensional image from the PDE when space is continuous.
- Describe Figure `fig-discrete-operators-tveps` as showing regularized gradient magnitude, rather than the gradient of the energy.
- Explain the distinction between choosing a stopping time for a flow and choosing a regularization parameter for a minimization problem.
- Explicitly retain **properness, lower semicontinuity, and lower boundedness** as sufficient conditions for a nonconvex denoising objective to attain a minimum. Properness was missing from the nonconvex existence sentence; uniqueness is still asserted only under the stated convexity assumptions.

##### Selected before/after examples

**Line 23. Before:**

```tex
The prior energy $J(f) \in \RR$ is intended to be low for images in a class $f \in \Theta$. The class of uniformly smooth functions detailed in Section \ref{subsec-smooth-class} corresponds to functions in Sobolev spaces. A simple prior derived from this Sobolev class is thus
```

**After:**

```tex
When finite, the prior energy $J(f) \in \RR$ should be small for images in the target class $f \in \Theta$. The smoothness classes described in Section \ref{subsec-smooth-class} motivate Sobolev priors. The simplest such energy is
```

**Line 27. Before:**

```tex
where $\nabla f$ is the distributional gradient and the integral is over $[0,1]^2$. This is the squared $H^1$ seminorm: it vanishes on constant functions. We set the energy to $+\infty$ when the gradient is not square integrable.
```

**After:**

```tex
where $\nabla f$ is the distributional gradient and the integral is over $[0,1]^2$. The energy is one half of the squared $H^1$ seminorm and vanishes on constant functions. We set it to $+\infty$ when the gradient is not square integrable.
```

**Line 407. Before:**

```tex
Instead of performing a gradient descent flow for denoising as detailed in Section \ref{subsec-pde-denoise} and selecting a stopping time, one can formulate an optimization problem. The estimator is then defined as a solution of this optimization. This setup has the advantage of being well defined mathematically even for nonsmooth priors such as the TV prior $\Jtv$ or the sparsity prior $\Ju$. Furthermore, this regularization framework is also useful to solve general inverse problems as detailed in Chapter \ref{chap-invpbm}.
```

**After:**

```tex
The flow-based approach in Section \ref{subsec-pde-denoise} defines the estimator by choosing a stopping time. Alternatively, we can define it as the minimizer of an objective that balances data fidelity and a prior. This formulation also accommodates nonsmooth priors such as total variation $\Jtv$ and sparsity $\Ju$. Chapter \ref{chap-invpbm} extends it to general inverse problems.
```

**Line 489. Before:**

```tex
Denoising by Sobolev regularization is a special case of the linear estimator considered in Section \ref{sec-linear-denoising}. The selection of the parameter $\la$ is related to the selection of an optimal filter as considered in Section \ref{subsec-optimal-selection-filter}, but with the restriction that the filter is computed in a parametric family.
```

**After:**

```tex
Sobolev denoising belongs to the class of linear estimators studied in Section \ref{sec-linear-denoising}. Choosing $\la$ is therefore a filter-selection problem, as in Section \ref{subsec-optimal-selection-filter}, restricted here to a one-parameter family.
```

#### Inverse problems

Source: `chapters/inverse-problems.tex`; 96 revised paragraphs or lines.

Reread every active paragraph, caption, theorem, and proof. Rewrote the acquisition examples, diagonal estimators, SVD construction, source-condition discussion, bias–variance explanation, regularization limit, and iterative solvers. Broke long proof sentences into the bound, its consequence, and the resulting conclusion. Improved the introduction of regression notation and the distinctions between the various reconstruction methods.

Clarifications affecting mathematical meaning:

- Call the domain `Ss` the **signal space**, rather than the data space.
- Present tomography and MRI through the stated idealized linear/Fourier models, avoiding categorical identification of physical acquisition with a simplified model.
- Explain the correspondence with ridge-regression notation and how the regression matrix is assembled.
- Describe the SVD eigendecomposition on the relevant range; call the matrix positive semidefinite and say that the **columns** of `U` form the basis.
- Explain the reduction to the orthogonal complement of the kernel and the role of rank. Tie the unbounded pseudoinverse statement to the compact, infinite-rank setting under discussion.
- Attribute frequency decay to the kernel's Fourier multiplier. Separate the convolution operator, its eigenfunctions, its multipliers, and its singular values in the exposition.
- Restrict the monotonicity description of bias and variance to **Tikhonov regularization**, rather than arbitrary spectral filters.
- Describe the zero-parameter result as convergence of minimizers under its stated assumptions, rather than literal equivalence of the penalized and constrained optimization problems.
- Explain coercivity using **sublevel sets**, not level sets.
- State that the conjugate-gradient quadratic differs from the original objective only by an additive constant; identify positive definiteness as the relevant hypothesis for the stated Cholesky factorization.
- Describe the smoothed-TV limit as a rescaled Sobolev energy, up to an additive constant.
- Specify a smooth prior when introducing the corresponding projected-gradient scheme.
- Describe the tomographic acquisition as using equally spaced orientations, rather than equally spaced rays.

##### Selected before/after examples

**Line 118. Before:**

```tex
We first consider a simpler setup where we seek an ``optimal'' inversion method which is assumed to be diagonal in some basis. We will see later that these estimators can also be seen as solving a multi-dimensional optimization problem (here, the optimization is just 1-dimensional over each basis coefficient).
```

**After:**

```tex
We first seek an optimal estimator among inversion methods that act diagonally in a fixed basis. The optimization then separates into one scalar problem per coefficient. Later we will obtain such estimators from a single optimization problem on the signal space.
```

**Line 257. Before:**

```tex
An example of such a setting which generalizes~\eqref{eq-diag-filter-finite} is when $\Phi f = f \star h$ on $\TT^d=(\RR/2\pi\ZZ)^d$, which corresponds to a translation-invariant kernel $k(x,y) = h(x-y)$, in which case $q_m(x)=(2\pi)^{-d/2}e^{\imath m\cdot x}$ diagonalizes the convolution, and $\sigma_m=|\hat h_m|$ when $\hat h_m=\int_{\TT^d}h(x)e^{-\imath m\cdot x}\,\mathrm{d}x$ denotes its multiplier.
```

**After:**

```tex
Convolution $\Phi f = f \star h$ on $\TT^d=(\RR/2\pi\ZZ)^d$ generalizes~\eqref{eq-diag-filter-finite}. Its kernel $k(x,y) = h(x-y)$ is translation invariant, and the Fourier functions $q_m(x)=(2\pi)^{-d/2}e^{\imath m\cdot x}$ diagonalize the operator. Its singular values are $\sigma_m=|\hat h_m|$, where $\hat h_m=\int_{\TT^d}h(x)e^{-\imath m\cdot x}\,\mathrm{d}x$ denotes the convolution multiplier.
```

**Line 567. Before:**

```tex
This shows that $J(f_\la) \leq J(h)$ so that since $J$ is coercive the set $(f_\la)_\la$ is bounded and thus one can consider an accumulation point $f_{\la_k} \rightarrow f^\star$ for $k \rightarrow +\infty$.  Writing $m=\inf J> -\infty$, since $\norm{\Phi f_{\la_k}-y}^2\leq2\la_k(J(h)-m)$ and $\la_k\to0$, one has in the limit $\Phi f^\star = y$, so that $f^\star$ satisfies the constraints in~\eqref{eq-ip-regul-noiseless}. Furthermore, by lower semicontinuity of $J$, passing to the limit inferior in $J(f_{\la_k}) \leq J(h)$, one obtains $J(f^\star) \leq J(h)$ so that $f^\star$ is a solution of~\eqref{eq-ip-regul-noiseless}.
```

**After:**

```tex
The inequality $J(f_\la) \leq J(h)$ and coercivity of $J$ make $(f_\la)_\la$ bounded. Choose a convergent subsequence $f_{\la_k} \rightarrow f^\star$ as $k \rightarrow +\infty$. Write $m=\inf J> -\infty$. The estimate $\norm{\Phi f_{\la_k}-y}^2\leq2\la_k(J(h)-m)$ and the limit $\la_k\to0$ imply $\Phi f^\star = y$, so $f^\star$ is feasible for~\eqref{eq-ip-regul-noiseless}. Lower semicontinuity of $J$, applied to $J(f_{\la_k}) \leq J(h)$, gives $J(f^\star) \leq J(h)$. Thus $f^\star$ solves~\eqref{eq-ip-regul-noiseless}.
```

**Line 644. Before:**

```tex
The crucial and remarkable fact is that this minimization can be computed in closed form at the cost of one matrix-vector product per iteration, for $\ell\geq1$ (posing initially $\itz{d} = \nabla \Ee(\itz{f}) = A\itz{f}-b$)
```

**After:**

```tex
This subspace minimization requires only one matrix-vector product per iteration. For $\ell\geq1$, with initialization $\itz{d} = \nabla \Ee(\itz{f}) = A\itz{f}-b$, the recurrence is
```

#### Sparse regularization

Source: `chapters/sparse-regularization.tex`; 49 revised paragraphs or lines.

Reread every active paragraph, caption, thresholding argument, and algorithm derivation. Refined the motivation for sparse coefficients and compressibility, thresholding conventions, dictionary notation, the positivity reformulation, ISTA, and the imaging/seismic examples. Replaced repeated “slight abuse of notation” and “so-called” constructions with direct definitions.

Clarifications affecting mathematical meaning:

- Distinguish approximate sparsity from exact zero coefficients in the text and wavelet caption.
- Qualify the BV approximation statement as applying to **suitable wavelet bases**, rather than arbitrary orthonormal bases.
- Make the `l0` nonconvexity example use disjoint **nonempty** coefficient supports, excluding the trivial zero-vector case.
- Describe `q=1` as the boundary of the convex range; remove an unqualified assertion that it achieves the “highest degree of sparsity.”
- Treat zero input explicitly in the hard-thresholding argument before comparing with the nonzero candidate.
- Explain that the dictionary symbol also denotes the synthesis operator and that its adjoint is the analysis operator.
- Clarify why the surrogate is an upper bound and is tangent at the current iterate. Explain the step-size restriction and ISTA's storage advantage over the positive/negative splitting.
- Identify the seismic unknown as **subsurface reflectivity** and the convolution kernel as the transmitted pulse; distinguish a seismic wavelet from an orthogonal wavelet basis.
- Distinguish synthesis ISTA for a redundant dictionary from simply applying analysis, thresholding, and synthesis once.

##### Selected before/after examples

**Line 95. Before:**

```tex
where $x_m \eqdef \dotp{f}{\psi_m}$, $y_m \eqdef \dotp{g}{\psi_m}$, and where we use the following slight abuse of notation for $q=0$
```

**After:**

```tex
Here $x_m \eqdef \dotp{f}{\psi_m}$ and $y_m \eqdef \dotp{g}{\psi_m}$. For $q=0$, we adopt the convention
```

**Line 189. Before:**

```tex
In the following, with a slight abuse of notation, we denote the ``synthesis'' and ``analysis'' operators as
```

**After:**

```tex
We also use the dictionary symbol for its synthesis operator; its adjoint is the analysis operator:
```

**Line 306. Before:**

```tex
A drawback of this projected gradient descent scheme is that it requires storing $2Q$ coefficients. A closely related method, with comparable convergence guarantees, is the so-called ``iterative soft thresholding algorithm'' (ISTA).
```

**After:**

```tex
Splitting positive and negative parts requires storing $2Q$ coefficients. The iterative soft thresholding algorithm (ISTA) avoids this duplication while retaining comparable convergence guarantees.
```

**Line 391. Before:**

```tex
In a simplified linearized one-dimensional model, ignoring multiple reflections, the acquisition of underground data $f_0$ is modeled as a convolution $y=h \star f_0 + w$, where $h$ is a so-called ``wavelet'' signal transmitted into the ground. This should not be confused with the construction of orthogonal wavelet bases detailed in Chapter \ref{chap-wavelets}, although the term ``wavelet'' originally comes from seismic imaging.
```

**After:**

```tex
In a linearized one-dimensional model that neglects multiple reflections, the subsurface reflectivity $f_0$ produces observations $y=h \star f_0 + w$. Here $h$ is the transmitted pulse, called a seismic wavelet. This use of the word differs from the orthogonal wavelet bases constructed in Chapter \ref{chap-wavelets}, although the terminology originated in seismic imaging.
```

#### Sparse recovery theory

Source: `chapters/sparse-theory.tex`; 88 revised paragraphs or lines.

Reread all active statements, proofs, explanatory paragraphs, and captions, including the duality, Bregman, support-stability, and continuous-domain sections. Reworked the sequence of definitions and conclusions, replaced ambiguous uses of “it” and “this,” and made the distinctions between a primal solution, dual solution, and certificate explicit.

Clarifications affecting mathematical meaning:

- Explain existence using continuity and coercivity; recall that the `l1` norm is not strictly convex.
- Describe polytope faces that are lost under projection, rather than an imprecise “small projected polytope.” Avoid suggesting that every signal is uniquely identifiable.
- Describe the empty cap as the cap bounded by the circumcircle, rather than the circle's interior.
- Say that certificates arise from Lagrange multipliers, keeping the data-space multiplier distinct from its image, the certificate.
- Explain that the certificate set is independent of the chosen primal solution, and point to the duality argument below rather than incorrectly to the next chapter.
- State the Lasso sparsity conclusion as existence of **a solution** supported on linearly independent columns, rather than as a property of every solution.
- Correct the direction of the Bregman-figure caption: under the stated certificate assumptions, the Bregman divergence controls the `l1` error. Explain how the strict saturation margin controls off-support coordinate error.
- Separate the nonlinear source condition, which concerns a subgradient, from the recovery method itself.
- Identify the limit of `p_lambda` as the dual solution of minimum norm; the associated `eta` is the certificate. Explain that the Fuchs candidate may fail dual feasibility.
- Say that dual certificates **may** be nonunique, rather than asserting they always are.
- Make the role of a strict off-support certificate bound in support stability explicit; distinguish it from ordinary feasibility used for recovery.
- Identify `T` as the smallest nonzero amplitude and explain the columnwise correlations of the noise perturbation.
- Use **Gram matrix** for the matrix of column inner products and **kernel** for its continuous analogue, rather than “covariance matrix/function.”
- Describe the discrete certificate as samples of a continuous linear combination, not as a continuous function itself.
- Explain the zero-derivative condition using magnitude at most one together with sign interpolation at an interior point; boundedness alone would not imply the condition.
- Explain that the coefficient `l2` norm on a refining grid is not an intrinsic distance between measures, replacing the unqualified statement that it is meaningless.

##### Selected before/after examples

**Line 129. Before:**

```tex
In the following, given an index set $I \subset \{1,\ldots,N\}$, denoting $A=(a_i)_{i=1}^N$ the columns of $A$, we denote $A_I \eqdef (a_i)_{i \in I} \in \RR^{P \times |I|}$ the extracted submatrix.  Similarly, for $x \in \RR^N$, we denote $x_I \eqdef (x_i)_{i \in I} \in \RR^{|I|}$.
```

**After:**

```tex
For an index set $I \subset \{1,\ldots,N\}$, write $A=(a_i)_{i=1}^N$ for the columns of $A$ and $A_I \eqdef (a_i)_{i \in I} \in \RR^{P \times |I|}$ for the corresponding submatrix. For a vector $x \in \RR^N$, write $x_I \eqdef (x_i)_{i \in I} \in \RR^{|I|}$ for its restriction.
```

**Line 181. Before:**

```tex
Although it looks like the definition of $\Dd_0(y,x^\star)$ depends on the choice of a solution $x^\star$, convex duality (studied in the next chapter) shows that it does not (it is the same set for all solutions).
```

**After:**

```tex
Although the notation $\Dd_0(y,x^\star)$ involves a particular solution $x^\star$, the set of certificates is the same for every primal solution. The duality argument below explains this independence.
```

**Line 431. Before:**

```tex
The issue with the control~\eqref{eq-bregman-rate} of the error in terms of Bregman divergence is that it is not ``distance-like'' for regularizers $J$ that are not strictly convex. This is in particular the case for the $\ell^1$ norm $J=\norm{\cdot}_1$ which we now study.
```

**After:**

```tex
For a regularizer $J$ that is not strictly convex, the Bregman bound~\eqref{eq-bregman-rate} need not control the error in norm. We now examine how additional certificate conditions restore such control for the $\ell^1$ penalty $J=\norm{\cdot}_1$.
```

**Line 554. Before:**

```tex
The following proposition shows that $p_\la$ converges to a specific dual certificate of the constrained problem, called ``minimal norm'' certificate.
```

**After:**

```tex
The next proposition shows that $p_\la$ converges to the dual solution of minimum norm for the constrained problem.
```

#### Compressed sensing

Source: `chapters/compressed-sensing.tex`; 39 revised paragraphs or lines.

Reread every active paragraph, caption, recovery statement, and proof, including the acquisition model, coherence and RIP estimates, concentration arguments, interpolation estimate, and structured sampling. Refined quantifiers and removed informal or unsupported comparisons between different guarantees.

Clarifications affecting mathematical meaning:

- Model the single-pixel scene as a light-intensity function; explain the roles of noise and reconstruction directly.
- Explain the relationship between a noise bound and a Lagrange multiplier, including limiting cases, without asserting a one-to-one correspondence.
- Distinguish the quantifiers in uniform/strong and nonuniform/weak recovery statements.
- Keep the distinction between coherence's entrywise control and the operator-norm control needed for recovery arguments.
- In the conditional concentration proof, explicitly say that the exponential quantity bounds the **probability** that a correlation reaches or exceeds one.
- Use Gram-matrix **eigenvalues** and measurement-matrix **singular values** consistently in the RIP explanation.
- Distinguish finite-sample concentration and the union bound from the asymptotic Marchenko–Pastur limit.
- Explain iterative certificate corrections without suggesting that the initial Fuchs certificate always suffices, or assigning an unsupported causal explanation to the logarithmic-factor improvement.
- Describe the interpolation tail as **at most** the threshold, consistently with the existing weak inequality.
- State that the coordinate change uses an **orthonormal** basis, as required by the existing adjoint-based formula.
- Describe structured sampling as selecting random coefficients in an orthonormal basis, removing the vague phrase “less random.”
- Keep uniformity, norm stability, and support stability distinct; remove the unsupported causal connection between logarithmic improvements and support stability.

##### Selected before/after examples

**Line 25. Before:**

```tex
The \guill{single-pixel} hardware performs the compressed sampling of an observed scene $\tilde f_0$ (the letter \guill{R} in Figure~\ref{fig-single-pixel}), which is a continuous function indicating the amount of light $\tilde f_0(s)$ reaching each point $s \in \RR^2$ of the focal plane of the camera.
```

**After:**

```tex
The scene $\tilde f_0$, represented by the letter \guill{R} in Figure~\ref{fig-single-pixel}, is modeled as a light-intensity function. Its value $\tilde f_0(s)$ gives the intensity at position $s \in \RR^2$ in the camera's focal plane.
```

**Line 192. Before:**

```tex
We first perform a crude analysis using the so-called coherence of the matrix $A=(a_j)_{j=1}^N$ where the $a_j \in \RR^P$ are the columns of $A$, which we assume to be normalized $\norm{a_j}=1$
```

**After:**

```tex
We begin with a bound based on the coherence of $A=(a_j)_{j=1}^N$. The columns of $A$, denoted $a_j \in \RR^P$, are normalized by $\norm{a_j}=1$, and the coherence is
```

**Line 303. Before:**

```tex
Conditional on $A_I$, the vector $p_F$ is fixed and the columns $a_j$, $j\notin I$, remain independent of it. Sub-Gaussian concentration bounds each $|\langle a_j,p_F\rangle|$ above $1$ by $2\exp(-cP/\norm{p_F}^2)$. A union bound over at most $N$ columns proves the stated scaling, with the strict inequality obtained by controlling the event $|\langle a_j,p_F\rangle|\geq1$.
```

**After:**

```tex
Condition on $A_I$. The vector $p_F$ is then fixed, while the columns $a_j$ for $j\notin I$ remain independent of it. The probability that $|\langle a_j,p_F\rangle|$ reaches or exceeds $1$ is bounded by $2\exp(-cP/\norm{p_F}^2)$. A union bound over at most $N$ columns gives the stated sample complexity and the strict off-support bound, by excluding $|\langle a_j,p_F\rangle|\geq1$.
```

**Line 584. Before:**

```tex
Implementing a matrix with fully independent random entries can be difficult in hardware. A more practical option is to use structured sampling operators, which are in some sense ``less random''. A possibility is to consider a random sub-sampling of orthogonal projections of the signal onto an orthonormal basis $\Xi = (\xi_\om)_{\om=1}^N$ of $\RR^N$ or $\CC^N$, so that
```

**After:**

```tex
Fully independent random measurements can be difficult to implement. A practical alternative is to sample a random subset of coefficients in an orthonormal basis $\Xi = (\xi_\om)_{\om=1}^N$ of $\RR^N$ or $\CC^N$:
```

#### Preservation checks and final validation

- Compared every changed line with the saved baseline: all inline mathematical expressions are preserved exactly, allowing their order within a paragraph to change.
- Compared all display equations and the custom equation blocks exactly: no formula changes.
- Compared label, reference, and citation keys, as well as figure commands: all preserved.
- Verified that inactive comment lines and source line counts are unchanged.
- Verified that all 329 recorded revised lines match the current sources and contain no trailing whitespace.
- The check results are saved in `build/pass2-reviews/inverse-checks.txt`. The final whole-book and standalone build results and visual checks appear in the validation section.

### Learning and optimization

#### Basics of Machine Learning

Reworked the introduction and transitions to explain the purpose of each model before its construction. Replaced informal language about axes of “gravity,” teleporting centroids, “fuzzy” probabilities, and the “perfect example” of kernelization with concrete descriptions of variation, reseeding, probability transitions, and feature-space computation.

Clarified the PCA feature/reconstruction interpretation; cluster assignments and centroid updates; the meaning of the k-means++ approximation guarantee; empirical versus population risk; the role of training, validation, and model complexity; conditional-mean and random-design regression language; the purpose of nearest-neighbor and logistic baselines; and multiclass likelihood estimation. Rewrote the log-sum-exp stabilization explanation, feature-matrix interpretation, kernel computational-cost discussion, and representer-theorem implications. Polished the PCA, clustering, nearest-neighbor, logistic, and digit-image captions and standardized “Probabilistic modeling.”

In `sec-pca-theory.tex`, made the proof strategy explicit: restrict to projections, derive a trace objective, bound it by a linear program, and exhibit attainment. Corrected vague antecedents, repeated “hence,” and the wording of the linear-program maximizing vector. All formulas and proof steps were preserved.

In `machine-learning-sec-pac.tex`, corrected “we select ... and minimizes,” tightened the Bayes-estimator and calibration discussion, distinguished approximation from sampling error without an informal bias/variance slogan, and clarified the progression from concentration to Rademacher complexity. Rewrote the norm-radius and kernel-trick explanations and replaced unsupported frequency claims about NP-hardness across “most” model classes with the narrower computational-difficulty statement used by the argument.

#### Smooth Optimization

Reorganized sentences around existence, coercivity, convexity, stationarity, and uniqueness so that their logical roles are explicit. Replaced calques such as “one has with the law,” “the first use of Hessian,” “re-manipulating,” and “decaying sequence.” Corrected agreement, the “minimizerss” heading, and minimum/minimizer terminology.

Improved the distinction between coordinate derivatives and differentiability, the role of the transpose as an adjoint, normal equations and pseudoinverses, the PCA/covariance connection, and the chain-rule explanation. Made the descent-direction and line-search discussion more direct, clarified geometric versus linear convergence terminology, and tightened the progression from quadratic spectral analysis to uniform Hessian bounds. Standardized positive-definite versus positive-semidefinite terminology and clarified the conditions under which a Newton step has a positive-definite metric.

#### Advanced Optimization

In `sec-regul.tex`, clarified the statistical role of regularization, the vanishing-regularization proof, primal/dual ridge systems, the cost of solving them, and the purpose of sparse feature selection. Improved the soft-thresholding proof and the surrogate-minimization derivation. The quadratic correction is now called nonnegative, consistent with the stated non-strict spectral inequality.

In `sec-stochastic-optim.tex`, rewrote the motivation around the cost of an update rather than blanket superiority claims. Clarified unbiased gradients, random-iterate convergence, computationally comparable plots, step-size schedules, averaging, and SAG's finite gradient table. Repaired sentence fragments across source-line boundaries, removed an unmatched quotation mark around “classical” SGD, and replaced an ambiguous “initial” article. The averaging discussion now describes reduced fluctuations and possible reduced tuning sensitivity without promising universal acceleration.

In `sec-autodiff.tex`, rewrote the opening motivation, finite-difference comparison, program/graph interpretation, local derivative notation, computational cost, and dual-number explanation. Tightened reverse-mode and backpropagation language, matrix-chain comparisons, trainable-layer derivatives, and the adjoint-state connection. Clarified recurrent weight sharing, reversible computation, equilibrium versus finite-iteration differentiation, argmin layers, and the envelope theorem. Preserved the verified Sinkhorn formulas and their gauge discussion. Standardized “nonlinear” and “nonlinearity.”

#### Shallow Learning / Perceptrons

Improved the introduction to arbitrary-depth networks, scalar outputs, activations, polynomial-degree limitations, and matrix gradient calculations. Replaced vague claims of constructiveness with the actual distinction between existence, quantitative width bounds, and efficient parameter selection. Tightened the transition from uniform approximation to mean-square approximation and the explanation of Fourier-defined Barron spaces.

Clarified the measure parameterization, bounded coefficients, first variations as derivative representations, the linear minimization oracle, and Frank–Wolfe's interpretation as adding neurons. Preserved the closed-convex-hull qualifications and explicitly labeled Fourier-to-dictionary proof sketch. Corrected the statement about linear activations with biases: their input-output map is affine, and is linear when the biases vanish.

#### Deep Learning

Refined the learned-feature interpretation, feedforward parameter indexing, training/stationarity discussion, dense-layer cost, and logistic-network equivalence. Shortened figure captions and made channel, location, impulse-response, and weight-sharing descriptions consistent. Clarified the distinction between activation gradients and filter gradients and the restriction to trainable filter support. Rewrote residual/bottleneck prose while preserving the batch-normalization, attention, and scattering content and citations.

#### Convex Analysis

Tightened definitions and the direct-method explanation. Replaced repeated informal “slope,” “size,” and “tricky problem” descriptions with supporting affine functions, subdifferential geometry, and domain qualifications. Clarified calculus rules and the normal cone at points of an affine constraint set.

Reworked the motivation for conjugacy, biconjugacy, smoothness duality, and Moreau regularization. Made the dependence of duality on the primal formulation explicit. Simplified the Lagrange/Fenchel derivations' prose, distinguished finite objectives from domain constraints, and clarified Slater qualification, primal/dual feasibility, stationarity, and complementary slackness. Formulas and qualification conditions were preserved.

#### Nonsmooth Optimization

Clarified the purpose and cost of subgradient, projected-gradient, barrier, proximal-point, forward–backward, Douglas–Rachford, ADMM, and primal-dual methods. Replaced loose claims about projections and “impossible” proximal computation with the practical issue of tractability. Improved the separation between majorization/touching and the additional conditions used for convergence.

Polished the soft-thresholding proof, resolvent and Moreau interpretations, product-space constructions, graph projection, and matrix-inversion discussion. Rewrote the explicit/implicit comparison to state first-order agreement, the stability distinction, and the computational tradeoff. Removed “Of course,” “we are lucky,” and similar qualifiers that obscured these points.

#### Mathematical and hypothesis clarifications

No inline or displayed mathematical fragment was changed. The following corrections concern the mathematical interpretation in the surrounding prose:

1. **Affine networks:** composition of affine layers with identity activation remains affine; it is linear only when the biases vanish. This follows directly from the displayed layer recurrence.
2. **Stationarity:** a zero gradient alone does not establish local minimality or convergence of a gradient algorithm. The chapter already contains counterexamples; the prose no longer promises that gradient methods generally reach local minima without further assumptions.
3. **Hessian terminology:** convexity gives a positive-semidefinite Hessian and nonnegative scalar second derivative. Positive definiteness is the stronger condition needed for the stated strict-minimum and descent-metric conclusions. A Newton metric is described under that positive-definiteness condition.
4. **Surrogate correction:** the quadratic term in the ISTA surrogate is nonnegative under the displayed non-strict spectral bound; it need not be strictly positive.
5. **Normal cones:** the normal cone of an affine set equals its orthogonal direction space at points belonging to that affine set. Outside it, the preceding empty-subdifferential convention still applies.
6. **Condition-number notation:** the prose now explicitly identifies the later epsilon notation with the already introduced inverse condition number. The preserved asymptotic formulas therefore have a defined interpretation.
7. **Logistic sharpening:** increasing parameter magnitude along a fixed direction sharpens the probability transition. This is no longer presented as an equivalence to SVM optimization.
8. **Rate and performance language:** initialization guarantees, nonparametric complexity, stochastic cost advantages, and averaging benefits are stated at the level justified by the existing formulas and assumptions; unsupported broad performance claims were narrowed.

#### Representative before/after examples

1. **PCA proof strategy** (`sec-pca-theory.tex`):
   - Before: “The next lemma provides an upper bound on the quantity being maximized as the solution of a convex optimization problem (a linear program). The proof of the theorem follows by showing that this upper bound is actually reached, which provides a certificate of optimality.”
   - After: “The next lemma bounds the objective by a linear program. We then show that PCA attains this bound, proving optimality.”
2. **Convex duality** (`convex-analysis.tex`):
   - Before: “Duality is associated to a particular formulation of the optimization problem, so that for instance making change of variables results in a different duality.”
   - After: “A dual problem depends on the chosen primal formulation. A change of variables or a different representation of the constraints can produce a different dual problem.”
3. **Second-order methods** (`sec-optim-smooth.tex`):
   - Before: “A second use, is to be used in practice to define second-order methods (such as Newton's algorithm), which converge faster than gradient descent, but are more costly.”
   - After: “The Hessian also yields second-order methods such as Newton's method. These can converge faster locally than gradient descent, at a greater cost per iteration.”
4. **Averaging** (`sec-stochastic-optim.tex`):
   - Before: “To improve somehow the convergence speed, it is possible to average the past iterates, i.e. run a ‘classical’ SGD on auxiliary variables.”
   - After: “Averaging reduces fluctuations in the estimate. Run SGD on auxiliary iterates.” The original variable notation is retained in the source.
5. **Proximal computation** (`optim-nonsmooth.tex`):
   - Before: “It is in general impossible to compute [the proximal map] so that the proximal point algorithm is not implementable.”
   - After: “For a general objective, evaluating [the proximal map] may be as difficult as solving the original problem.” The bracketed phrase substitutes for the identical mathematical notation in both versions.

### Final editorial and visual checks

- Standardized the active spelling of “nonlinear,” “nonlinearity,” and “nonsmooth,” including the approximation chapter/section titles, image-panel labels, and Gaussian-kernel classification caption. Standardized “Probabilistic modeling.” Historical asset filenames and inactive drafts retain their original spellings.
- Distinguished functions defined on a continuous domain from continuous functions in the Fourier and approximation introductions. Described the wavelet refinement equation using a function on the real line and a discrete measure; scaling functions such as Haar need not be continuous.
- Shortened the spherical-harmonics opening to remove the only overfull line introduced by paragraph reflow. Clarified the eigenfunction/eigenvalue sentences and used present tense for spherical coordinates. Removed a leading space before a tab in the Fourier theorem statement.
- Corrected the first panel label of Figure 4.22 from `a_j` to `a_{j-1}`, matching both the image's internal annotation and the input of the displayed forward wavelet step.
- Moved the horizontal-axis label below the lower Gaussian-maximum plot in Figure 7.16. The original empirical image and plotting commands are preserved in `figures/denoising/max-gaussian-std-corrected.pdf`. The original asset is retained. Optional regeneration uses `scripts/correct_max_gaussian_std.py` with `pypdf`; this adds no dependency to the normal LaTeX build. The README documents the distinction.
- Checked all active displayed equation blocks, labels, and citation keys against the saved pre-pass sources. Reviewed newly introduced spelling flags; they are mathematical terminology or proper names. The final source diff passes whitespace checks.
- Rebuilt the book and all standalone chapters, audited PDF links and embedded fonts, compared chapter numbering, and inspected every book page plus the first and last pages of every standalone chapter. Enlarged the corrected figure pages and dense mathematical pages. Results are recorded below.

## Typography and presentation

- Created `book-preamble.tex` as the shared design for the book and standalone chapters. Replaced the mixed legacy typography with Palatino-style New PX text and mathematics, paired with Source Sans Pro headings.
- Set consistent A4 margins, slightly increased line spacing, microtypographic spacing, restrained navy/teal colors, clearer chapter openings, readable captions, and unobtrusive running headers and page numbers.
- Redesigned the title page around the existing wave illustration, with a clear title hierarchy, author information, and working project links.
- Gave the preface a contents entry and rewrote it to match the material actually included. Added the references to the contents and standardized their heading as “References.”
- Removed blank verso pages from the digital edition and corrected the cover's paragraph indentation.
- Unified theorem, proposition, lemma, corollary, definition, procedure, remark, and example numbering within chapters. Definitions use upright text; theorem-like statements retain their appropriate styles. The shared style also supports documents without a chapter counter.
- Corrected norm delimiters and their spacing, preserved the book's transpose notation when loading the new mathematics font, and restored differential, theta, and sphere commands after font/language initialization.
- Removed an unused notation package that overwrote active mathematical and accent commands. Removed the duplicate listing-float declaration and standardized code-listing typography.
- Declared the scalable script/symbol font sizes used by the new mathematics font, eliminating font-substitution warnings.
- Replaced fragile wrapped illustrations by balanced floating illustrations with side captions. Adjusted overwide image grids, split long formulas and derivations where necessary, and corrected mathematical heading text used in PDF bookmarks.
- Repaired duplicate and missing labels, broken references, empty citations, and links to inactive material. Unique PDF destinations now also work for manually tagged equations.

## Complete and standalone builds

- Added `scripts/build_book.py` and a `Makefile`. `make` builds the complete book and all active standalone chapters; the include list is the single source of truth.
- Kept the requested output paths: `FundationsDataScience.pdf` and `chapters-pdf/<chapter>.pdf`.
- Added a standalone mode to the main document. Each chapter preserves its book chapter, equation, and theorem numbers, while its page numbers begin at 1.
- Each standalone chapter has a local bibliography containing its own citations. Chapters without citations omit an empty bibliography.
- Cross-chapter references in standalone PDFs link to `../FundationsDataScience.pdf`, with the numbers and destinations imported from the current complete build. Only active chapters' labels are imported; stale auxiliary files from excluded chapters are ignored.
- Isolated compilation in `build/`. BibTeX uses the current build's chapter auxiliaries, avoiding stale source-folder bibliography data. Existing source-folder auxiliary files are neither deleted nor required by this workflow.
- Added repeated LaTeX passes until references stabilize, explicit compilation failures, and final diagnostic checks. Warnings, missing characters, and overfull or underfull boxes prevent publication of PDFs by default; no warning-suppression package or blanket tolerance setting is used to hide them.
- Added `build/build-report.json` and retained detailed per-document logs. A development-only `--allow-warnings` option permits intermediate previews.
- Replaced the obsolete MATLAB chapter-building script with a wrapper around the same Python workflow. Documented installation, output layout, commands, and the need for a full build before distributing a synchronized set of PDFs in `README.md`.
- Converted the formerly Latin-1 active chapter sources to UTF-8 without dropping accented characters. Added build/cache exclusions to `.gitignore`.

## Repaired figures and bibliography

- Replaced the incorrect classification-loss graphic with an editable vector plot. It now uses the negative margin `s = −y f(x)`, hinge loss `max(0,1+s)`, logistic loss `log(1+exp(s))`, and a clearly marked zero–one loss with ties counted as errors. The old hinge curve was shifted incorrectly.
- Replaced the legacy variance-stabilization figure containing the visible note “Formula is false.” The new vector plot computes the variance of the square-root, Anscombe, and Freeman–Tukey transforms directly from Poisson probabilities for means from 0 to 20. The probability sums are checked, and the unit-variance reference is explicit.
- Replaced the wavelet-support sketch with a vector diagram whose horizontal coordinate is `2^j n_1` and vertical coordinate is `2^j n_2`; the original labels were swapped.
- Added editable figure sources and `scripts/variance_stabilization.py`. The build regenerates the corrected plot PDFs when their source or data changes.
- Removed unknown SNR/acquisition placeholders and corrected swapped or mismatched captions as itemized in the chapter records. Missing numerical acquisition parameters were not guessed. Most historical experimental graphics were retained.
- Corrected the Foucart–Rauhut book entry by removing spurious volume/number fields that caused a BibTeX warning; corrected Ciarlet's book entry type and publisher.
- Updated the existing Bach book entry to the published 2024 MIT Press edition, retaining its citation key. Verified against [MIT Press](https://mitpress.mit.edu/9780262381369/learning-theory-from-first-principles/).
- Added the missing Duarte et al. single-pixel-camera article information, including IEEE Signal Processing Magazine 25(2), 83–91 (2008); checked the [author's publication record](https://mdav.ece.gatech.edu/publications/).
- Added the original batch-normalization and transformer references, verified from [PMLR](https://proceedings.mlr.press/v37/ioffe15.html) and [NeurIPS](https://papers.neurips.cc/paper_files/paper/2017/hash/3f5ee243547dee91fbd053c1c4a845aa-Abstract.html).
- Added Tropp's *Just Relax* and Candès–Tao's *Decoding by Linear Programming*. The ERC claim is proved directly, with awareness of the corrigendum listed in the [Tropp publication record](https://authors.library.caltech.edu/records/v3yk2-6z279).
- Corrected the year and journal details of Candès's RIP article to Comptes Rendus Mathématique 346(9–10), 2008, verified against [Numdam](https://archive.numdam.org/articles/10.1016/j.crma.2008.03.014/).

## Additional cross-chapter review

- Cross-checked Fourier normalization, finite-difference signs, adjoints, noise scaling, and source-exponent conventions across the signal, inverse-problem, and learning chapters.
- Corrected both squared Tikhonov bias denominators to `sigma_m^(2 beta)`; the exponent must agree with the squared source-condition sum.
- Qualified Barron approximation with the continuous Fourier-inversion representative, so that the result also applies to sampling measures with atoms.
- Qualified Douglas–Rachford convergence when the minimizer is nonunique, used infima in the intermediate Fenchel derivation, and specified block shrinkage for isotropic TV.
- Performed a further English pass after mathematical changes, including the compressed-sensing acquisition model, random-matrix explanations, quantifier order, and RIP proof transitions. Removed duplicate explanations and repaired introduced typographical errors before the final builds.
- In the compressed-sensing review, explicitly treated zero coherence without division by zero, defined the support-sign vector in the coherence proof, specified symmetric normalized Bernoulli entries, distinguished the minimum-norm data vector from its certificate, and made the restricted-orthogonality constant the smallest admissible bound.
- The final visual audit found equation labels attached to unnumbered displays, which silently inherited unrelated numbers despite compiling without warnings. Those displays were numbered explicitly and their uses checked.
- Refined the remaining deep-learning exposition on shared-parameter differentiation, channel dimensions, periodic convolution, downsampling, receptive fields, and output layers.

## Final validation and delivered PDFs

The final full `make` completed successfully after the second editorial pass, with no development warning override. All 20 documents passed the final LaTeX/BibTeX diagnostic checks: no unresolved references or citations, duplicate destinations, font warnings, missing characters, or overfull/underfull boxes. References stabilized in each complete and standalone build.

The complete book has **248 pages**, compared with 252 before this second pass. All chapters and displayed equations remain present. The PDF audit checked **1815 internal links** and **156 links from chapters back to the book**; every named destination exists. All font entries are embedded. Text extraction found no unresolved-reference or editorial TODO markers in the delivered PDFs.

The source audit verified **1276 unchanged active displayed equation blocks**, together with preserved label definitions and citation keys. Inline dimensional, caption, and explanatory corrections are recorded above. The auxiliary-file audit compared **725 labels** and found that every chapter label has the same number in the complete book and its standalone PDF.

Every book page was rendered and reviewed for layout. Dense mathematical pages and the corrected figure labels received enlarged checks. The first and last pages of every standalone chapter were also inspected. The two figure issues found during the review were corrected and the affected pages checked again.

| Chapter | Source / standalone PDF | PDF pages | Final diagnostics |
|---|---|---:|---:|
| 1 | `chapters/shannon.tex` → `chapters-pdf/shannon.pdf` | 8 | 0 |
| 2 | `chapters/fourier.tex` → `chapters-pdf/fourier.pdf` | 21 | 0 |
| 3 | `chapters/shannon-coding.tex` → `chapters-pdf/shannon-coding.pdf` | 5 | 0 |
| 4 | `chapters/wavelets.tex` → `chapters-pdf/wavelets.pdf` | 21 | 0 |
| 5 | `chapters/approximation.tex` → `chapters-pdf/approximation.pdf` | 19 | 0 |
| 6 | `chapters/compression.tex` → `chapters-pdf/compression.pdf` | 7 | 0 |
| 7 | `chapters/denoising.tex` → `chapters-pdf/denoising.pdf` | 22 | 0 |
| 8 | `chapters/variational-priors.tex` → `chapters-pdf/variational-priors.pdf` | 10 | 0 |
| 9 | `chapters/inverse-problems.tex` → `chapters-pdf/inverse-problems.pdf` | 18 | 0 |
| 10 | `chapters/sparse-regularization.tex` → `chapters-pdf/sparse-regularization.pdf` | 11 | 0 |
| 11 | `chapters/sparse-theory.tex` → `chapters-pdf/sparse-theory.pdf` | 15 | 0 |
| 12 | `chapters/compressed-sensing.tex` → `chapters-pdf/compressed-sensing.pdf` | 11 | 0 |
| 13 | `chapters/machine-learning.tex` → `chapters-pdf/machine-learning.pdf` | 22 | 0 |
| 14 | `chapters/optim-ml-smooth.tex` → `chapters-pdf/optim-ml-smooth.pdf` | 15 | 0 |
| 15 | `chapters/optim-ml-advanced.tex` → `chapters-pdf/optim-ml-advanced.pdf` | 15 | 0 |
| 16 | `chapters/perceptrons.tex` → `chapters-pdf/perceptrons.pdf` | 9 | 0 |
| 17 | `chapters/deep-learning.tex` → `chapters-pdf/deep-learning.pdf` | 5 | 0 |
| 18 | `chapters/convex-analysis.tex` → `chapters-pdf/convex-analysis.pdf` | 9 | 0 |
| 19 | `chapters/optim-nonsmooth.tex` → `chapters-pdf/optim-nonsmooth.pdf` | 13 | 0 |

Build details are recorded in `build/build-report.json`; link/font checks are in `build/pdf-audit.json`; source preservation checks are in `build/pass2-source-audit.json`; number comparisons are in `build/label-number-audit.json`; rendered review sheets are in `build/qa/pass2/`. These are generated validation artifacts. Run `make` to regenerate a synchronized set of PDFs after further edits.

## Chapters 1–7: sampling, Fourier analysis, coding, wavelets, approximation, compression, denoising

### Reviewed coverage and validation

Read all active prose, displayed equations, theorem statements, proofs, algorithms, and captions in `chapters/shannon.tex`, `chapters/fourier.tex`, `chapters/shannon-coding.tex`, `chapters/wavelets.tex`, `chapters/approximation.tex`, `chapters/compression.tex`, and `chapters/denoising.tex`. No additional chapter input files occur in these seven files. Inactive `\if 0` drafts and commented MATLAB listings were not treated as finished book content. Both Shannon chapters received an additional prose and mathematical-consistency pass, alongside the visual and build checks.

All seven chapters received sentence-level editing: agreement, articles, punctuation, terminology, captions, transitions, and clearer separation of assumptions, conclusions, and experimental observations. The entries below identify the substantive changes and grouped writing corrections by section. Source section names, equation labels, and theorem labels are supplied because line numbers move as the book is edited.

Independent numerical checks passed for the repaired FFT recursion, first/second finite-difference convolution and Fourier symbols, Haar analysis/synthesis and energy conservation, odd-length spectral zero-padding interpolation, dead-zone quantization error bounds, and firm-threshold endpoint continuity. These checks used direct sums and the displayed formulas rather than a library implementation of the repaired algorithms. Compilation, reference, typography, and rendering results are recorded in the final validation section above.

### `chapters/shannon.tex`

- **Opening and analog/discrete signals:** rewrote the three-topic introduction, removed repeated “it studies,” corrected the channel-coding terminology, clarified continuous models versus sampled data, and corrected dimensional examples and captions. A second prose pass refined the final text.
- **Linear sampler (`eq-linear-sampling`):** replaced the incorrectly fixed integration interval by convolution over the extended signal domain. Explained that the kernel width controls spatial resolution and that exact reconstruction needs bandlimiting rather than smoothness alone.
- **Fourier reminders:** corrected Plancherel from the reciprocal constant to `||hat f||² = 2π ||f||²` for the stated unnormalized transform. Specified the `L²` meaning of truncated-transform convergence. Replaced insufficient differentiation/decay assumptions by `W^{p,1}` hypotheses. Added the missing `1/(2π)` in the Fourier expression for the squared Sobolev norm and clarified weak derivatives.
- **Fourier series (`eq-fourier-series`):** repaired the inconsistent summation index, the coefficient-decay index, and the prose on normal convergence. Restricted the jump-convergence assertion to piecewise regular functions, for which the one-sided midpoint conclusion is justified.
- **Poisson formula:** removed the duplicate `eq-poisson-formula` label from its proof, clarified the Fourier-sign convention, and explicitly selected continuous representatives. The continuity hypothesis prevents arbitrary modification of sampling values on null sets.
- **Sampling theorem (`thm-shannon-sampling`):** corrected the statement that the *Fourier transform* becomes analytic; it is the bandlimited signal that has an entire extension. Defined sinc at zero, corrected the missing factor `2π` in the absolute-integrability calculation, and supplied the uniform-convergence argument. Removed its accidentally duplicated sentence during final checking.
- **Oversampling and splines:** replaced the unjustified assertion of exponential decay for a smooth compactly supported frequency window by faster-than-any-power decay. Clarified cardinal B-spline definitions, interpolation coefficients, and the cubic example.
- **Aliasing and finite storage:** assigned distinct labels to the spline and aliasing figures. Clarified that a sinusoid is a spectral-line example rather than an `L²(R)` example. Replaced the blanket “finite storage makes perfect reconstruction impossible” claim by nonidentifiability for arbitrary bandlimited signals from finitely many samples; finite-dimensional model classes are different.
- **Quantization:** distinguished integer quantizer indices from reconstruction levels, stated the worst-case midpoint property accurately, clarified finite-level clipping/endpoints, and described quantization loss relative to the already sampled vector. Repaired the missing quantizer reference and figure wording.
- **References:** sampling and source coding no longer share the same chapter label. Added/retained the sampling references needed by the active book.

### `chapters/shannon-coding.tex`

- **Coding definitions:** replaced infinite-sequence notation for codewords by finite binary words, repaired the prefix-tree definition and tree traversal prose, and made the singleton-alphabet empty-word convention explicit. A message length must then be known separately.
- **Probability model:** consistently names the symbols, distinguishes an empirical histogram from an assumed source law, and notes that an estimated coding model must be available to the decoder. Corrected grammar and the entropy illustrations.
- **Entropy:** kept logarithms consistently in base two; gave correct derivatives of `-u log₂u`; distinguished differential entropy from discrete entropy and noted its dependence on reference measure and possible negativity.
- **Entropy bounds:** rewrote the proof using `ln u ≤ u−1`, correctly restricting sums to positive-probability symbols. Corrected the erroneous upper-bound sign (`log(1/K)` became `log₂K`) and stated the equality cases.
- **Relative entropy:** corrected “Kulback” to “Kullback–Leibler,” “triangular” to “triangle,” and “positive” to “nonnegative.” Clarified zero-probability conventions and distinguished ideal cross-entropy penalties from integer-code rounding effects.
- **Kraft inequality:** made integer nonnegative lengths an explicit hypothesis; corrected subtree *height* from `2^{m−ell}` to `m−ell`, while retaining `2^{m−ell}` as the number of leaves. Explained the dyadic alignment argument in the constructive converse. Renamed the Kraft illustration so it no longer collides with the Fourier/wavelet-atom illustration.
- **Shannon source-coding theorem:** replaced the incomplete constrained-optimization proof by the relative-entropy proof using `q_k=2^{-ell_k}/Z`. Removed zero-probability singularities by working on the support. Corrected the upper bound to the strict `L < H+1`, including the singleton case.
- **Algorithms:** called the rounded-length construction a Shannon code; described Huffman as a greedy merging algorithm rather than dynamic programming, with `O(K log K)` complexity using a priority queue.
- **Block coding:** replaced the unnecessary large-sample argument by an exact expected-length statement under an i.i.d. model. Specified equal marginals for `H(P) ≤ rH(p)`, and corrected the dependent-case equality to an inequality. Distinguished block prefix codes from symbol-by-symbol prefix codes.
- **Dependence:** explained that reversible differencing preserves joint entropy while changing marginal statistics. Removed the unsupported claim that run-length coding automatically reaches the entropy rate of every Markov source. Gave the stationary Markov entropy-rate formula and described conditional/block coding as the appropriate route to that rate.
- **Presentation:** corrected the entropy figure’s column wrapper and malformed prose/escape introduced during intermediate editing; the coordinating review checked both Shannon chapters again.

### `chapters/fourier.tex`

- **Hilbert spaces:** made the inner-product convention explicit, replaced “orthogonal” by “orthonormal” where unit norms are needed, corrected Parseval to use absolute squares, and repaired the Gram–Schmidt coefficient from `<bar phi_i,phi_i>` to `<bar phi_k,phi_i>`. Stated independence/dense-span assumptions.
- **Polynomial examples:** distinguished conventional Legendre normalization from unit `L²` normalization. Replaced the inconsistent Gaussian/Hermite formulas by the probabilists’ Hermite convention with weight `exp(−x²/2)`, and explained the corresponding Hermite functions for Lebesgue measure.
- **Fourier series (`eq-defn-fourier-coeffs`):** repaired the malformed `hat f_kk` and the summation index. Clarified the distinction between a countable circle basis and the continuous real-line transform.
- **Convolution:** stated normalized Haar measure on the circle, so the product formula has the claimed constants. Added the exponent range in Young’s inequality, restricted the continuity statement to a justified case, and added bounded-derivative assumptions for differentiating a convolution. Specified an approximate identity, including the circle’s `2π` normalization. Replaced the incorrect `ell¹` claim for Fourier coefficients of arbitrary `L¹` functions by bounded coefficient sequences.
- **Translation-invariant operators:** repaired the domain in the continuous-output proposition, used `C_b` for the supremum norm, restored missing equalities and the integration variable, and correctly conjugated the Riesz kernel. Removed the unlicensed application of `H` to a Dirac mass from the rigorous argument. Added the converse for bounded Fourier multipliers.
- **Distributions:** distinguished multiplication of measures by continuous functions from multiplication of general distributions by smooth functions. Corrected the transform-duality sign for bilinear distribution pairing. Repaired the Dirichlet-kernel variable and clarified its distributional, rather than pointwise, limit.
- **DFT and FFT:** corrected zero-based inner-product indices, the FFT twiddle denominator (`N`, not `N/2`), and complex-valued intermediate spaces. Explained that zero-padding changes the transform. Corrected the target algebra to complex vectors and removed a stray equality in convolution diagonalization.
- **Finite convolution and polynomial products:** restored missing equalities, fixed translation notation, repaired polynomial powers/coefficient indices, and increased the zero-padding length to `A+B+1` to avoid wraparound.
- **Spatial/spectral zero-padding:** clarified fixed-frequency quadrature accuracy and Nyquist versus spectral-grid spacing. Rebuilt the warned-about spectral interpolation derivation with inverse factor `1/Q`, negative forward-transform sign, correct `Q/N` rescaling, signed-frequency placement modulo `Q`, and removable-singularity conventions. It exactly evaluates the trigonometric interpolant.
- **Multiple dimensions:** corrected the inverse-transform integrand to `hat f(omega)`, restored the missing `2π` phases in the multidimensional DFT, and distinguished an index pairing from a signal inner product. Corrected tensor indices/conjugates and the 2-D matrix formula to `F1 f F2^T` rather than `F2*`.
- **PDE discretization:** corrected `R/NZ` to `Z/NZ`, the sign and normalization of the forward-difference filter, `D2=−D1^T D1`, and the missing frequency index in its symbol. Clarified periodization of the heat kernel.
- **Group theory:** repaired the missing exponential, product group order, unit-circle requirement for infinite-group duals, and `2π` in the distributional character pairing. Corrected the symmetric-group order description, transposition parity, and the distinction between one-dimensional homomorphisms and traces of representations. Restricted invariant-complement arguments to finite groups, used complex spaces and unitary representatives, repaired representation equivalence, and described the nonabelian dual as a set rather than a group. Restored `|G|` in the regular-character identity and corrected its use in inversion. Convolution becomes block-diagonal, not scalar-diagonal.
- **Manifold Laplacian:** replaced the zero-valued averaging limit by the correctly normalized `2(d+2)/epsilon²` mean-value formula. Replaced the false “Laplacian is compact” claim by an unbounded self-adjoint operator with compact resolvent. Added compact/no-boundary hypotheses, multiplicities, the constant zero mode, and the correct unit-torus `4π²` eigenvalue factor.
- **Spherical harmonics:** repaired spherical coordinates and `S²` versus `S³`, described associated Legendre functions accurately, and stated normalization requirements.
- **Graph Laplacian:** required nonnegative symmetric weights, corrected the cycle-graph scaling, rebuilt the quadratic-form proof, and corrected the connected-graph conclusion to a simple zero eigenvalue and strictly negative second eigenvalue. Removed references to the excluded meshes chapter and supplied the Dirichlet-energy interpretation locally.
- **Presentation:** added safe bookmark text for mathematical headings, split a long derivation, adjusted overwide figure rows, and polished prose throughout without renaming existing misspelled asset paths.

### `chapters/wavelets.tex`

- **Multiresolution:** repaired nested-space indices, complete/orthonormal terminology, Shannon sample coefficients (`2^{j/2} f(2^j n)`), and the corresponding reconstruction normalization.
- **Haar illustration:** replaced the graphic whose caption explicitly reported an incorrect amplitude by native formulas for scaling and wavelet functions, both with amplitude `2^{−j/2}` and unit norm.
- **Spectral orthogonalization:** corrected the reversed Kronecker-delta cases, reflected-conjugate kernel, missing Fourier hat in Poisson summation, and closed-span formulation.
- **Detail spaces:** corrected the direct-sum identity to use `W_j`, reversed the wrongly directed truncated convergence limit, and repaired the Haar detail intervals and mother-wavelet normalization.
- **Bounded domains and sampling:** specified integer coarse levels on the unit torus, fixed the translation range, inserted the sample-normalization factor, and stated the role of `integral phi=1`. Added chapter/section labels needed by other active chapters.
- **Forward FWT:** repaired level ranges and output sets, input to the elementary step, vector dimensions, downsampling endpoints, refinement-space membership, and the wavelet expansion in *scaling* functions. Standardized reflected-filter notation and explicitly treated real filters. Corrected finite-filter support assumptions and the geometric-series operation count.
- **Haar algorithm:** fixed the high-pass filter index and the use of level `j`, rather than `j−1`, in the update to `j+1`.
- **Inverse FWT:** gave descending level order, fixed upsampling endpoints, repaired the second input of the paired convolution operator and the adjoint `S2*`, and separated forward/inverse equation labels.
- **2-D bases:** restored scaling and mixed scaling–wavelet products in the anisotropic periodic basis; wavelet–wavelet products alone are incomplete. Corrected support scales, independent tensor indices, the three distinct detail spaces, two-dimensional translation indices, and approximation-coefficient notation.
- **2-D algorithms:** replaced the erroneous second low-pass filter by a high-pass filter, fixed level ranges, and reversed the incorrect synthesis order to *upsample then convolve*. Repaired vertical-index definitions and horizontal/vertical mixups in the inverse algorithm.
- **Filter design:** added finite-support/summability assumptions for the stated converse construction; fixed the infinite product to frequencies tending to zero (`omega/2^k`, `k≥1`). Corrected QMF hypotheses to require the energy-complement condition, separated time/Fourier labels, and removed the false assertion that all high-pass phase choices yield the same wavelet.
- **Moments and support:** restricted moment indices to nonnegative integers, used integrability of moments rather than unspecified regularity, restored Fourier hats in derivative statements, fixed the constant in the auxiliary factor, and corrected the repeated low-pass condition reference. Defined support length as the smallest enclosing interval, qualified the symmetry statement with the Haar exception, and softened the claim that smoothness is merely cosmetic.
- **Presentation:** polished algorithm descriptions, terminology, and captions; retained original asset filenames; removed the dangling discrete-wavelet display reference.

### `chapters/approximation.tex`

- **Basis and thresholding:** made orthonormality explicit, corrected the coarse-wavelet indexing description, and replaced the false bijection between threshold and retained count by the correct interval/tie convention. Standardized ordered magnitudes to one-based indexing and supplied the missing forward implication in the tail-decay equivalence for `alpha>0`.
- **Smoothness models:** defined integer and noninteger Hölder regularity correctly, included mixed derivatives in several dimensions, made periodic boundary conditions explicit for Fourier derivative formulas, and distinguished the homogeneous Sobolev seminorm from a full norm. The model now controls both `L²` energy and the Sobolev seminorm.
- **Piecewise regularity:** fixed the partition endpoints and restriction notation, made the edge set a union rather than a tuple, specified one-sided regularity of regions, and corrected “modem” and related prose.
- **Bounded variation (`eq-tv-coaera`):** replaced the false level-set-length formula for discontinuous images by the coarea formula using perimeters of superlevel sets. Distinguished relative-domain perimeter from whole-space perimeter and kept the indicator-function example with the correct interpretation.
- **Cartoon model:** added uniform curve-length control and clarified arc-length parametrization, blur, and edge-localization assumptions.
- **Fourier approximation:** stated integer regularity for the elementary integration-by-parts bound, repaired the missing negative-frequency tail in the Sobolev proof, used the periodic domain, and distinguished uniform minimax rates from rates of individual functions. Corrected the interval-indicator formula; its Fourier magnitudes are not pointwise asymptotic to `1/m`. Repaired the multidimensional Sobolev norm.
- **Wavelet coefficient decay:** stated support/extension and multidimensional moment assumptions; replaced an unjustified `||psi||_1` constant by the appropriate weighted moment integral.
- **1-D approximation proof:** made wavelet/boundary assumptions explicit, restricted scale sums to the fixed coarsest level, retained scaling coefficients, rounded cutoff scales, corrected the step numbering, and explained how to choose a threshold from a prescribed coefficient budget.
- **2-D approximation proof:** added `alpha≥1` for the asserted `M^{-1}` rate and corrected the regular-coefficient count from exponent `−1/(alpha+1)` to `−2/(alpha+1)`. Stated the different rate for `0<alpha<1`.
- **Adaptive triangulation:** corrected near-edge dimensions to length `M^{-1}` and width `M^{-2}` and fixed the corresponding aspect-ratio caption. Removed the outdated unsupported claim that no algorithm can attain the rate.
- **Curvelets:** corrected the rotation convention, separated schematic scaling/sampling from exact frame construction, and explicitly included the required frequency windows and low-frequency family in the Parseval identities. An arbitrary mother function and the displayed sampling pattern do not by themselves define a tight frame. Corrected coefficient/reconstruction captions and qualified claims about compression suitability.
- **References and presentation:** repaired inactive duplicate labels and active missing references, added the chapter alias, renamed the hard-threshold label to avoid collision with sparse regularization, narrowed an overwide figure row, and preserved the actual `curvelet-spacial` asset path. Converted the original Latin-1 source to UTF-8.

### `chapters/compression.tex`

- **Quantizer:** corrected the dead-zone interval to `(-T,T)`, distinguished quantizer indices from dequantized values, defined the nonzero support, and reconstructed with `D_T(Q_T(a))` rather than the integer `Q_T(a)`. Made threshold endpoints and the zero reconstruction explicit.
- **Error bound (`prop-error-compression`):** repaired endpoint hypotheses and related the quantized approximation to the best approximation with the same retained support size.
- **Rate theorem (`thm-compression`):** rebuilt the theorem and proof using an upper `M`-term approximation bound and polynomial finite-access size in `1/T`. The old proof incorrectly treated two-sided tail asymptotics as equivalent to two-sided coefficient asymptotics. The new proof uses the scale `K=T^{-2/(alpha+1)}`, a valid bound on the sum of quantization errors, and `R=O(K log K)`.
- **Finite computation:** corrected the wavelet scale sign; explained fixed computable basis ordering versus sorting by unknown coefficient magnitude. Corrected the former self-referential “N polynomial in N” hypothesis.
- **Bit counts:** added ceiling/integer requirements, bounded coefficient values using `||f||_2`, included count/header costs, and corrected inverse `r log r` asymptotics (the previous `+o(1)` assertion was false). Stated what information the decoder must know.
- **Huffman discussion:** removed the false pointwise Huffman length inequality and wrong logarithm sign. Applied the Shannon construction to bound Huffman’s *expected* total length; an arbitrary realized code length need not obey that expectation bound.
- **JPEG-2000:** removed “latest standard”; distinguished irreversible CDF 9/7 from reversible 5/3 lifting, and qualified orthonormal error identities for the biorthogonal transform. Corrected the comparison caption to top/bottom, “Steam” to “Stream,” significance tests to use magnitudes, the refinement bit to the next binary digit, the coding scan to four-sample vertical stripes, and the neighborhood to all offsets in `{-1,0,1}²`.
- **References/layout:** repaired wavelet-transform and sampling-consistency references and used proof-compatible bold subheadings. Rewrote explanatory prose throughout.

### `chapters/denoising.tex`

- **Noise model and risk:** corrected white-noise covariance to `sigma² I`, replaced a continuous-noise point-probability formula by a density, distinguished generalized Gaussian noise from impulses, made the real-signal setting explicit, and replaced use of the denoising operator as an expectation symbol by `E`. Added the missing chapter label.
- **Filtering:** constant preservation is a filter-design condition, not a consequence of zero-mean noise. Corrected Fourier-log captions and the filter/image labels, the optimal SNR direction (maximize), and the bias/noise discussion. Increasing the number of retained coefficients decreases the denoising strength.
- **Oracle/Wiener estimators:** normalized Fourier coordinates by `sqrt N`, repaired `hat h=lambda`, restored the missing equality and summation parentheses, and gave gains without division by zero for zero signal energy. Distinguished real orthogonal Gaussian invariance from dependent conjugate Fourier coefficients. Rebuilt the averaged autocorrelation formula and its normalization. Clarified the bias/variance problems of a raw noisy periodogram.
- **Linear denoising theorem:** replaced a noninteger unconstrained choice of `M` by `min(N,ceil(t))`, stated parameter assumptions, kept the exact bias–variance identity, and proved a valid bound including rounding. Removed the incorrect claim that every finite-energy random vector must be finite-dimensional; white noise is the exceptional infinite-energy process here.
- **Discrete Sobolev interpretation:** required consistent normalization, boundary conditions, and frequency ordering rather than claiming automatic equivalence for arbitrary discrete samples.
- **Hard/soft thresholding:** corrected soft-threshold definition at zero, repaired coarse/fine level indices, and qualified image-specific empirical threshold comparisons. Corrected Donoho–Johnstone’s names and grammar.
- **Universal threshold:** replaced the claimed exact mean of the Gaussian maximum by an asymptotic statement, removed its claim to exactly minimize realized risk, restored the missing factor `sigma`, and replaced the misreferenced theorem by an oracle inequality for the actual nonlinear estimator. Included a residual-noise term and showed how approximation rates imply the stated risk rate.
- **Translation invariance:** corrected “computes” to “commutes,” repaired quantified indices, stated a sufficient basis condition instead of an unjustified equivalence, and required a finite group of shifts for exact cycle-spinning equivariance. Explicitly noted that a selected 4×4 shift subset gives approximate, not guaranteed exact, invariance.
- **Undecimated frames:** corrected 1-D versus 2-D coefficient counts (`log₂ N` versus `3 log₂ N0`), included the coarse array, explained multiplicity weights when identical translated atoms are merged, and replaced inconsistent array indexing by direct pixel indexing. Corrected “à trous.”
- **Firm/Stein/block shrinkage:** restored the missing factor `mu` in firm thresholding so endpoints join continuously; described endpoint limits correctly. Reversed the erroneous attenuation ratios to `T/|x|` and `T²/|x|²`; defined values at zero; distinguished block sets from energies and used the square root of mean energy as the attenuation input. Replaced “unbiased” implications by the correct large-amplitude shrinkage statement.
- **Poisson model:** allowed any nonnegative real mean, distinguished observed counts from latent intensities, corrected variance notation, and replaced a misleading relative-error expectation statement by relative standard deviation `1/sqrt(lambda)`. Added inverse-transform domain/bias qualifications. Removed unverified numerical count labels from the four legacy Poisson images (the original captions repeated 50 for visibly different examples); no missing acquisition parameters were invented. Replaced the variance graphic bearing “Formula is false” by an exact Poisson-expectation plot; the source/caption now use that corrected graphic.
- **Multiplicative model:** specified independent looks, fixed sample indices and density notation, explained that averaging reduces rather than exactly cancels additive noise, and gave Gamma shape/rate, density, mean, and variance. Corrected the log-transformed target from `phi(f0)` to `log f0`, stated digamma/trigamma moments, and noted exponentiation bias. Repaired the missing figure reference and swapped stabilization caption order to match the images.
- **Presentation:** polished prose, formulas, captions, and empirical claims throughout; removed references to excluded/commented code listings.

### External checks and limits

The JPEG-2000 filter and coding-pass corrections were checked against the primary [ITU-T T.803 conformance specification](https://www.itu.int/epublications/publication/itu-t-t-803-v3-2024-02-information-technology-jpeg-2000-image-coding-system-conformance-testing) and [T.812 entry-level encoder specification](https://www.itu.int/rec/dologin_pub.asp?id=T-REC-T.812-200708-I%21%21PDF-E&lang=e&type=items). The Gaussian thresholding oracle discussion was cross-checked with Donoho and Johnstone’s [original paper](https://doi.org/10.1093/biomet/81.3.425) and its [Stanford technical-report record](https://statistics.stanford.edu/technical-reports/ideal-spatial-adaptation-wavelet-shrinkage). Bibliography repairs are listed in the shared bibliography record above.

Historical SNR measurements and most legacy illustrations were retained as experiments, not recomputed numerical benchmarks. Curvelet and wavelet existence results remain cited theorems rather than fully developed proofs; their stated assumptions and the distinction between schematic constructions and exact tight frames were corrected. The inactive unfinished compressibility draft remains excluded from the book. No claim is made that the legacy MATLAB implementations were audited or updated.

### Final English refinement

A final sentence-level pass polished the active prose in Fourier, wavelets, approximation, compression, and denoising after the mathematical review. It corrected subject–verb agreement, articles, prepositions, misplaced adverbs, dangling constructions, repeated explanations, and awkward figure references. The main groups were Fourier algorithm and convolution explanations; wavelet refinement, filter-bank, and vanishing-moment discussions; approximation-model comparisons and decay interpretations; transform coding and JPEG-2000 context coding; and denoising models, filter selection, cycle spinning, and variance stabilization. Examples include “recover efficiently missing information” → “recover missing information efficiently” and “the error decay of uniformly smooth signal” → “the error decay for uniformly smooth signals.” All modified sentences were reread; inline mathematical expressions and citation/reference targets were checked unchanged. No inactive content, formulas, or figure asset paths were changed in this final prose pass. Both Shannon chapters also received an additional prose pass.

### Final equation-label and reference audit

The final rendered-page review exposed a silent numbering problem: labels inside unnumbered displays inherited unrelated section or theorem counters while compiling without warnings. Audited all 19 active chapters, the preface, and the six recursively included section files, excluding commented and disabled content. The audit checked `\eq`, `\eqm`, starred equation/alignment/gather/multline environments, bracket and double-dollar displays, labels immediately after unnumbered displays, and labeled alignment rows whose number was suppressed. It also compared equation and named structural references with their destinations in the current book auxiliary files.

- `wavelets.tex`, “Low-pass Filter Constraints”: changed the labeled filter Fourier-series definition `eq-fourier-tr-infi` to a numbered display so its later reference identifies that formula.
- `approximation.tex`, “Decay of Wavelet Coefficients”: numbered `eq-wavcoef-decay-2`, the bounded-function coefficient estimate. Its proof references previously displayed the unrelated hard-thresholding number.
- `compression.tex`, “Decoding”: numbered `eq-decompression-approx`, the complete decoded-image reconstruction formula; its later citation had inherited the subsection number.
- `optim-nonsmooth.tex`, “Gradient Descent”: replaced “Equation” and an equation-style citation of `eq-example-grad` with a section citation, since the destination is “Examples of Gradient Computation.”
- `perceptrons.tex` introduction and `sec-autodiff.tex`, “Residual recurrent networks”: replaced numbered section citations to the unnumbered “Multilayer perceptrons” paragraph with named hyperlinks to that paragraph. This avoids presenting the inherited parent subsection number as the paragraph's own number.

No other active labeled unnumbered displays or labeled suppressed equation rows were found. Formula contents and mathematical hypotheses were unchanged; the final rebuild regenerates the corrected numbers and standalone external references.

## Chapters 8–12: variational priors, inverse problems, sparsity, compressed sensing

### Reviewed scope

Read and reviewed all active text, equations, theorem statements, proofs, captions, and cross-references in `chapters/variational-priors.tex`, `chapters/inverse-problems.tex`, `chapters/sparse-regularization.tex`, `chapters/sparse-theory.tex`, and `chapters/compressed-sensing.tex` (3,700 original source lines in total). These five chapters have no further active TeX inputs. Commented-out MATLAB listings were not treated as active book content, and the MATLAB implementations and original numerical experiments have not been rerun. Existing image assets were retained. The source labels below provide stable locations even when pagination changes.

All five sources are UTF-8; the three formerly Latin-1 files were converted without discarding characters. The review included an English editing pass, a second mathematical consistency pass, and checks for introduced misspellings, duplicate labels, empty references, and editorial placeholders. Global and standalone compilation results, cross-reference validation, and final visual checks are recorded centrally in `corrections.md` and the build report.

### Variational Priors and Regularization (`chapters/variational-priors.tex`)

1. **Chapter introduction and continuous priors:** corrected subject–verb agreement, plural forms, and the description of priors as seminorms on function spaces. Added the missing `chap-variational` chapter label.
2. **Sobolev energy (`eq-dfn-sob-energy`):** replaced the misleading norm notation by the homogeneous `H^1` seminorm. Specified the integration domain, distributional gradient, vanishing on constants, and extended value `+∞` when the gradient is not square integrable.
3. **TV definition (`eq-dfn-tv-energy`):** corrected “bounded variations” to “bounded variation,” “indicator functions,” and “seminorm.” Completed the Rudin–Osher–Fatemi citation sentence.
4. **Discrete gradient:** distinguished an `n\times n` grid from its total pixel count `N=n^2`. Periodic indices now wrap modulo `n`, not the total number of pixels. Corrected the finite-difference terminology and boundary-condition explanation.
5. **Discrete divergence:** repaired the malformed space in the vector-field dimension and corrected the integration-by-parts wording. Kept and checked the consistent adjoint convention `\operatorname{div}=-\nabla^*`.
6. **Discrete Laplacian (`eq-discrete-lapl`):** stated its negative semidefiniteness and corrected the continuous scaling to `n^2\Delta f_k` at `x=k/n`, with both continuous derivatives evaluated at `x`.
7. **Laplacian Fourier symbols (`eq-discrete-lapl-fourier`):** corrected the positive continuous multiplier to `-\|\omega\|^2`, and the discrete multiplier to `-\rho_\omega^2`, where `\rho_\omega^2=4\sin^2(\pi\omega_1/n)+4\sin^2(\pi\omega_2/n)`. The original discrete expression lacked both the sign and the factor four. Stated the continuous Fourier convention explicitly.
8. **Discrete Sobolev energy (`eq-discr-sobolev`):** grouped both squared directional differences inside the sum and its factor `1/2`.
9. **Gradient flows (`eq-pde-flow`):** distinguished an ODE for a finite-dimensional image from a PDE in continuous space. Polished the first-order interpretation and time-discretization explanation.
10. **Heat-flow gradient:** replaced the incorrect `o(\|\epsilon\|^2)` quadratic remainder by the exact remainder `\tfrac12\|\nabla\epsilon\|^2`.
11. **Heat kernel (`eq-heat-gaussian`):** removed the double negative in the exponential. The original kernel grew like `e^{+\|x\|^2/(4t)}`; it now has the integrable Gaussian exponent `-\|x\|^2/(4t)`.
12. **Discrete heat iteration:** used the current iterate in every neighboring pixel and in the convolution; specified the filter weights `1-4\tau` at the center and `\tau` at the four neighbors; included positivity in the sufficient step-size condition.
13. **TV gradient and nonsmooth flow:** required nonvanishing gradients at every pixel for the classical gradient formula. Replaced the false claim that TV flow cannot be defined by the correct subgradient inclusion; smoothing is presented as one numerical option.
14. **Smoothed TV (`eq-tv-smoothed`, `eq-grad-tv-eps`):** repaired the malformed large-ε formula to `\nabla J_{\rm TV}^\epsilon(f)=-\Delta f/\epsilon+O(\epsilon^{-3})` for fixed `f`. Explained the smaller explicit time steps required as ε decreases, rather than calling the flow intrinsically unstable.
15. **Denoising flows:** standardized the Gaussian noise variance as `\sigma^2`, corrected “initial image,” stopping-time wording, and plural captions, and restricted convergence to a constant to the connected periodic-grid setting.
16. **Regularization (`eq-regul-denoising`):** removed the contradictory unconditional claim of convexity. Distinguished uniqueness for proper lower semicontinuous convex priors from nonconvex penalties, whose existence and uniqueness require additional assumptions.
17. **Regularized gradient descent (`eq-grad-desc-regul`):** replaced the missing-data and wrong-sign update by `f^{k+1}=f^k-\tau[f^k-f+\lambda\nabla J(f^k)]`, and closed its delimiter.
18. **Sobolev update:** replaced the undefined `\Delta J(f^k` expression by `(1-\tau)f^k+\tau f+\tau\lambda\Delta f^k`. Restored the `1/2` convention in the quadratic energy and described conjugate gradients as a linear-system solver.
19. **Sobolev Fourier denoising (`eq-sobol-reg-fourier`):** made its denominator consistent with the corrected Laplacian symbol. Corrected the heat-kernel variance from `t^2` to `2t` in each coordinate.
20. **Smoothed-TV denoising (`eq-grad-tveps-regul`):** replaced continuous-time `f_t` by `f^{(k)}` in the iteration and supplied the sufficient step bound `0<\tau<2/(1+8\lambda/\epsilon)`. Clarified that it minimizes the smoothed TV objective.
21. **Figures:** exchanged the two SNR image assets whose Sobolev/TV regularization column captions were reversed. Reduced widths in the two four-column flow figures that produced overfull boxes.

### Inverse Problems (`chapters/inverse-problems.tex`)

1. **Chapter structure and English:** added `chap-invpbm`; polished acquisition, denoising, deblurring, regression, and algorithm descriptions; corrected spelling, agreement, and mathematical terminology throughout.
2. **Forward models:** separated noninjectivity, inconsistent data, and an unbounded or ill-conditioned inverse. Defined the inpainting region consistently as a spatial set or a subset of pixel indices.
3. **Medical imaging examples:** corrected “Medical resonance” to “Magnetic resonance.” Separated the duplicated inpainting/Fourier equation labels by naming the latter `eq-partial-fourier-operator`. Replaced an overly crude inverse-square scalar MEG/EEG model with the general linearized kernel formulation, whose geometry and lack of translation invariance are explicit.
4. **Regression comparison:** corrected the chapter/section reference and the dimensional correspondence. Removed categorical assertions that statistical models must be misspecified or that design matrices must be random. Qualified the `n^{-1/2}` empirical-fluctuation statement by independence, fixed dimension, and suitable moment assumptions.
5. **Oracle inversion:** defined the signal coefficients and noise variance; restricted ordinary white-noise vectors to finite dimension. Added the nonzero-coefficient condition for the zero-noise oracle inverse limit.
6. **Oracle polynomial example:** corrected the filter from the erroneous decaying numerator and inverse power in the denominator to `\lambda_k=k^{\beta/2}/(1+\sigma^2k^{\alpha+\beta})`. Recomputed its maximum as `k^*=[\beta/((2\alpha+\beta)\sigma^2)]^{1/(\alpha+\beta)}`. Explained the discrete-frequency boundary case and removed the false claimed limit at `k=0`.
7. **Finite SVD (`prop-svd`):** corrected left singular vectors to `\mathbb R^P` and right singular vectors to `\mathbb R^N`. Corrected “right singular values” to vectors and stated the real sign/complex phase ambiguity under distinct singular values.
8. **Convolution SVD (`eq-diag-filter-finite`):** replaced undefined/inconsistent Fourier-vector notation with a unitary Fourier matrix, convolution multipliers, and correct phase-adjusted left singular vectors. Distinguished the multiplier from a unitary-normalized Fourier coefficient.
9. **Compact operators (`eq-svd-operators`):** expressed rank-one operators as `\sigma_m\langle\cdot,v_m\rangle u_m`, explicitly stated orthonormality, and corrected the operator norm from a minimum to a supremum. Distinguished finite-rank sums and Hilbert-space finite-rank approximation from the general Banach-space case.
10. **Integral-kernel examples:** corrected the kernel domain to `\Omega\times\Omega`, the convolution singular value from the unknown signal's Fourier coefficient to the kernel's multiplier, and the integration kernel from `1_{x\le y}` to `1_{y\le x}`.
11. **Pseudoinverse (`eq-pseudo-inverse`):** added the finite-dimensional/closed-range hypotheses needed for orthogonal projection onto the range and the bounded inverse formulas. Explained the restricted domain of an unbounded pseudoinverse in the infinite-rank compact case. Corrected the surjective rank from `R=N` to `R=P`.
12. **Tikhonov regularization (`sec-tikhonov`):** corrected the name, added the missing argument `f_0` of the kernel-complement projection, used an absolute spectral-filter bound, and distinguished noise amplification from impossibility of using the pseudoinverse in every noisy problem.
13. **Variational formulation:** replaced unnecessary continuity requirements on general regularizers by the relevant lower-semicontinuity/coercivity discussion. Corrected the normal-equation identity from `\mathrm{Id}_P` to the identity on the source space.
14. **Spectral proof (`prop-var-regul-svd`):** replaced inversion of a rank-deficient reduced matrix and the erroneous `V^*` output by a componentwise proof of `f_\lambda=\sum_m\sigma_m\langle y,u_m\rangle v_m/(\sigma_m^2+\lambda)`.
15. **Convolution regularization:** corrected the denominator from `|\hat h_m|^2+\sigma^2` to `|\hat h_m|^2+\lambda`.
16. **Range and kernel statements:** distinguished `\operatorname{Im}\Phi^*` from its closure in infinite dimension; removed malformed `\ker(\Phi^\perp)`. Restricted claims that the estimator lies in the kernel complement to the squared-norm/spectral penalty rather than all quadratic methods.
17. **Source condition (`eq-source-cond-init`, `eq-source-cond-quad`):** adopted one consistent convention `f_0\in\operatorname{Im}((\Phi^*\Phi)^{\beta/2})`. The coefficient bound remains `\sum_m\sigma_m^{-2\beta}|\langle f_0,v_m\rangle|^2\le\rho^2`, which now agrees with the proof and rate `\rho^{1/(1+\beta)}\delta^{\beta/(1+\beta)}`. Corrected the source-vector space and the inverse powers, added the missing kernel-complement condition, and used absolute squares. A final independent cross-check also restored the exponent `2\beta` in both squared-coefficient denominators of the bias proof.
18. **Source interpretation and saturation:** explained decay along small-singular-value directions instead of an imprecise distance from the kernel. Distinguished a supremum over all positive singular values from the bounded spectrum of a fixed operator. Qualified saturation and the impossibility of a uniform linear rate; a fixed finite-dimensional inverse can have a linear noise bound.
19. **Zero-regularization limit (`prop-gamma-conv-regul`):** supplied lower semicontinuity and the limit `\lambda\downarrow0`; corrected coercivity's quantifiers and removed an erroneous absolute value. Repaired the residual estimate to `\|\Phi f_\lambda-y\|^2\le2\lambda(J(h)-\inf J)`. Used lower semicontinuity for the limit and replaced the incorrect kernel-complement coercivity claim by joint-sublevel boundedness/no-common-nullspace conditions.
20. **General quadratic penalties:** attached uniqueness to the actual regularized objective. Corrected the observation coefficient from `\langle y,v_m\rangle` to `\langle y,u_m\rangle`, used a full right-singular basis for simultaneous diagonalization, and corrected the Laplacian from `-GG^*` to `-G^*G`.
21. **Conjugate gradients (`eq-quad-func`, `eq-iter-cgs`):** supplied the missing `1/2` in the quadratic, made the gradient notation consistent, corrected search directions from Euclidean-orthogonal to `A`-conjugate, stated exact-arithmetic termination in at most `N` steps and stopping at zero gradient, and replaced the claim of unchanged applicability to arbitrary convex functions by the need for nonlinear line searches/safeguards.
22. **TV and coarea:** distinguished a set's volume from its perimeter; defined the total variation of the distributional gradient over the correct ambient domain. Replaced the false BV-level-set formula by `\int\operatorname{Per}(\{f>t\})\,dt`, retaining the smooth level-set interpretation with its proper qualification.
23. **Smoothed TV:** restored the square inside `\sqrt{\epsilon^2+\|u\|^2}`; corrected its asymptotic remainder; required `\ker\Phi\cap\ker\nabla=\{0\}` for strict convexity and uniqueness after composition.
24. **Differentiability:** corrected the first-order expansion from `\langle f,r\rangle` to `\langle\nabla E(f),r\rangle`. Replaced the false claim that `C^1` implies a quadratic remainder by the correct locally Lipschitz-gradient condition.
25. **Energy gradients:** corrected both occurrences of `\Phi^*(\Phi y-y)` to `\Phi^*(\Phi f-y)`, distinguished gradients from differentials, and supplied the missing gradient on the TV functional.
26. **TV smoothness and strong convexity:** corrected `\|\nabla\|\le2\sqrt d`, included the missing `\lambda` in `L=\|\Phi\|^2+\lambda\|\nabla\|^2/\epsilon`, computed the pointwise Hessian, and removed the false global `\epsilon`-strong-convexity claim. Explained the actual lower eigenvalue and the fidelity-based strong-convexity condition.
27. **Inpainting:** corrected sufficient step bounds to positive bounds using `\|\Delta\|\le8`, rather than asserting equality for every grid size. Removed unrecoverable “SNR=?dB” typesetting placeholders rather than inventing measurements.
28. **Tomography:** corrected evenly spaced angles from the singular expression `\pi/k` to `\pi k/K`. Stated the Fourier convention, fixed filtered-backprojection parentheses, used complex Fourier measurements, and distinguished the idealized partial-Fourier model and its unitary normalization from the Radon operator itself. Corrected the pseudoinverse notation, invalid TV reference, and the spatial-inpainting comparison figure.
29. **Cross-references:** renamed this chapter's local gradient section/equation to `sec-ip-grad-descent`/`eq-ip-grad-desc` and the bias figure to `fig-ip-bound-regul`, avoiding collisions while preserving references to the general optimization chapter.

### Sparse Regularization (`chapters/sparse-regularization.tex`)

1. Replaced the unfinished “Ref” opening with a complete references sentence and polished terminology, grammar, and captions throughout.
2. Made orthonormality explicit wherever Parseval's identity, coefficientwise optimization, or thresholding requires it. Corrected the ideal sparsity of an approximation to `M=J_0(f_M)`, and explained that the threshold determines `M`.
3. Distinguished pointwise convergence `J_q\to\|\cdot\|_0` from convergence of bounded unit balls to the unbounded zero-pseudonorm “ball.” Replaced the categorical impossibility of minimizing a nonconvex penalty by the correct combinatorial difficulty.
4. Corrected the denoising regularizer's argument from `J_q(f)` to `J_q(g)` and used set membership for potentially nonunique minimizers.
5. Rewrote the hard/soft-threshold proof for both signs. Documented both hard-threshold minimizers at equality `|u|=\sqrt{2\lambda}`, while preserving the displayed convention that selects the nonzero value. Soft thresholding is uniquely minimizing.
6. Qualified the universal Gaussian threshold and empirical threshold choices; removed the unconditional asymptotic-optimality claim independent of a signal class.
7. Corrected nonorthogonal analysis/synthesis terminology, operator order, coefficient dimension `Q`, and the empty generic-formulation reference. Explained dictionary spanning and the extended-value synthesis prior, rather than claiming analysis recovery theory is mostly an open problem.
8. Corrected the noisy feasible-set caption from `\|\Phi f=y\|` to `\|\Phi f-y\|`.
9. **Linear program:** changed the nonexistent `a=x_+-x_-` variable to `x`; distinguished the optimizing pair from its difference `x^*=x_+^*-x_-^*`; completed the Douglas–Rachford reference and qualified computational cost claims.
10. **Projected gradient:** used `A=\Phi\Psi`, not `\Phi`, in coefficient-space objectives and gradients. Corrected the factor-two error in the allowable step: the split-variable smooth gradient has Lipschitz constant `2\|A\|^2`, giving `\tau<1/\|A\|^2`. Clarified componentwise scalar additions.
11. **ISTA surrogate (`eq-ista-surrog`):** replaced `E(x,x)=0` by `E_\tau(x,x)=E(x)`; restored the missing square and the identity dimension `Q`; stated nonnegativity of the correction and tangency at `x=x'`.
12. **ISTA proof (`eq-ista`):** corrected `AA^*x'` to `A^*Ax'`, restored `1/2` in the completed square, removed a mixed fixed/varying step, and corrected the proximal threshold to `\rho\lambda`. Distinguished the surrogate condition `\tau\le1/\|A\|^2` from the broader forward–backward convergence condition `\tau<2/\|A\|^2`.
13. **Sparse deconvolution:** replaced the empty historical citation and inaccurate density-change wording with the sparse acoustic-impedance/reflectivity model. Used the adjoint kernel `h^\vee[n]=\overline{h[-n]}` in the gradient; explained when symmetry permits using `h` itself.
14. Distinguished convergence of objective values from convergence of iterates and cited the role of the forward–backward theorem.
15. Corrected the effective Gaussian-blur bandwidth count to dimension-dependent `\mu^{-d}`, removed the reference to a commented-out implementation table, and replaced invalid Sobolev/ridge equation references.
16. Explained the difference between synthesis ISTA in a redundant frame and the generally invalid assumption that analysis–threshold–synthesis is the proximal map of an analysis penalty.
17. Corrected inpainting to a small-positive-λ approximation of the exactly constrained problem; made the orthonormal-basis assumption explicit in the pixel-space threshold iteration.
18. Removed unknown SNR placeholders from figures, corrected remaining iteration references, and supplied a plain-text PDF bookmark alternative for the `\ell^1` heading.

### Theory of Sparse Regularization (`chapters/sparse-theory.tex`)

1. **Objective convention (`eq-lasso-lagr`):** removed the erroneous second factor of `\lambda`, obtaining `\|y-Ax\|^2/(2\lambda)+\|x\|_1`, consistent with the rest of the book, its dual, and its threshold parameter.
2. Corrected the noisy recovery problem from `\mathcal P_0` to `\mathcal P_\lambda`, and kept the constrained `\lambda=0` problem separate from a division-by-zero formula.
3. **Polytope characterization (`prop-polytope-proj`):** required `x_0\ne0`; completely repaired the contrapositive proof, its right inverse, incorrect norms, radii, and inequality. Treated `Ax_0=0` separately. Clarified spherical subdivision/general-position assumptions and signed columns in the empty-cap property.
4. **Optimality (`prop-first-order-lagr`):** corrected the off-support condition to `\|\eta_{I^c}\|_\infty\le1`, corrected the subdifferential scaling, and described the normal space of an affine constraint correctly.
5. **Sparse minimizer existence:** repaired the support-reduction proof by moving along a kernel direction to an endpoint of a constant-sign interval; corrected the reversed stopping condition and distinguished penalized/constrained feasibility.
6. **Uniqueness:** corrected the predictor from `\Phi x` to `Ax`, the convexity proof's `\|x_2\|_2` to `\|x_2\|_1`, and duplicate local-shape labels. Stated the general-position hypothesis for almost-sure noisy Lasso uniqueness.
7. **Duality (`eq-dual-lasso`):** rewrote the derivation with a separate residual variable and correctly dimensioned Lagrangian; fixed `Ax` versus `A^*x`, the ℓ1 conjugate, missing dual feasibility, and primal–dual residual relations. Distinguished strict concavity of the dual maximum from convexity and stated feasible zero-parameter dual attainment. Corrected the Fenchel–Rockafellar reference and spelling.
8. **Bregman divergence:** supplied the missing equality `D_\eta(x_0\mid x_0)=0`, replaced an editorial drawing TODO with the existing figure, corrected the squared-distance description, and specified positive reference entries and the zero convention for relative entropy.
9. **Bregman rate proof (`thm-bregman-rates`):** restored the missing `1/2` in the completed square. The corrected argument yields the improved residual bound `\|Ax_\lambda-y\|\le\|w\|+2\lambda\|p\|`, replacing the earlier looser `1+\sqrt2` constant. Propagated the revised constant into the subsequent norm-rate proof.
10. Corrected the quadratic source-condition comparison to β=1 under the book's repaired `(A^*A)^{\beta/2}` convention, and resolved the explicit inconsistency TODO. Qualified the converse multiplier/source-condition claim by a subdifferential qualification.
11. **Linear ℓ1 rates:** corrected missing coordinate subscripts; stated nonnegative rather than strictly positive Bregman summands; adopted the empty-complement norm convention; corrected the input/output order of mixed operator norms; introduced `p`, repaired a wrong `x` argument, and made the constants' dependence and positive-noise/positive-λ assumptions explicit.
12. Replaced the claim that a fixed-reference error estimate proves a pairwise Lipschitz solution map by its actual linear error conclusion.
13. **Minimum-norm certificate:** corrected the noiseless problem in its assumption and supplied a direct boundedness/cluster-point proof. Explained that the minimized norm is that of the data-space vector `p`.
14. Distinguished the true minimum-norm certificate, Fuchs precertificate, and vanishing-derivative precertificate by separate labels (`eq-minnorm-certif`, `eq-fuchs-precertif`, `eq-vanishing-precertif`); corrected their constraints, pseudoinverse notation, and prerequisites.
15. **Support stability (`thm-support-stable`):** replaced an assumption that incorrectly permitted arbitrarily large `\lambda` by small `\lambda` and small `\|A^*w\|_\infty/\lambda`. Stated exact support/sign equality, the correct error `O(\|A^*w\|_\infty+\lambda)`, and the additional proportional choice needed for a linear noise rate.
16. Corrected the sign of the orthogonal noise residual, malformed delimiters, and a displayed inequality that bounded a nonnegative quantity by a negative quantity. Derived the actual condition `R\delta<S\lambda`.
17. Corrected the reciprocal error in choosing `\lambda`: use a strict multiple of `R\delta/S`, not `S\delta/R`. Propagated the strict slack into the explicit error bound and retained the valid triangular admissible region.
18. **Arbitrary-noise sparsistency:** completed the unfinished ERC subsection. Defined `\mathrm{ERC}(I)=1-\max_{j\notin I}\|A_I^+a_j\|_1`, stated an explicit noise-dependent lower bound on `\lambda`, and proved uniqueness/support containment by restricted minimization and strict off-support dual feasibility. Distinguished containment from equality of supports. This directly proved formulation avoids relying on more delicate assertions affected by the published Tropp corrigendum.
19. **Grid-based super-resolution:** corrected the convolution-kernel orientation, finite sampling indices, number `2f_c+1` of Fourier measurements, covariance matrix indices/dimensions, and the number of kernel functions.
20. Replaced the false assertion that one finite grid detects all continuous overshoots by a statement about increasingly dense grids. Explained that a zero derivative at a support point is necessary but not sufficient for a bounded certificate.
21. Repaired the vanishing-derivative interpolation problem, undefined kernel symbol, missing minimization equality, support-location variables, and grid/continuous adjoint distinction. Added the needed independent-constraint assumption.
22. Corrected the convolution certificate figure reference and the injectivity condition needed to apply the linear-rate theorem. Qualified grid-free and grid-support stability claims instead of presenting them as automatic.
23. Polished the prose and captions, corrected “Sparsitency” to “Sparsistency,” removed active TODOs, and added the plain-text PDF bookmark alternative for the ℓ1 heading.

### Compressed Sensing (`chapters/compressed-sensing.tex`)

1. Polished the motivation and replaced dated categorical assertions about physical realizability with an explanation of the distinction between ideal random matrices and hardware constraints.
2. **Single-pixel camera:** corrected the scrambled mirror/pixel definition, measurement summation indices, and compression ratio from `P/Q=6` to `Q/P=6`. Described off mirrors as directing light away from the sensor. Explained how centered signed measurements differ from raw 0/1 mirror patterns.
3. Corrected the dictionary composition from `\Psi\Phi` to `\Phi\Psi`, specified the dictionary dimension, and stated the correct right-orthogonal Gaussian invariance.
4. Replaced the false universal bijection between penalty and residual parameters by the qualified Lagrange-multiplier correspondence, including flat segments and limiting endpoints.
5. **Random polytopes:** corrected an `s`-sparse vector's face dimension to `s-1`, the dimensionally impossible sparsity threshold (P/N) to a proportional count, and the negative logarithmic asymptotic to inverse-logarithmic decay. Removed unsupported precise strong-threshold numbers while retaining the appropriately approximate weak-threshold example.
6. **Random matrices:** specified the (1/P) Gaussian entry variance; corrected the inverse Gram matrix and law-of-large-numbers terminology.
7. **Marchenko–Pastur law:** replaced an undefined probability of “eig” by convergence of the empirical eigenvalue measure; stated `s/P\to\beta\in(0,1)`, the correct density/support, and why an atom at zero requires separate treatment for `s>P`. Split the duplicated figure labels.
8. **Concentration (`eq-talagrand`):** replaced the incorrect Gram-matrix concentration threshold and exponential scale by the correct Gaussian singular-value bounds with probability `1-2e^{-Pt^2/2}`. Derived the corresponding Gram bound `2a+a^2`.
9. **Coherence:** separated entrywise coherence from the induced infinity norm; included the endpoint μ=1; removed division-by-zero in the `s=1` assumption; restored absolute values in correlation sums and corrected the Neumann-series identity and missing comparison with 1.
10. Qualified support containment under coherence by the explicit ERC/noise-dependent penalty condition. Stated the Welch bound in its `N>P` regime and replaced an overly precise Gaussian coherence equivalence by a high-probability order bound for normalized columns.
11. **Sub-Gaussian model:** specified independent centered normalized entries and corrected the equation numbering so the tail definition is referenceable.
12. **Nonuniform certificate theorem (`thm-cs-etaf`):** replaced an unusable formula with an unrestricted/undefined δ and constants by the standard valid `P\ge Cs\log(2N/\epsilon)` sufficient condition. Explicitly required the vector to be fixed independently of the matrix, obtained injectivity and strict off-support bounds, and linked to the corrected small-λ support theorem.
13. Rewrote its proof using conditioning on `A_I`, an inverse-Gram bound, and a union bound. Separated general sub-Gaussian tails from the special Gaussian distribution; corrected the count to `N-s` off-support columns and the distinction between `p_F` and its correlations.
14. Corrected the inverse-Wishart/chi-square argument and its direction: `Ps/\|p_F\|^2\sim\chi^2_{P-s+1}`, so the upper bound uses a factor `1+\delta`, not `1-\delta`. Clearly identified `2s\log N` as a Gaussian leading-order asymptotic, not a universal finite-sample constant.
15. **RIP/RO lemma:** normalized vectors before polarization to obtain the product of norms, rather than incorrectly deducing it from their squared-norm sum; supplied the previously omitted proof of the reverse inequality.
16. **Random RIP theorem:** supplied parameter ranges and the correct sample count `P` in its failure probability. Qualified the required tail assumptions; finite covariance alone does not give sub-Gaussian concentration. Corrected the smallest eigenvalue versus its deficit from 1.
17. **RIP certificate theorem (`thm-rip-dualcertif`):** repaired the iterative proof using `\delta_{2s}` on interpolation supports of size `2s`. The displayed sufficient hypothesis is correspondingly strengthened to `\delta_{2s}+\theta_{s,s}+\theta_{s,2s}<1`. This is deliberately a valid slightly more conservative theorem, not a claim to have proved the original sharper statement with the flawed argument.
18. Corrected the undefined `\delta_{s,2s}`/`\delta_{s,s}` quantities, support cardinality, initial term in the bound on `p`, convergence ratio, pointwise cancellation argument, and index sets. Removed the drawing TODO and retained the implication `\delta_{3s}<1/3`.
19. **Interpolation lemma:** corrected the inverse from `A_JA_J^*` to `A_J^*A_J`, the coordinate `c_J`, and the tail denominator from `\sqrt s` to `\sqrt{s'}`. Allowed fewer off-support coordinates than `s'`.
20. **Uniform stability (`thm-rip-final`):** replaced a noisy uniqueness claim unsupported for general sub-Gaussian ensembles by the standard constrained recovery guarantee for every minimizer, with error `C_0\epsilon+C_1\|x_0-x_{0,s}\|_1/\sqrt s`. Used `\log(eN/s)` and explicit simultaneous quantifiers over signals and noise. Exact sparse noiseless recovery is unique; noisy minimizers need not be.
21. Retained the connection to penalized recovery with its proper parameter choice and qualifications, and replaced the crude operator-norm tail amplification by the standard compressibility term.
22. **Subsampled Fourier/orthogonal measurements:** supplied the missing `\sqrt{N/P}` normalization, without which the stated RIP claim is false. Restored `N^{-1/2}` in the Fourier basis vectors, allowed complex measurements, clarified the bounded-orthonormal-system parameter, corrected “Dirac” and “Vershynin,” and stated the dependence of sample-size constants on the RIP tolerance and failure exponent.
23. Completed invalid threshold/sampling/RIP references, updated the Fuchs precertificate reference, corrected figure captions, and removed unsupported precision in claims about practical constants.

### Final prose refinement

A final read of the active prose in the variational-priors, inverse-problems, sparse-regularization, and sparse-theory chapters refined articles, subject–verb agreement, sentence transitions, proof introductions, and captions. It shortened repetitive descriptions of gradient flows, SVD, conjugate gradients, synthesis priors, and certificate constructions; clarified oracle parameter selection and the relation between ISTA and forward–backward splitting; and replaced awkward phrases such as “but it is rarely the case,” “This requires to decide,” and “this offers an alternative point of view.” The mathematical expressions and source labels were preserved. Compressed-sensing prose received its final refinement in the main editorial pass.

### Mathematical verification and consulted primary sources

Independent numerical sanity checks, with a fixed random seed, passed for the corrected periodic Laplacian multiplier and discrete adjointness; the smoothed-TV derivative against a centered directional difference; the ISTA completed-square identity and soft-threshold optimality conditions; the low-noise support formula and dual residual; and the oracle-filter maximizer. These are consistency checks, not substitutes for the analytical proofs.

The corrections to the more specialized compressed-sensing claims were checked against primary mathematical sources, with proofs written in the book's notation:

- Joel A. Tropp, [Just Relax: Convex Programming Methods for Identifying Sparse Signals in Noise](https://users.cms.caltech.edu/~jtropp/papers/Tro06-Just-Relax.pdf), 2006. The basic ERC support-containment statement used here is proved directly; the [publication record](https://authors.library.caltech.edu/records/v3yk2-6z279) also identifies a later corrigendum.
- Emmanuel Candès and Terence Tao, [Decoding by Linear Programming](https://candes.su.domains/software/l1magic/downloads/papers/DecodingLP.pdf), 2005, for the iterative certificate method. The book's repaired proof states the slightly stronger hypothesis actually established by its interpolation lemma.
- Emmanuel Candès, [The Restricted Isometry Property and Its Implications for Compressed Sensing](https://candes.su.domains/publications/downloads/RIP.pdf), 2008, for the robust constrained RIP estimate. The main editor corrected the existing bibliography entry's year and journal fields.
- Roman Vershynin, [Non-asymptotic random-matrix discussion](https://www.math.uci.edu/~rvershyn/papers/conjectures-2006.pdf), for Gaussian singular-value concentration and its scaling.

### Scope limits and final validation

No active mathematical TODO, empty citation, or knowingly false statement remains in these five reviewed sources. This is a detailed editorial and mathematical pass, not a formal verification of every background theorem or an independent reproduction of all numerical figures. General random-matrix and sampling theorems intentionally retain proof outlines/references where the original exposition did so. Historical image assets, their plotted empirical curves, and numerical SNR values have not been regenerated; unknown SNR placeholders were removed instead of guessed. The final layout, reference, bibliography, and PDF checks are recorded above and in the build report.

## Chapters 13–19: machine learning, optimization, and convex analysis

This review covers all active text, definitions, statements, proofs, examples, equations, and captions in `machine-learning.tex`, `perceptrons.tex`, `deep-learning.tex`, `convex-analysis.tex`, and `optim-nonsmooth.tex`, together with the active inputs `sec-pca-theory.tex`, `machine-learning-sec-pac.tex`, `sec-optim-smooth.tex`, `sec-regul.tex`, `sec-stochastic-optim.tex`, and `sec-autodiff.tex`. The two optimization chapter wrappers were inspected to establish this coverage. `sec-optim-mirror.tex` is not actively included and was not reviewed. Commented-out and `\if 0` material is outside the active-content review.

Across these files, prose was corrected for grammar, spelling, agreement, punctuation, mathematical terminology, and misleading generalizations. Repeated explanations were tightened; statements were separated from the assumptions needed to justify them. Mathematical changes are itemized below. This records substantive changes and groups repeated copyedits rather than listing each punctuation edit separately.

### Machine learning

#### PCA and clustering

- Made the data convention consistent: observations are rows of the `n × p` matrix, centering uses the `n`-vector of ones, and covariance is `XᵀX/n`. Stated the finite-second-moment assumption behind covariance estimation.
- Used the empirical mean in compression and reconstruction, indexed the retained components from `1` to `d`, and distinguished covariance eigenvalues from singular values. Explained extension to an orthonormal basis when the requested dimension exceeds the data rank.
- Corrected the PCA optimality proof's matrix dimensions and transposes. The unnormalized matrix is `C = XᵀX = n Ĉ`, with eigenvectors in `V`. Restored the missing minus sign in the projection residual and the constant Frobenius-energy term in the trace identity.
- Replaced the picture-based conclusion of the eigenvalue bound with a proof valid in arbitrary dimension: the constraints `0 ≤ βᵢ ≤ 1`, `Σβᵢ ≤ k` imply `Σλᵢβᵢ ≤ Σᵢ≤k λᵢ`. Kept the diagram as an illustration.
- Specified how ties and empty clusters affect Lloyd iterations and qualified finite-termination claims. Corrected the direction of cluster reallocation toward regions with high distortion.
- Corrected k-means++ sampling from inverse-squared distances to `dᵢ²/Σⱼdⱼ²`; the distance now includes every already selected center. Added the zero-distance stopping case and corrected the variables in the expected-cost theorem.
- Corrected Voronoi-cell notation and distinguished finite configurations from the asymptotic planar hexagonal-lattice picture.

#### Risk, regression, and model selection

- Corrected the conditional-risk integrand, expectation parentheses, and target-function notation. Added integrability and positive-density conditions to conditional-mean formulas.
- Distinguished validation data from final test data and clarified how cross-validation estimates held-out performance.
- Corrected a class-label example in which both healthy and pathological cases had been assigned the same label.
- Distinguished an uncentered second-moment matrix from a covariance matrix, corrected `E[XᵀX]` to `E[XXᵀ]` for column-valued observations, and corrected the asymptotic direction from `n → 0` to `n → ∞`.
- Corrected the regression/inverse-problem correspondence: the design operator is `X/√n`, not the covariance matrix. Qualified the claim that sampling noise behaves like deterministic inverse-problem noise.
- Corrected the ridge system to `(Ĉ + λI)⁻¹ û` under the chapter's averaged-loss convention; removed an inconsistent extra factor of `n`. Corrected the comparison between the primal `p × p` and dual `n × n` systems.
- Replaced an unsupported ridge estimation-rate theorem by a proved bound with explicit bounded-feature, moment, and source assumptions. For `β̄ = C^γz`, `0 < γ ≤ 1`, `‖z‖ ≤ ρ`, the bound is `E‖β̂−β̄‖² ≤ 2ρ²λ^(2γ) + 4B²/(nλ²)`, where `B² = κ²(EY² + κ²‖β̄‖²)` and `‖X‖ ≤ κ`. Balancing the two terms gives the stated `n^(−γ/(γ+1))` rate. Explained the qualification limit of ordinary Tikhonov regularization for larger source exponents.
- Corrected nearest-neighbor consistency to require the number of neighbors to increase while its ratio to sample size tends to zero. Corrected the class-histogram numerator and validation normalization.

#### Logistic classification, SVMs, and kernels

- Made binary probabilities and labels consistent. Corrected the Bernoulli likelihood and its negative logarithm to obtain `log(1 + exp(−ys))`.
- Corrected the logistic gradient to contain `−y`, not an unrelated score variable, and retained the `1/n` factor from the empirical mean.
- Replaced the assertion that strict convexity guarantees a minimizer by the correct conclusion of at most one minimizer. Explained nonattainment for separable data and existence/uniqueness after positive quadratic regularization.
- Corrected the hinge slack constraint to `yᵢsᵢ ≥ 1−uᵢ` and qualified broad comparisons between logistic regression, SVMs, and sparsity penalties.
- Switched to the corrected vector loss plot generated during the main review. Its caption now defines the negative margin, explicitly counts zero margin as an error, and explains that unscaled logistic loss is not an upper bound for 0–1 error; division by `log 2` supplies that normalization.
- Made multiclass softmax logits, normalization, averaged objective, and gradient consistent. Distinguished scalar log-sum-exp from rowwise application and explained the stable shift used in numerical evaluation.
- Corrected Gram-operator orientation, polynomial-feature dimension and products, and the normalization of the Gaussian feature lift.
- Distinguished positive semidefinite from positive definite kernels. Stated Bochner's condition using a nonnegative finite spectral measure, rather than strictly positive density. Added the nonzero-diagonal requirement for kernel normalization.
- Replaced an incorrect representer-theorem derivation by an orthogonal-decomposition proof. Stated finite nonnegative convex loss and positive regularization hypotheses, uniqueness of the Hilbert-space minimizer, and possible nonuniqueness of coefficients for a singular Gram matrix.
- Corrected computational claims: a dense direct solve is cubic; quadratic cost is a dense matrix-vector product or iteration, not a universal total solve complexity.

### Statistical learning / PAC theory

- Stated i.i.d. sampling and measurability/existence assumptions. Corrected the conditional-risk variable and the 0–1 loss from equality to misclassification.
- Distinguished real-valued scores from class labels and stated when convexity of the model class makes empirical risk minimization convex.
- Handled the possibility that logistic conditional optima are attained only at infinite scores when a conditional class probability is zero or one, and allowed either Bayes label at a tie.
- Corrected the natural-log logistic calibration bound from `√s` to `√(2s)` and stated assumptions on the convex margin loss. Expressed the final calibration conclusion in terms of actual 0–1 risk with fixed tie breaking.
- Removed claims that confuse distribution-free rates with universal consistency, or imply that a random estimation error must decrease monotonically with sample size.
- Replaced an invalid bound on the absolute supremum deviation by two one-sided deviations. Applied symmetrization to each side and a union bound to McDiarmid's inequalities. The resulting excess-risk bound has the correct constants and logarithmic confidence factor.
- Specified independent symmetric Rademacher signs, restored the expectation in the linear-class dual-norm bound, and retained the dimension-independent Hilbert-space second-moment estimate.
- Rewrote the SVM margin example as a bound valid uniformly for every score in a norm ball. The comparison uses the bounded ramp loss and empirical hinge loss; it no longer assumes that a hinge minimizer also minimizes the ramp loss. The feature moment, radius, margin, and confidence terms are explicit.
- Corrected the assertion of strict convexity in kernel coefficients when the Gram matrix can be singular.

### Smooth optimization

- Corrected coercivity to concern `‖x‖ → ∞`; separated existence from strict convexity; stated distinct-point requirements for strict convexity.
- Corrected the signs in the supporting-hyperplane argument and distinguished the differential, directional derivative, and gradient. Clarified Fermat's rule and the role of differentiability; the one-dimensional cubic example has an inflection point rather than a local extremum.
- Restored the square in least-squares objectives and corrected rank/injectivity conditions and covariance scaling.
- Corrected the steepest-descent direction and the nonzero-gradient qualification. Restored missing perturbation variables in Taylor and line-search formulas and corrected the Armijo small-step argument.
- Standardized quadratic objectives with the factor `1/2`, corrected minimizer notation and spectral formulas, and distinguished a condition number from its reciprocal.
- For rank-deficient quadratic problems, corrected the claim that all iterates are orthogonal to the kernel: the initial kernel component is preserved. Added the existence condition and used the positive part of the spectrum.
- Corrected Hessian hypotheses and second-order necessary conditions; repaired the Newton Taylor-expansion endpoint.
- Corrected both smoothness and strong-convexity quadratic bounds to use the gradient at the expansion point and displacement `x−x′`.
- Added convexity, smoothness, and minimizer-existence assumptions to convergence statements; corrected iterate indices in telescoping and contraction inequalities and a constant in the `1/(k+1)` bound.
- Replaced a conflation of heavy-ball and Nesterov acceleration. The accelerated scheme, momentum parameter, squared-distance constant, strong-convexity qualification, and restart/tuning discussion are now consistent.
- Corrected the accelerated continuous-time scaling to time step `√s`, repaired the damping and gradient coefficients, and removed a false uniqueness claim about the damping value `3`.

### Regularization and stochastic optimization

- Added lower semicontinuity, coercivity, nonnegativity, and feasibility assumptions to the vanishing-regularization argument. Corrected its residual estimate to a bound on the squared residual.
- Made the ridge penalty's `1/2` normalization agree with its linear system and distinguished the minimum-norm least-squares limit from the exactly feasible case.
- Corrected the Lasso input and norm notation and the ISTA surrogate's squared norm, touching value, matrix size, and transpose.
- Corrected empirical-distribution weights to `1/n` and distinguished pointwise laws of large numbers from convergence of minimizers.
- Corrected the stochastic chain rule and logistic derivative, and stated conditional unbiasedness with respect to the iterate history.
- Replaced the SGD convergence proof with a valid induction. Strong convexity is combined with a conditional second-moment bound along the iterates, avoiding an impossible global bounded-gradient assumption on a strongly convex function over all of Euclidean space. With step `1/[μ(k+2)]`, the explicit mean-square bound has denominator `k+2`; smoothness then gives an objective bound.
- Corrected stochastic averaging indices and weights. Qualified statements that averaging always accelerates every problem.
- Corrected SAG's stored-gradient initialization, incremental average factor `1/n`, and evaluation iterate. Distinguished individual-gradient smoothness from objective strong convexity and qualified the convex-case averaging guarantee.

### Automatic differentiation

- Clarified that automatic and symbolic differentiation use the same chain rule but differ in representation and reuse of intermediate computations. Stated operation-cost and domain assumptions and corrected finite-difference evaluation counts.
- Corrected computational-graph indices, parent/child direction, and reverse traversal. Added the missing multiplier in the feedforward weight adjoint.
- Corrected recurrent-network dimensions, terminal adjoint data, parameter-gradient integrand, and time-step scaling for controls.
- Corrected staggered inverse updates and qualified claims of symplecticity and stable numerical reversibility.
- Corrected fixed-point differentiation to `(I−Dₓh)⁻¹Dθh`, with differentiability and invertibility assumptions. Distinguished an equilibrium derivative from the derivative of finitely many solver iterations.
- Corrected the argmin vector field, Hessian, and mixed-derivative signs.
- Rewrote the Sinkhorn differentiation example using positive marginals, positive Gibbs kernel, genuinely alternating updates, and consistent dimensions. Corrected the dual objective to a maximization with the entropy factor, the derivatives through both scaling steps, and the marginal gradients to `ε log u` and `ε log v`. Explained additive gauge freedom and the need to fix it before implicit inversion.

### Convex analysis

- Corrected convex-set quantifiers, strict-convexity conditions, and indicator-domain notation.
- Repaired the minimizer-existence proof so it constructs a nonempty compact sublevel set without assuming a minimum exists first. Separated properness, lower semicontinuity, coercivity, existence, and uniqueness.
- Stated the subdifferential outside the domain and the interior condition for differentiability/singleton equivalence. Distinguished dual pairing from its Hilbert-space identification.
- Added the continuity/domain qualification for the subdifferential sum rule and the relative-interior qualification for a linear composition. Restricted positive scaling to avoid the extended-valued `0 × ∞` ambiguity. Completed the Lasso optimality conditions.
- Corrected the conjugate of `f(x)=xᵀAx/2−bᵀx` to `(u+b)ᵀA⁻¹(u+b)/2`, including symmetry/positive-definiteness assumptions. Corrected translation, scaling, and star notation.
- Clarified the lsc convex-envelope and affine-minorant qualifications for biconjugation. Corrected smoothness/strong-convexity conjugacy and its reciprocal constant.
- Corrected infimal convolution from a supremum to an infimum, with the necessary conjugacy/closure qualifications. Restored the `1/μ` factor in the Moreau-envelope gradient.
- Completed the affine-duality subsection with consistent signs.
- Corrected the primal domain, use of suprema for multipliers, and the weak-duality feasibility argument. Corrected KKT's nonnegative inequality multipliers and complementary slackness.
- Corrected affine projection to `z−A⁺(Az−y)` and used the pseudoinverse to cover non-full-row-rank matrices.
- Stated properness, lower semicontinuity, and relative-interior assumptions in Fenchel–Rockafellar duality. Used infima in the inner dual derivation: domain qualification does not ensure attainment of every fixed-multiplier inner problem.

### Nonsmooth optimization

- Added the proper/lsc/convex and existence assumptions used by convergence claims. Corrected Taylor bounds that require a Lipschitz gradient rather than mere differentiability.
- Corrected subgradient step-size conditions, the best-iterate bound, and horizon-dependent `1/√N` tuning; distinguished step-size convergence from iterate convergence.
- Stated the assumptions and step range for projected gradient and qualified its objective and strong-convexity rates.
- Corrected the log barrier from concave to convex, its strict domain, dimensions, and data variables. Qualified self-concordance and complexity claims.
- Restored the missing proximal parameter, threshold, and input. Corrected nonexpansiveness versus strict contraction and made hard thresholding set-valued at a tie.
- Corrected affine perturbation and translation identities for proximal maps and repaired their signs and proof points.
- Corrected the proximal-point/resolvent sign and assumptions, and the forward-backward surrogate's linear term and touching value. Stated the conditions supporting the convergence rates.
- Corrected the dual forward-backward gradient, TV projection notation, discrete Laplacian sign and norm bound, and the hinge conjugate proximal map to clipping `u+τ1` into `[0,1]`.
- Added the subdifferential solvability assumption for Douglas–Rachford. Its shadow sequence now converges to a minimizer that may depend on initialization, rather than to an arbitrarily prescribed witness minimizer.
- Corrected the graph-projection signs and completed the proof.
- Corrected ADMM's penalty factor, subproblem signs, and twisted-prox argument. Explained when a pseudoinverse recovers the primal variable and why a nullspace component cannot generally be discarded.
- Corrected the constrained-TV Lagrange multiplier optimization to a supremum. Specified block soft thresholding for isotropic TV and componentwise soft thresholding for anisotropic TV.
- Gave the TV primal-dual splitting that actually avoids a linear-system inverse. Corrected Chambolle–Pock's parentheses, extrapolation index, operator norm, saddle-point assumption, and the standard convergence choice `θ=1`.

### Perceptrons and approximation theory

- Corrected the layer recurrence, dimensions, and polynomial-degree dependence on depth. Clarified that linear networks can still be trained with supervision.
- Corrected the squared loss, current-iterate gradient, perturbation parameter, remainder, and both shallow-network matrix gradients.
- Repaired the uniform-approximation proof: shifted cosines separate points; sigmoid staircases approximate each cosine through small jumps on compact intervals. This avoids the false assertion that one sigmoid approximates an entire cosine and avoids unjustified uniform convergence of Fourier partial sums.
- Corrected the Fourier transform variable, sign, normalization, and inverse convention. Selected the continuous Fourier-inversion representative so the theorem remains valid for sampling measures with atoms.
- Distinguished the Fourier-defined Barron class from activation-dependent spaces, and ordinary Fourier densities from the measure-valued transforms of ridge functions.
- Corrected the Gaussian dimension dependence and the Sobolev sufficient condition from `s>d/2` to `s>d/2+1`, including all weak derivatives up to order `s`. Added decay requirements to ridge examples.
- Stated a valid Barron theorem for a bounded nondecreasing sigmoid and arbitrary probability measure on a bounded ball, with an output bias, `2R‖f‖B/√n` infimum bound, coefficient bound, and finite-parameter approximation up to arbitrary extra `ε`.
- Replaced an unjustified exact probability-measure representation on finite sigmoid parameters by a closed-convex-hull representation of centered neurons. The Fourier-to-dictionary argument is explicitly labeled a proof sketch.
- Corrected the Monte Carlo argument using the Hilbert-space variance identity, distinguishing squared error `1/n` from norm error `1/√n`, and extended it to closure by approximation.
- Qualified claims about mean-field global convergence and the difficulty/attainment of the single-neuron optimization oracle.
- Corrected the target and measure in the convex optimization formulation and the first-variation perturbation notation.
- Repaired the Frank–Wolfe gap sign, step-size indexing, diameter-versus-radius constant, base case, and induction. Corrected the kernel supremum domain and the Lipschitz bound `L≤M²` for the centered-neuron dictionary.

### Deep learning

- Corrected layer dimensions and index conventions; separated nonconvex SGD stationarity claims from the convex convergence theorem. Clarified shared-parameter gradient accumulation and appropriate regression output layers.
- Corrected the translation-equivariance theorem and convolution representation: the impulse response itself is the forward filter; reversal belongs in the adjoint. Stated periodic boundary assumptions and the loss of full translation equivariance under striding.
- Replaced incorrect CNN backpropagation equations by separate activation adjoints and filter gradients, with the proper reversed convolution factors, derivative multiplier, initialization, and backward indices. Explained support restrictions and adjoint upsampling.
- Clarified residual dimensions and the ODE interpretation.
- Replaced the incomplete batch-normalization description by the mean, variance, stabilizer, and learned scale/shift formula; distinguished training from inference and convolutional sharing. Added the original Ioffe–Szegedy citation.
- Replaced the unnormalized transformer expression by scaled dot-product attention with query/key/value dimensions, row-normalized weights, and the `1/√dₖ` factor. Explained dense quadratic sequence cost, positional information, and permutation equivariance; added the original Vaswani et al. citation.
- Qualified scattering stability by its wavelet/deformation assumptions and clarified its fixed-feature role.

### Cross-references, sources, and verification

- Resolved duplicate labels for the two classification figures, kernel formulas, regularization/ISTA equations and figures, nonsmooth gradient/problem equations, recurrent backpropagation, and the deep-network comparison figure. Added missing labels for reverse mode, feedforward differentiation, multiclass softmax, and the strong-convexity theorem alias. Updated local references consistently.
- Supplied text alternatives for mathematical section titles used in PDF bookmarks and replaced the PCA proof heading by plain text.
- Checked the original batch-normalization paper, transformer paper, Barron approximation paper, and Chizat–Bach mean-field result when revising specialized claims. Relevant sources: [Ioffe and Szegedy](https://proceedings.mlr.press/v37/ioffe15.html), [Vaswani et al.](https://papers.neurips.cc/paper/2017/file/3f5ee243547dee91fbd053c1c4a845aa-Paper.pdf), [Barron](https://pages.cs.wisc.edu/~brecht/cs838docs/93.Barron.Universal.pdf), and [Chizat and Bach](https://arxiv.org/abs/1805.09545).
- Independently cross-checked the conjugacy, Fenchel/KKT, projection, ADMM, Barron, and Frank–Wolfe arguments with the other chapter reviewers. Incorporated their final qualifications on inner infima, initialization-dependent Douglas–Rachford limits, TV conventions, and atomic sampling measures.
- Ran numerical checks of the quadratic conjugate on random positive-definite systems; affine and graph projections for tall and wide matrices; binary and multiclass logistic gradients; both shallow-network weight gradients; convolutional activation/filter gradients; the implicit fixed-point Jacobian; primal/dual ridge identities; and the SGD/Frank–Wolfe induction inequalities. All passed. Gradient/Jacobian checks used central finite differences.
- Complete and standalone builds, fonts, layout, corrected graphics, and PDF diagnostics passed the checks recorded above. These checks support the corrected formulas but do not constitute a formal verification of the entire book.

## Coverage limits

Commented-out chapter includes, archived material, and inactive drafts were excluded. In particular, the mesh, optimal-transport, diffusion, and mirror-descent material is not part of the active book. Existing experimental illustrations were not all regenerated, and legacy MATLAB numerical implementations were not independently audited. Background results that were originally cited or sketched remain cited or sketched where appropriate; this revision is a detailed mathematical and editorial review, not a formal verification of every theorem.

## Hand-drawn figure reconstruction pass

Completed the requested original/TikZ review edition on 2026-09-05. The active chapter include list was used to inventory hand-drawn artwork. There are **91 reconstructions in 14 chapters**, each accompanied by its original on a dedicated landscape comparison page in both the complete book and the affected standalone chapter. Existing photographs, computed plots, and already typeset diagrams were excluded from this reconstruction pass.

### Presentation and build changes

- Added editable standalone TikZ sources and matching vector PDFs under `figures/tikz/`, with the same New PX text/mathematics and Source Sans Pro design as the book.
- Added a shared drawing preamble, four source manifests, and generated per-chapter review sections. Each review records the original path, mathematical context, reconstruction ID, and interpretation notes.
- Added automatic compilation of changed TikZ sources and strict diagnostic checks to the normal book build. Review sections are regenerated from the manifests, so review text has one authoritative source.
- Added contents/bookmark links to the comparison sections and clickable equation, figure, and proposition references in the review context. Original figures and their numbering are retained in the chapter text for review; no original asset was overwritten.
- Used the surrounding formulas to resolve signs, normalizations, matrix dimensions, endpoints, orthogonality, and limiting cases. Reconstructed point locations and schematic profiles are identified as illustrative. Fourteen interpretations are explicitly marked for checking where details could not be recovered confidently.

### Further text corrections found while reconstructing figures

- In `chapters/machine-learning.tex`, repaired the sentence that had interchanged the estimator and sample-size notation: the estimator is written `f-hat_n` to emphasize dependence on the sample size `n`.
- In `chapters/sec-pca-theory.tex`, removed an obsolete tilde from the eigenvalue matrix in the first trace identity, consistently using the full matrix `Lambda` from `C = V Lambda V^T`.

### Figure-by-figure reconstruction record

#### Shannon Sampling Theory (5)

- **Sampling and spectral overlap** (`shannon-sampling-aliasing`). Original: `figures/1-shannon/sampling-aliasing.png`; editable source: `figures/tikz/shannon-sampling-aliasing.tex`. Unit sample spacing is used. These are exact compactly supported cubic B-spline spectra of sinc(x/a)^4, with a=4 and a=2. Replicas are spaced by 2 pi; the shaded Nyquist interval is [-pi, pi]. Ideal reconstruction retains the sum inside that interval. The right-hand result differs from the original spectrum.
- **The cardinal sine** (`shannon-sinc`). Original: `figures/1-shannon/sinc.png`; editable source: `figures/tikz/shannon-sinc.tex`. The value at zero is the continuous extension, equal to one. Every nonzero integer is a zero. The plotted normalization is sin(pi x)/(pi x), as in the theorem.
- **Cardinal B-splines of degrees zero, one, and two** (`shannon-splines`). Original: `figures/1-shannon/spline.png`; editable source: `figures/tikz/shannon-splines.tex`. Supports and heights follow the convolution definition exactly. The quadratic spline peaks at 3/4, so its integer translates are not an interpolating basis with the samples as coefficients. Endpoint values of the box follow the closed-interval convention in the text.
- **Two frequencies with identical samples** (`shannon-aliasing`). Original: `figures/1-shannon/aliasing.png`; editable source: `figures/tikz/shannon-aliasing.tex`. For unit sample spacing, cos(pi x/2) and cos(3 pi x/2) agree at every integer. Frequencies are identified modulo 2 pi; the high-frequency pair folds to the low-frequency pair. This example is periodic and is not an L2 signal on the real line.
- **Uniform scalar quantization** (`shannon-quantizer`). Original: `figures/1-shannon/quantizer.pdf`; editable source: `figures/tikz/shannon-quantizer.tex`. The horizontal axis is u/T and the vertical axis is the integer v. Filled left endpoints and open right endpoints implement v-1/2 <= u/T < v+1/2 exactly. The diagonal is the unquantized value; the discrepancy is at most one half in these normalized coordinates.

#### Fourier and Convolution (16)

- **Convolution with a box window** (`fourier-convolution`). Original: `figures/2-fourier/convolution.png`; editable source: `figures/tikz/fourier-convolution.tex`. The shaded interval is [x-epsilon,x+epsilon]. The original box is unnormalized, so its convolution is an integral over that interval, not an average. The drawn signal is illustrative; no numerical data are asserted.
- **Convolution becomes multiplication** (`fourier-convolution-fourier`). Original: `figures/2-fourier/convolution-fourier.png`; editable source: `figures/tikz/fourier-convolution-fourier.tex`. The diagram uses the chapter convention F(f*g)=f-hat times g-hat. Both vertical arrows apply the Fourier transform, making the square commute; the original reverse Fourier arrow is equivalently written in the forward direction.
- **Cardinal splines by convolution** (`fourier-splines`). Original: `figures/2-fourier/splines.pdf`; editable source: `figures/tikz/fourier-splines.tex`. **Interpretation to check.** The box, triangle, and quadratic spline are exact centered cardinal B-splines B0, B1=B0*B0, and B2=B1*B0. The original final handwritten symbol is unclear; explicit degree and convolution labels replace it.
- **Translation equivariance** (`fourier-translation-inv`). Original: `figures/2-fourier/translation-inv.png`; editable source: `figures/tikz/fourier-translation-inv.tex`. The original square is made explicit with T_tau f=f(.-tau). The bottom-right value is H(T_tau f)=T_tau(Hf); this is equivariance, called translation invariance in the chapter.
- **Fourier modes align at the Dirac comb** (`fourier-poisson-distrib`). Original: `figures/2-fourier/poisson-distrib.png`; editable source: `figures/tikz/fourier-poisson-distrib.tex`. Three real Fourier modes cos(n omega), n=1,2,3, are shown. They all equal one at multiples of 2pi. The distributional comb identity is stated separately; a finite sum is not drawn as literal Dirac masses.
- **One radix-two FFT step** (`fourier-fft`). Original: `figures/2-fourier/fft.pdf`; editable source: `figures/tikz/fourier-fft.tex`. The original butterfly wiring is replaced by an equivalent block diagram for the chapter's decimation-in-frequency recursion. Sum and difference are formed from the two halves, twiddle factors multiply the difference before its half-size FFT, and outputs are interleaved into even and odd frequencies.
- **Spatial zero padding refines frequency sampling** (`fourier-padding-spatial`). Original: `figures/2-fourier/padding-spacial.png`; editable source: `figures/tikz/fourier-padding-spatial.tex`. **Interpretation to check.** The spatial extent is T=Q/N while the physical sample step remains 1/N. The frequency interval is bounded by the unchanged Nyquist frequency pi N, and its step is 2pi/T. This corrects the handwritten reciprocal Nyquist labels. Curves are schematic and do not assert empirical sample values. The spectrum is shown by its magnitude; the two continuous curve profiles are schematic rather than an asserted exact transform pair.
- **Spectral zero padding interpolates the samples** (`fourier-padding-fourier`). Original: `figures/2-fourier/padding-fourier.png`; editable source: `figures/tikz/fourier-padding-fourier.tex`. **Interpretation to check.** The retained signed-frequency coefficients are padded with zeros, and the inverse transform is multiplied by Q/N. The illustrative trigonometric polynomial and coefficient magnitudes are consistent for N=5 and Q=15; dark markers identify the five original samples and smaller colored markers the finer grid. These are analytic example values, not measurements. Odd N avoids the even-length Nyquist splitting convention.
- **A plane wave and its normal direction** (`fourier-wave`). Original: `figures/2-fourier/wave.png`; editable source: `figures/tikz/fourier-wave.tex`. The parallel lines are constant-phase sets, normal to the wave vector omega. They indicate geometry rather than a propagation direction. Integer wave vectors apply on the periodic domain; the general figure uses omega in R squared.
- **A square with opposite edges identified** (`fourier-torus`). Original: `figures/2-fourier/torus.png`; editable source: `figures/tikz/fourier-torus.tex`. Opposite sides of the square carry matching oriented arrows. Identifying each pair produces the two-dimensional torus. The pictured doughnut is a topological visualization, not an isometric embedding of the flat quotient.
- **A periodic two-dimensional grid** (`fourier-torus-discrete`). Original: `figures/2-fourier/torus-discr.png`; editable source: `figures/tikz/fourier-torus-discrete.tex`. A four-by-four representative grid illustrates the two periodic coordinate directions. A step beyond the last displayed index wraps back to zero; the first and last displayed vertices remain distinct. The identifications occur after N1 or N2 steps, so the rectangle has no physical boundary. The grid size is illustrative.
- **Heat diffusion is Gaussian convolution** (`fourier-heat`). Original: `figures/2-fourier/heat.png`; editable source: `figures/tikz/fourier-heat.tex`. The left curves show an illustrative signal with two Fourier modes at time zero and after positive diffusion times. The right curves are Gaussian heat kernels for the same two times: increasing time broadens the kernel and lowers its peak. The normalization and variance follow the chapter heat equation.
- **Continuum and discrete Laplacian spectra** (`fourier-finite-difference-spectrum`). Original: `figures/2-fourier/finite-diff.png`; editable source: `figures/tikz/fourier-finite-difference-spectrum.tex`. The original positive curves are explicitly labeled as negative eigenvalues divided by N squared. The continuum curve is (2pi k/N) squared, reaching pi squared at Nyquist; the discrete curve is 4 sin squared(pi k/N), reaching four. Actual Laplacian eigenvalues are nonpositive.
- **A small geodesic ball on a surface** (`fourier-laplacian-surface`). Original: `figures/2-fourier/laplacian-surf.png`; editable source: `figures/tikz/fourier-laplacian-surface.tex`. The sketch is a local surface patch containing a geodesic ball centered at x. The ellipse represents the projected ball, not a claim that geodesic balls are Euclidean circles. The displayed expansion retains the necessary epsilon-squared scaling and dimension-dependent constant.
- **Colatitude and azimuth on the sphere** (`fourier-spherical-coordinates`). Original: `figures/2-fourier/spherical.png`; editable source: `figures/tikz/fourier-spherical-coordinates.tex`. **Interpretation to check.** Theta is measured from the positive north axis, matching the chapter colatitude convention. The original angle appears measured from the equatorial projection, which would instead be latitude; this ambiguity is deliberately resolved from the formula. Phi is the equatorial azimuth.
- **A weighted graph and its local Laplacian** (`fourier-weighted-graph`). Original: `figures/2-fourier/graph.png`; editable source: `figures/tikz/fourier-weighted-graph.tex`. **Interpretation to check.** The graph patch preserves the intended local adjacency picture, highlighting an edge between vertices i and j with weight w_ij. Unreadable numerical marks are not assigned values. The Laplacian sign matches the chapter nonpositive convention.

#### Shannon Coding Theory (3)

- **Entropy on the probability simplex** (`coding-entropy-extrema`). Original: `figures/1-shannon/entropy-extrema.pdf`; editable source: `figures/tikz/coding-entropy-extrema.tex`. **Interpretation to check.** The left sketch has unclear level-curve labels; it is recast as the exact binary entropy graph on p=(t,1-t). The right panel retains the three-symbol simplex, with numerically computed entropy contours at 1.05, 1.30, and 1.50 bits. Cropped prose in the scan is omitted.
- **Kraft inequality: disjoint descendant blocks** (`coding-kraft-necessity`). Original: `figures/1-shannon/kraft-ineq-1.png`; editable source: `figures/tikz/coding-kraft-necessity.tex`. The code is exactly the one specified in the caption: (0,10,110,111). The corresponding descendant blocks contain 4, 2, 1, and 1 leaves of the full depth-three tree. Highlighted circles are codeword roots; the colored leaf blocks are disjoint.
- **Kraft inequality: packing prescribed lengths** (`coding-kraft-sufficiency`). Original: `figures/1-shannon/kraft-ineq-2.png`; editable source: `figures/tikz/coding-kraft-sufficiency.tex`. For lengths (1,2,3), the blocks have 4, 2, and 1 leaves. Packing them from the left gives codewords (0,10,110), leaving 111 unused. The tree demonstrates both alignment and the inequality 7/8 <= 1; it does not assume equality.

#### Wavelets (7)

- **Orthogonalizing a cardinal spline** (`wavelets-spline-orthogonalization`). Original: `figures/wavelets/spline-ortho.pdf`; editable source: `figures/tikz/wavelets-spline-orthogonalization.tex`. The formula identifies the handwritten example as the orthogonal linear spline. Its triangular generator has compact support, while the orthogonalized function has alternating exponentially decaying tails. The replacement uses the correct Fourier multiplier sqrt(3/(2+cos omega)); its cosine coefficients were evaluated by 8192-point midpoint quadrature to draw the illustrative profile. Four small degree examples recover the original spline-family sketches.
- **Sharp and smooth frequency partitions** (`wavelets-shannon-spline-bands`). Original: `figures/wavelets/spline-vs-shannon.pdf`; editable source: `figures/tikz/wavelets-shannon-spline-bands.tex`. **Interpretation to check.** The Shannon row shows the exact low-pass band |omega|<=pi and the next two dyadic wavelet bands. The spline row deliberately shows schematic smooth overlap, not the measured spectrum of a named spline order. Only a finite frequency window is drawn; the partition identity refers to all scales. Cropped neighboring prose in the scan is omitted. The displayed normalized outer-band energy is twice the energy of psi at scale minus one, compensating for its L2 dilation factor.
- **Wrapping a function around the unit interval** (`wavelets-periodization`). Original: `figures/wavelets/periodize-1.pdf`; editable source: `figures/tikz/wavelets-periodization.tex`. **Interpretation to check.** The drawing uses an illustrative localized function f, its left-shifted copy, and their periodic sum on [0,1]. It makes the matching boundary values explicit. The profile is not presented as a particular orthogonal wavelet; the construction applies to each basis function individually.
- **One representative of each periodic translate** (`wavelets-periodic-translates`). Original: `figures/wavelets/periodize-2.pdf`; editable source: `figures/tikz/wavelets-periodic-translates.tex`. **Interpretation to check.** The selected translation centers lie in the half-open unit interval. The example j=-2 has four distinct centers; the center at one is the same periodic translate as zero. Smooth pulse profiles are schematic and do not assert that this particular shape generates an orthogonal wavelet basis.
- **Haar filter coefficients from overlap integrals** (`wavelets-haar-inner-products`). Original: `figures/wavelets/haar-inner.pdf`; editable source: `figures/tikz/wavelets-haar-inner-products.tex`. The original inner products are +1 for the scaling function and +1 or -1 for the wavelet. The chapter filter definition adds a factor 1/sqrt(2); the resulting h=(1,1)/sqrt(2) and g=(1,-1)/sqrt(2) are therefore explicitly distinguished from the raw overlap integrals.
- **Complementary low-pass and high-pass energies** (`wavelets-filter-constraints`). Original: `figures/wavelets/filter-constr.pdf`; editable source: `figures/tikz/wavelets-filter-constraints.tex`. The solid and dashed curves are the exact Haar examples |h-hat| squared=2 cos squared(omega/2) and its pi shift=2 sin squared(omega/2). Their sum is two, and the shifted energy equals that of the QMF high-pass filter. The nearby cropped text in the scan is replaced by the complete relevant equality.
- **Vanishing moments and flat filter zeros** (`wavelets-vanishing-moments`). Original: `figures/wavelets/vanmoments.pdf`; editable source: `figures/tikz/wavelets-vanishing-moments.tex`. The drawing separates the low-pass zero at pi from the high-pass zero at zero, avoiding the ambiguous overlaid labels in the scan. The curves are magnitudes for a two-moment QMF example; the derivative statement concerns the smooth filter symbols themselves, including their phases. The displayed zeros have order two, while the formula states the general p-moment equivalence.

#### Linear and Nonlinear Approximation (1)

- **Approximation rates across signal models** (`approximation-rates`). Original: `figures/approximation/approx-comparaison.png`; editable source: `figures/tikz/approximation-rates.tex`. **Interpretation to check.** The sketches are reorganized into a rate table. L and NL mean linear and nonlinear approximation. The one-dimensional piecewise-smooth rates are shown for alpha at least one; wavelets require the stated regularity, moments, and boundary assumptions. The BV class also has the chapter's uniform bound. The cartoon column assumes C2 pieces and contours, illustrated by a smooth closed curve. The curvelet rate retains the logarithmic factor given in the chapter, which the handwritten summary omitted. The entries are upper rates, not assertions that every signal attains them.

#### Inverse Problems (2)

- **Capped inversion and Tikhonov filtering** (`inverse-problems-variance-bound`). Original: `figures/inverse-problems/bound-variance.png`; editable source: `figures/tikz/inverse-problems-variance-bound.tex`. **Interpretation to check.** The left plateau is interpreted explicitly as the capped inverse 1/max(sigma,tau), rather than truncated SVD. A separate cutoff tau avoids confusing its parameter with the squared-scale Tikhonov parameter lambda. The right maximum is attained at sigma=sqrt(lambda) with height 1/(2sqrt(lambda)). Curves are rescaled schematically.
- **The source-order dependence of the bias bound** (`inverse-problems-bias-bound`). Original: `figures/inverse-problems/bound-bias.png`; editable source: `figures/tikz/inverse-problems-bias-bound.tex`. Corrects the handwritten maximizer: sigma*=sqrt(beta*lambda/(2-beta)), including the missing factor beta. For beta=2 the supremum lambda is approached asymptotically. For beta>2 the supremum over all sigma>=0 is infinite; this does not say that the bias on a fixed bounded operator spectrum is infinite. Curves show representative beta values with rescaled axes.

#### Sparse Regularization (5)

- **Hard thresholding as a scalar minimization** (`sparse-regularization-hard-objective`). Original: `figures/sparse-reg/var-l0.pdf`; editable source: `figures/tikz/sparse-regularization-hard-objective.tex`. Uses exactly the caption normalization F(x)=(x-y)^2+T^2*1_{x!=0}. The objective has an isolated value y^2 at zero and minimum T^2 among nonzero candidates, attained at y. The displayed example has y>T, and the decision rule also states the equality case.
- **Soft thresholding: lambda=0** (`sparse-regularization-soft-objective-1`). Original: `figures/sparse-reg/var-l1-1.pdf`; editable source: `figures/tikz/sparse-regularization-soft-objective-1.tex`. Reconstructs the exact piecewise-quadratic objective one-half*(x-y)^2+lambda*abs(x). The minimizer is max(y-lambda,0). At lambda=y the right derivative at zero is zero; at lambda>y it is positive, while the left derivative remains negative. The curves are illustrative numerical instances of the symbolic cases.
- **Soft thresholding: 0<lambda<y** (`sparse-regularization-soft-objective-2`). Original: `figures/sparse-reg/var-l1-2.pdf`; editable source: `figures/tikz/sparse-regularization-soft-objective-2.tex`. Reconstructs the exact piecewise-quadratic objective one-half*(x-y)^2+lambda*abs(x). The minimizer is max(y-lambda,0). At lambda=y the right derivative at zero is zero; at lambda>y it is positive, while the left derivative remains negative. The curves are illustrative numerical instances of the symbolic cases.
- **Soft thresholding: lambda=y** (`sparse-regularization-soft-objective-3`). Original: `figures/sparse-reg/var-l1-3.pdf`; editable source: `figures/tikz/sparse-regularization-soft-objective-3.tex`. Reconstructs the exact piecewise-quadratic objective one-half*(x-y)^2+lambda*abs(x). The minimizer is max(y-lambda,0). At lambda=y the right derivative at zero is zero; at lambda>y it is positive, while the left derivative remains negative. The curves are illustrative numerical instances of the symbolic cases.
- **Soft thresholding: lambda>y** (`sparse-regularization-soft-objective-4`). Original: `figures/sparse-reg/var-l1-4.pdf`; editable source: `figures/tikz/sparse-regularization-soft-objective-4.tex`. Reconstructs the exact piecewise-quadratic objective one-half*(x-y)^2+lambda*abs(x). The minimizer is max(y-lambda,0). At lambda=y the right derivative at zero is zero; at lambda>y it is positive, while the left derivative remains negative. The curves are illustrative numerical instances of the symbolic cases.

#### Theory of Sparse Regularization (6)

- **Interior points and failure of exact recovery** (`sparse-theory-polytope-proof`). Original: `figures/sparse-theory/proofs/polytope-proof-1.png`; editable source: `figures/tikz/sparse-theory-polytope-proof.tex`. Uses the normalized polytope AB and q=Ax0/r, r=||x0||1, consistently with the corrected proof. Replaces the original ambiguous label z on an observation-space boundary point by (1+t)q. Both panels illustrate failure of optimality through interior membership; no uniqueness claim is implied.
- **Removing a dependent active coefficient** (`sparse-theory-injectivity`). Original: `figures/sparse-theory/proofs/injectivity.png`; editable source: `figures/tikz/sparse-theory-injectivity.tex`. The original draws affine coefficient trajectories. Three representative coordinates are drawn with a marked sign-preserving interval. At its left endpoint the first coefficient vanishes. The constant-fidelity and affine-l1 statements follow the chapter proof; the drawn slopes also sum to zero within the positive orthant.
- **Bregman divergence as a vertical gap** (`sparse-theory-bregman`). Original: `figures/sparse-theory/proofs/bregman.png`; editable source: `figures/tikz/sparse-theory-bregman.tex`. The left panel uses a strictly convex quadratic and its tangent; the positive vertical gap is the Bregman divergence. The right panel uses J=abs(x), with x0 and x on the same positive affine branch, so the divergence vanishes even when x differs from x0.
- **A strict subgradient margin controls the error** (`sparse-theory-bregman-l1`). Original: `figures/sparse-theory/proofs/bregman-l1.png`; editable source: `figures/tikz/sparse-theory-bregman-l1.tex`. Explicitly sets x0=0, as required on nonsaturated coordinates for J=abs(x). For abs(eta)<1 the supporting line eta*x lies strictly below abs(x) away from zero. The lower bound is D_eta(x|0)>= (1-abs(eta))*abs(x), not a norm bound on saturated coordinates.
- **Parameter region for sign consistency** (`sparse-theory-recovery-region`). Original: `figures/sparse-theory/proofs/recov-zone.png`; editable source: `figures/tikz/sparse-theory-recovery-region.tex`. Reconstructs the two strict inequalities and the contained simpler triangle. The vertical cutoff is lambda_max=TR/[K(R+S)], using R=KL+1. Shading indicates interiors; boundary lines are not included in the sufficient conditions. The axes follow the caption: horizontal lambda, vertical delta.
- **Convolution of a discrete measure** (`sparse-theory-convolution-spikes`). Original: `figures/sparse-theory/proofs/spikes.png`; editable source: `figures/tikz/sparse-theory-convolution-spikes.tex`. Two shifted kernels and two point masses are illustrative; their locations are reconstruction-grid points z_j and z_i. The output is their weighted sum, written as the continuous operator applied to m_x and as the discretized dictionary product Ax. The diagram distinguishes the index-space weights x from the physical locations z.

#### Basics of Machine Learning (13)

- **PCA directions and projected variance** (`ml-pca-variance`). Original: `figures/ml/variance.png`; editable source: `figures/tikz/ml-pca-variance.tex`. This centered illustrative cloud has exact reflection symmetries in its principal coordinates. The solid ellipse has semiaxes equal to the empirical standard deviations, using the chapter normalization 1/n. Variance is sigma squared, whereas the indicated geometric lengths are sigma. The point locations are illustrative, not recovered measurements.
- **Nearest-centroid assignment and Voronoi cells** (`ml-voronoi`). Original: `figures/ml/voronoi.png`; editable source: `figures/tikz/ml-voronoi.tex`. Cell edges are actual perpendicular bisectors of the displayed centroids. Their meeting point is equidistant from all three. Colored sample points lie in the corresponding nearest-centroid cells. The diagram illustrates an assignment step and does not assert that these centroids are already cluster means.
- **A joint probability model for input and output** (`ml-joint-model`). Original: `figures/ml/proba-model.pdf`; editable source: `figures/tikz/ml-joint-model.tex`. Contours represent an illustrative bivariate Gaussian density, not an asserted fit to the original hand-drawn points. The dashed projections identify the two coordinates of one observation. The general discussion in the chapter is not restricted to Gaussian models.
- **Conditional expectation as the mean of a slice** (`ml-conditional-expectation`). Original: `figures/ml/cond-expect.pdf`; editable source: `figures/tikz/ml-conditional-expectation.tex`. An illustrative model Y=m(X)+epsilon with centered Gaussian noise is used. The left curves show constant conditional-density levels around m(x); the right density is normalized at x=x0. Its mean is m(x0), which need not be linear in x. No distributional assumption is added to the theorem.
- **An affine fit and a nonlinear regression function** (`ml-linear-fit`). Original: `figures/ml/linear-fit.png`; editable source: `figures/tikz/ml-linear-fit.tex`. The red line is the least-squares affine fit to the actual displayed illustrative points. The dashed curve is the model conditional mean used to construct those points. Residuals are vertical output errors, not orthogonal distances to the line.
- **Four nearest neighbors of a query point** (`ml-nearest-neighbors`). Original: `figures/ml/clustering.png`; editable source: `figures/tikz/ml-nearest-neighbors.tex`. The arrows identify the four nearest observations in increasing distance order. The dashed circle has radius equal to the fourth distance. All remaining points are farther away. This is a nearest-neighbor query, not the centroid assignment used in k-means.
- **Logistic probability and transition width** (`ml-logistic-one-dimension`). Original: `figures/ml/logistic-1d.png`; editable source: `figures/tikz/ml-logistic-one-dimension.tex`. The increasing sign convention follows the current equation. Two positive slopes show the same probability-one-half boundary at x=0; the larger slope gives a narrower transition. The vertical axis is a probability in [0,1], not a class label in {-1,+1}.
- **The same decision boundary at two parameter norms** (`ml-logistic-two-dimensions`). Original: `figures/ml/logistic-2d.png`; editable source: `figures/tikz/ml-logistic-two-dimensions.tex`. Both panels use the same unit normal v and beta=t v. The shaded bands delimit probabilities from 0.1 to 0.9, so their normal width is 2 log(9)/t. Changing t leaves the decision boundary unchanged; it does not change whether a data set is linearly separable.
- **Rows are observations; columns are features** (`ml-pca-data-matrix`). Original: `figures/ml/pca-maths/pca-maths-1.pdf`; editable source: `figures/tikz/ml-pca-data-matrix.tex`. The reconstruction follows the current n-observations, p-features convention. An observation x_i is a column vector in R^p; the corresponding matrix row is its transpose. Centering is explicit, and the covariance is X transpose X divided by n.
- **Linear compression followed by reconstruction** (`ml-pca-compression-reconstruction`). Original: `figures/ml/pca-maths/pca-maths-2.pdf`; editable source: `figures/tikz/ml-pca-compression-reconstruction.tex`. Both R and S have p rows and k columns. Compression applies S transpose; reconstruction applies R, so the product is a p-by-p map of rank at most k. Orthogonality is not assumed for arbitrary R and S; it is justified later when reducing to orthogonal projections.
- **The polytope bounding the PCA objective** (`ml-pca-linear-program`). Original: `figures/ml/pca-maths/pca-maths-3.pdf`; editable source: `figures/tikz/ml-pca-linear-program.tex`. The left feasible set is the triangle for p=2, k=1. The right is the unit cube with the corner (1,1,1) cut off for p=3, k=2. The marked maximizers are (1,0) and (1,1,0), respectively, for ordered nonnegative eigenvalues. The three-dimensional shaded face is sum beta=2.
- **Row norms of the rotated basis** (`ml-pca-row-norms`). Original: `figures/ml/pca-maths/pca-maths-4.pdf`; editable source: `figures/tikz/ml-pca-row-norms.tex`. The matrix B has p rows and k orthonormal columns. The vector b_i is in R^k and appears transposed as a row. Its squared row norms sum to k, the squared Frobenius norm of B. These dimensions follow the current proof rather than ambiguous letters in the scan.
- **Completing the columns to an orthogonal matrix** (`ml-pca-orthogonal-extension`). Original: `figures/ml/pca-maths/pca-maths-5.pdf`; editable source: `figures/tikz/ml-pca-orthogonal-extension.tex`. The first k columns are B; the remaining p-k columns form B perpendicular. Every full row of the square orthogonal matrix has norm one. Restricting a row to the first k entries cannot increase its norm. When k=p, the complementary block is empty.

#### Optimization \& Machine Learning: Smooth Optimization (24)

- **Convexity of functions and sets** (`smooth-convex-examples`). Original: `figures/optim-smooth/convex-examples.png`; editable source: `figures/tikz/smooth-convex-examples.tex`. Red segments distinguish failure of convexity, non-strict convexity with boundary segments, and strict convexity. The middle function has a flat interval, and the middle set has flat boundary segments.
- **Linear regression** (`smooth-linear-regression`). Original: `figures/optim-smooth/ml-exemp-1.pdf`; editable source: `figures/tikz/smooth-linear-regression.tex`. The line represents the scalar-feature illustration of y=〈a,x〉; point positions are schematic and are not a numerical dataset.
- **A linear classifier** (`smooth-linear-classifier`). Original: `figures/optim-smooth/ml-exemp-2.pdf`; editable source: `figures/tikz/smooth-linear-classifier.tex`. Positive labels are placed on the positive-score side and negative labels on the negative-score side. The normal vector is perpendicular to the separating line.
- **Classification losses** (`smooth-classification-losses`). Original: `figures/optim-smooth/classif-loss.pdf`; editable source: `figures/tikz/smooth-classification-losses.tex`. The sketch labels a logistic loss divided by log(2), whereas this section defines the unscaled log(1+exp(u)); the reconstruction follows the section. The hinge is max(0,1+u). Zero margin is counted as an error. The unscaled logistic loss is not claimed to majorize the 0–1 loss.
- **Existence and uniqueness of minimizers** (`smooth-minimizer-examples`). Original: `figures/optim-smooth/uniqueness.pdf`; editable source: `figures/tikz/smooth-minimizer-examples.tex`. Five representative functions retain the original cases. The two no-minimizer examples distinguish an unbounded infimum from a finite unattained infimum. Vertical offsets and coordinate scales are schematic.
- **Coercive least squares** (`smooth-least-squares-coercive`). Original: `figures/optim-smooth/least-square-1.pdf`; editable source: `figures/tikz/smooth-least-squares-coercive.tex`. The bowl depicts the full-column-rank case. Its height may include a nonzero residual at the minimizer; no exact fit is assumed.
- **Least squares with a nontrivial kernel** (`smooth-least-squares-flat`). Original: `figures/optim-smooth/least-square-2.pdf`; editable source: `figures/tikz/smooth-least-squares-flat.tex`. The trough is flat exactly along the kernel direction. The highlighted affine line represents x-star+ker(A), rather than an isolated minimizer.
- **Failure of the secant inequality** (`smooth-secant-nonconvex`). Original: `figures/optim-smooth/cvx-vs-noncvx-1.pdf`; editable source: `figures/tikz/smooth-secant-nonconvex.tex`. The interior point is compared with the secant at the same abscissa. The convex example also includes a supporting tangent; the nonconvex example explicitly violates the secant inequality.
- **Convexity and the secant inequality** (`smooth-secant-convex`). Original: `figures/optim-smooth/cvx-vs-noncvx-2.pdf`; editable source: `figures/tikz/smooth-secant-convex.tex`. The interior point is compared with the secant at the same abscissa. The convex example also includes a supporting tangent; the nonconvex example explicitly violates the secant inequality.
- **Strict convexity** (`smooth-strict-convexity`). Original: `figures/optim-smooth/strictly-cvx-2.pdf`; editable source: `figures/tikz/smooth-strict-convexity.tex`. The reconstruction follows the definitions and notation in the accompanying section.
- **Convexity with an affine interval** (`smooth-nonstrict-convexity`). Original: `figures/optim-smooth/strictly-cvx-1.pdf`; editable source: `figures/tikz/smooth-nonstrict-convexity.tex`. The highlighted interval is exactly affine, so equality holds in the secant inequality there. The surrounding curve is chosen with nondecreasing slope.
- **The gradient is a vector field** (`smooth-gradient-field`). Original: `figures/optim-smooth/gradient-vector-field.pdf`; editable source: `figures/tikz/smooth-gradient-field.tex`. The radial field is shown for f(x)=||x||²/2, so the gradient points outward and vanishes at the minimizer. The accompanying bowl is schematic; its axes are not numerically calibrated.
- **Stationarity at local extrema** (`smooth-stationary-extrema`). Original: `figures/optim-smooth/first-order-1.pdf`; editable source: `figures/tikz/smooth-stationary-extrema.tex`. The representative double-well function has two minima and a maximum, all with horizontal tangents. Stationarity alone does not distinguish these cases.
- **A stationary point need not be an extremum** (`smooth-stationary-inflection`). Original: `figures/optim-smooth/first-order-2.pdf`; editable source: `figures/tikz/smooth-stationary-inflection.tex`. The reconstruction uses f(x)=x³ as in the accompanying text and labels the point as a stationary inflection point. The old caption calls the one-dimensional example a saddle point; the more precise terminology is used here.
- **Stationarity for a convex objective** (`smooth-stationary-convex`). Original: `figures/optim-smooth/first-order-3.pdf`; editable source: `figures/tikz/smooth-stationary-convex.tex`. The tangent plane touches the bowl at its minimum and supports the graph from below. The plane is horizontal in the depicted three-dimensional coordinates.
- **PCA and least-squares level sets** (`smooth-pca-quadratic-geometry`). Original: `figures/optim-smooth/link-pca.pdf`; editable source: `figures/tikz/smooth-pca-quadratic-geometry.tex`. Both panels share the eigenvectors u₁,u₂. The twelve reflection-symmetric points are centered and have empirical covariance eigenvalues λ₁=1 and λ₂=0.325². Their covariance ellipses and the quadratic level sets are derived from these same eigenvalues, with reciprocal axis ratios. The caption inside the drawing specifies the common multiplier used for each family of semiaxes; axis labels follow the covariance normalization C/n in the text.
- **A first-order Taylor approximation** (`smooth-taylor-line`). Original: `figures/optim-smooth/taylor-exp-1.pdf`; editable source: `figures/tikz/smooth-taylor-line.tex`. The drawn line is the exact tangent to the representative quadratic at x; the vertical gap at z shows the Taylor remainder.
- **A tangent plane** (`smooth-taylor-plane`). Original: `figures/optim-smooth/taylor-exp-2.pdf`; editable source: `figures/tikz/smooth-taylor-plane.tex`. The tangent plane is computed for the displayed representative paraboloid and touches it at the marked point. The graph and its supporting plane are distinct objects.
- **The gradient is orthogonal to a level curve** (`smooth-level-set-tangent`). Original: `figures/optim-smooth/level-sets-1.pdf`; editable source: `figures/tikz/smooth-level-set-tangent.tex`. At the leftmost point of an elliptical level curve, the tangent is vertical and the outward gradient is horizontal. The right-angle marker records orthogonality, not a descent step.
- **Gradients and nested level sets** (`smooth-nested-level-sets`). Original: `figures/optim-smooth/level-sets-2.pdf`; editable source: `figures/tikz/smooth-nested-level-sets.tex`. The arrows are true normals to elliptical level curves, rather than radial arrows except on the principal axes. All level values increase outward.
- **Step-size effects in gradient descent** (`smooth-gradient-step-size`). Original: `figures/optim-smooth/grad-desc-1.pdf`; editable source: `figures/tikz/smooth-gradient-step-size.tex`. Both trajectories are exact gradient-descent iterates for the displayed quadratic. The larger step oscillates but remains stable because both step sizes lie in (0,2/L), with L=5. The sketch is schematic; no numerical values were legible in it.
- **Exact line search** (`smooth-exact-line-search`). Original: `figures/optim-smooth/grad-desc-2.pdf`; editable source: `figures/tikz/smooth-exact-line-search.tex`. The path uses exact line search for a positive definite quadratic. The right-angle marker is computed from consecutive actual steps; it is not merely a schematic zigzag.
- **The optimal quadratic step size** (`smooth-quadratic-contraction`). Original: `figures/optim-smooth/proof-quadr.pdf`; editable source: `figures/tikz/smooth-quadratic-contraction.tex`. The highlighted curve is the maximum of the two absolute-value branches. Its minimum occurs at 2/(L+μ), and h(τ)<1 exactly for 0<τ<2/L when 0<μ≤L. A representative ratio L/μ=3 determines the drawing; the labels remain symbolic.
- **Quadratic upper and lower bounds** (`smooth-curvature-bounds`). Original: `figures/optim-smooth/up-low-bounds.pdf`; editable source: `figures/tikz/smooth-curvature-bounds.tex`. Both bounding quadratics have the same value and gradient as f at x. The lower curvature μ and upper curvature L are ordered correctly on both sides of the expansion point.

#### Optimization \& Machine Learning: Advanced Topics (2)

- **An unbiased stochastic gradient** (`advanced-unbiased-gradient`). Original: `figures/ml/sgd/unbiased-grad.png`; editable source: `figures/tikz/advanced-unbiased-gradient.tex`. The red dashed vector is the exact arithmetic mean of the three drawn individual gradient vectors. This preserves unbiasedness geometrically, rather than placing the mean arrow arbitrarily.
- **Stochastic gradient trajectories** (`advanced-sgd-trajectory`). Original: `figures/ml/sgd/sgd-schematic.png`; editable source: `figures/tikz/advanced-sgd-trajectory.tex`. **Interpretation to check.** The shrinking ellipses and path retain the qualitative intent of the sketch. They are illustrative neighborhoods, not proved confidence regions or a simulated dataset. Convergence requires the assumptions and step-size conditions developed in the section; it is not asserted for arbitrary constant steps.

#### Deep Learning (2)

- **A fully connected feedforward network** (`deep-learning-fully-connected`). Original: `figures/deep-learning/fc.png`; editable source: `figures/tikz/deep-learning-fully-connected.tex`. Retains the original three inputs, two hidden layers of four units, and three outputs. Every adjacent-layer pair is connected; line crossings are not nodes. The output is labelled as scores rather than class labels because a separate logistic map converts scores into probabilities in the chapter.
- **Feature maps in a convolutional network** (`deep-learning-convolutional`). Original: `figures/deep-learning/cnn.png`; editable source: `figures/tikz/deep-learning-convolutional.tex`. **Interpretation to check.** The low-resolution original dimension labels are replaced by the chapter notation: n_l=bar n_l*d_l and RGB input d_0=3. The five tensor blocks show an input, convolutional features, their downsampled version, a second feature stack, and a second downsampling. Depth illustrates channel count and face size illustrates spatial resolution; exact numerical widths or filter sizes cannot be recovered and are not invented.

#### Convex Analysis (4)

- **Convexity of functions and sets** (`convex-analysis-convex-examples`). Original: `figures/convexity/convex-examples.png`; editable source: `figures/tikz/convex-analysis-convex-examples.tex`. The same original also appears in the smooth-optimization chapter. This reconstruction is shared in design with that chapter: segments exhibit failure of convexity, an affine or flat boundary segment, and strict convexity. The epigraph is shaded above each function. The functions and sets are schematic examples rather than named formulas from the source.
- **Supporting affine functions at a kink** (`convex-analysis-subdifferential`). Original: `figures/convexity/subdifferential.png`; editable source: `figures/tikz/convex-analysis-subdifferential.tex`. The original fan of supporting lines is reconstructed at a nonsmooth point of a convex function. Every drawn slope belongs to the interval between the one-sided derivatives. The diagram labels the supporting expression with a separate evaluation variable z, avoiding confusion with the base point x.
- **Graphs of two subdifferentials** (`convex-analysis-subdiff-l1`). Original: `figures/convexity/subdiff-l1.png`; editable source: `figures/tikz/convex-analysis-subdiff-l1.tex`. **Interpretation to check.** At each kink the entire closed interval of intermediate slopes is drawn, not just its endpoints. The second function uses representative slopes m1<m2<0<m3 at knots a<0 and 0; the original does not specify numerical slopes. The lower panels are set-valued graphs, so their vertical segments are essential.
- **Normal cones at boundary, interior, and exterior points** (`convex-analysis-normal-cones`). Original: `figures/convexity/normal-cone.png`; editable source: `figures/tikz/convex-analysis-normal-cones.tex`. The convex body is rendered schematically as a polygon so its normal directions are exact. Edge interiors have outward normal rays, corners have cones, interior points have only the zero normal, and exterior points have an empty normal cone. The arrows depict directions and are translated to their base points; they do not bound the cone lengths.

#### Nonsmooth Convex Optimization (1)

- **Projection and proximal maps** (`optim-nonsmooth-prox-projection`). Original: `figures/convexity/prox-proj.png`; editable source: `figures/tikz/optim-nonsmooth-prox-projection.tex`. The left panel uses a quarter disk; all displayed input-to-projection segments lie in the appropriate outward normal cones. The right panel uses nested quadratic sublevel sets and one exact proximal pair for f(u,v)=(u^2+2v^2)/2. It states that the relevant sublevel set depends on the input through z=prox_f(x); a single fixed set is not a projection representation of the entire proximal map.


### Verification of the figure review edition

- Rebuilt and published `FundationsDataScience.pdf` (**339 pages**) and all **19 standalone chapter PDFs** with **zero final LaTeX/BibTeX warnings, missing characters, or overfull/underfull boxes**.
- Verified all **91 standalone TikZ PDFs** are single-page vector documents, have embedded fonts, contain no raster image objects, and have no text outside their page bounds. Their PDFs are current with respect to their editable sources.
- Matched all **91 comparison pages** in the book to the manifests and found exactly the corresponding **91 comparisons** across the standalone chapters. Every comparison has both panel labels and the expected landscape orientation.
- Verified embedded fonts throughout all 20 published documents, **1,857 internal links**, **156 links from chapters to the full book**, and **725 shared label numbers**, with no unresolved destinations or numbering mismatches.
- Rendered and visually reviewed the comparison pages, 33 portrait transition pages, and 60 standalone chapter samples. Refined overlapping sampling annotations, the quadratic-spline plotting expression, the aliasing-panel spacing, the conditional-mean alignment, the FFT connector label, and the PCA dimension labels.
- Independently cross-checked the mathematical reconstructions, including Fourier normalization and padding, Haar coefficients and spectral bands, Tikhonov extrema, approximation assumptions, Bregman and normal-cone geometry, matrix dimensions, line search, curvature bounds, and stochastic-gradient averaging. The illustrated PCA clouds now use exact covariance constructions rather than unmatched arbitrary points.
- Build and validation evidence is retained in `build/build-report.json`, `build/tikz-pdf-audit.json`, `build/pdf-audit.json`, `build/label-number-audit.json`, `build/tikz-reviews/`, and `build/qa/tikz/`. The comparison notes retain the qualifications needed to distinguish exact mathematical examples from uncertain details of the originals.

## Third editorial and mathematical polishing pass

Completed another full reading of the preface and all 19 active chapters, including their section inputs, captions, definitions, statements, and proofs. This pass revised 23 of the 26 active source files; the preface and the two optimization chapter wrappers were reviewed and retained. The existing visual style, original figures, editable TikZ sources, and all 91 original/TikZ comparisons were preserved. Existing labels and citation commands were checked against the start-of-pass sources.

The entries below record the actual changes, including the assumptions and notation needed to make the exposition mathematically precise. Earlier correction records remain unchanged.

### Shannon sampling — `chapters/shannon.tex`

- **Continuous and discrete signals:** connected the continuous model directly to the analysis of signal-processing methods and to physical sensors. Split the opening into two short, connected sentences and avoided repeating “examples.”
- **Translation-invariant sampling:** rewrote the opening as a sequence of operations—extend, convolve, sample—so that the introductory clause no longer dangles from “sampler.” The sampling formula and all boundary conventions are unchanged.
- **Poisson formula:** stated the concrete correspondence between sampling in the signal domain and periodization in the Fourier domain. Removed a redundant formula name and clarified the introduction to the commutative diagram.
- **Spline interpolation:** put the interpolation conditions before the resulting linear system. Qualified the statement that coefficients equal samples: this is an identity for general sample data only for degree 0 or 1; special higher-degree data may also happen to give equality. No interpolation kernel or formula changed.
- **Quantization:** separated the loss introduced by quantization from the exact preservation of quantized symbols by source coding and decoding.

All 22 displayed blocks and every inline formula were preserved, apart from moving an unchanged interpolation condition earlier within its sentence. No mathematical theorem was altered.

### Fourier and convolution — `chapters/fourier.tex`

#### Editorial changes

- Replaced “expressing their input in a basis” by “expanding the input signal in a basis.” Clarified that the conventional Legendre normalization differs from unit L² norm.
- Described normalized Haar integration without defining `dt` in terms of itself; introduced Young's inequality as a complete statement. Replaced the vague “smallest space” by explicit containment within the Lᵖ scale.
- Explained translation invariance as applying the same rule at every location. Used a consistent operator name in the proof.
- Standardized the section title to “Discrete Orthonormal Bases,” identified the displayed distance as a **squared** Euclidean distance, and replaced “decomposition … in such an ortho-basis” with a direct expansion statement.
- Shortened the zero-padding introduction, removed “discrete finite,” and supplied its missing punctuation. Identified the spectrum as the object being padded.
- Improved the tensor-product introduction, named Fubini's theorem and the rectangular partial sums, and removed “We detail … but this generalizes.”
- Replaced the conversational warning about noncommutative convolution by a direct statement that convolution need not commute.

#### Mathematical and notation clarifications

1. **Shannon basis example:** specified the subspace of `L²(R)` and the continuous representative used for samples. Called the sinc translates orthonormal. The preceding chapter proves the pointwise formula under decay assumptions stronger than general bandlimited L² membership, so the attribution now says that explicitly. The L²-basis and pointwise expansion statements themselves are unchanged in content.
2. **Translation-invariant operator proof:** changed its local name `T` to `H`, matching the next display, and its stated codomain from `C⁰` to `C_b` with the sup norm, matching the theorem. Young's inequality supplies the required bounded map.
3. **FFT recursion:** changed the types of `f_e` and `f_o` from `R^{N'}` to `C^{N'}`. The DFT is defined on `C^N`; even a real initial input produces complex recursive inputs in the odd branch after multiplication by twiddle factors. The sums, factors, interleaving, and complexity formulas were preserved.
4. **Finite/continuous comparison:** corrected the referenced circle from `R/Z` to `R/(2πZ)`, the normalization used by the preceding translation-invariance results.
5. **Spatial padding endpoint:** continuity of a function on R supported in `[0,1]` gives `f(1)=0`; therefore padded samples vanish for `n≥N`, including the endpoint. The Riemann-sum estimate is unchanged. Spectral padding now explicitly assumes `Q≥N`.
6. **Plane-wave geometry:** described constant-phase hyperplanes and excluded the zero wave vector from the orthogonality description. A fixed complex value of a periodic exponential can correspond to a union of parallel hyperplanes, rather than one hyperplane.
7. **Multidimensional sampling:** explicitly uses the continuous representative of the bandlimited function. Decay and L² Fourier support alone do not determine values on null sets, so a statement about samples and every spatial point needs this convention.
8. **Laplacian on the torus:** corrected `T` to `T^d` when the frequency index lies in `Z^d`.
9. **Compact-group paragraph:** restricted the statement about an infinite irreducible family to **infinite** compact groups. Finite groups are compact too and have finitely many irreducibles.
10. **Graph matrix and gradient:** made `W` an N×N matrix, with zero weights on nonedges. Restricted the gradient components to `(i,j)∈E`; the old `i<j` indexing included all vertex pairs while declaring a P-dimensional codomain with P equal to the number of edges. The quadratic form and Laplacian identities are unchanged.

All 100 displayed blocks, all Fourier equation labels, and all Fourier reference uses are identical to the baseline. The mathematical changes above are in assumptions, notation explanations, or inline types.

### Shannon coding — `chapters/shannon-coding.tex`

- Distinguished signal samples from a continuous signal in the alphabet example.
- Explained that tree traversal starts at the root and that codewords correspond to leaves.
- Made the cross-entropy/likelihood equivalence explicit: the data distribution is empirical, and minimizing its cross-entropy against a model is equivalent to maximizing likelihood and minimizing relative entropy. The empirical entropy is the constant difference.
- Explained the equality conditions in `L−H=KL−log₂Z`: equality requires both `Z=1` and `q=p`, which yields ideal lengths `−log₂p_k`.
- Identified the entropy rate as an asymptotic lower bound on the **expected number of bits per symbol**, and explained that modeling dependence reduces expected length.
- Replaced “Integer rounding can mask this excess” by the precise statement that the difference between expected lengths of rounded prefix codes need not equal the divergence.
- Restricted the entropy-preservation statement to invertible transformations of **discrete symbols**. This matters because the chapter also mentions differential entropy; an invertible continuous change of variables generally changes differential entropy by an expected log-Jacobian term.

All 15 display blocks, all labels and reference uses remain identical. The only added inline equations are the already-used equality conditions `Z=1` and `q=p`.

### Wavelets — `chapters/wavelets.tex`

#### Spectral orthogonalization: corrected statement and completed proof

The old statement assumed a vaguely “sufficiently regular” and “sufficiently fast” generator, asserted the spectral identity at every frequency, and invoked the earlier pointwise Poisson proposition. That proposition has compact spectral support and cubic spatial decay assumptions unavailable for general spline generators. The old proof also established orthonormality of the normalized family without establishing that its closed span equals the original translate span.

The final proposition is the natural L² statement:

- Let θ belong to L²(R).
- Its integer translates are orthonormal iff `A(ω)=Σ_k|θ̂(ω−2πk)|²=1` **almost everywhere**.
- If `0<a≤A≤b<∞` almost everywhere, define `φ̂=θ̂/√A`; its translates are an orthonormal basis of the original closed translate span.

Compact derivation, using the chapter's Fourier and inner-product conventions:

1. Tonelli and Plancherel give `A∈L¹(T)` and `(2π)⁻¹∫₀^{2π}A=‖θ‖²₂`.
2. Partitioning R into periods gives the nth Fourier coefficient
   `(2π)⁻¹∫₀^{2π}A(ω)e^(−inω)dω = (2π)⁻¹∫R|θ̂(ω)|²e^(−inω)dω = ⟨θ,θ(·+n)⟩`.
   The **plus** translation is correct because the inner product is linear in its first argument.
3. Thus these coefficients are `δ_{n,0}` iff the translates are orthonormal. Uniqueness of Fourier coefficients for L¹ functions gives `A=1` almost everywhere.
4. Periodicity of A gives `Σ_k|φ̂(ω−2πk)|²=1` almost everywhere.
5. To prove equality of spans, approximate the periodic bounded functions `A^(−1/2)` and `A^(1/2)` by trigonometric polynomials in L²(T). For the first direction the squared product error on R is the period integral of the squared multiplier error weighted by A, bounded by b times the ordinary error. For the reverse direction the periodized spectral energy of φ is exactly 1. Fourier inversion turns each polynomial multiplier into a finite linear combination of translates. Thus θ and φ belong to each other's closed translate spans.

This replaces the old criterion display and the two Poisson proof displays with the a.e. criterion and one aligned Fourier-coefficient calculation: three old displayed blocks become two. No label is removed or renamed; the now-unneeded Poisson and convolution reference uses disappear from this proof.

#### Other mathematical clarifications

- **Reused Haar panels:** `multires-1d-haar-1`, `-2`, and `-3` were labeled −8, −7, −6 in one figure and −7, −6, −5 in the approximation/detail figure. No generating source identifying the absolute calibration was found. With parent agreement, the latter now uses the consistent relative levels `j,j+1,j+2`; detail labels are `W_{j+1},W_{j+2}`. Its caption explains that increasing j **coarsens** the approximation. No figure asset, plotted data, or figure ID changed, and no absolute sample scale was guessed.
- **FWT algorithm:** explained why the real-filter derivation also applies to the algorithms' complex-valued inputs.
- **Anisotropic supports:** restricted rectangular compact-support claims to compactly supported wavelets and allowed the fixed support-size constant. Shannon and orthogonal spline wavelets need not have compact support.
- **2-D spaces:** clarified that the displayed span is closed in L² on unbounded domains and that the torus uses periodized atoms.
- **Low-pass condition C\*:** described nonvanishing on the central half-period as the sufficient nondegeneracy condition used here, rather than implying it is necessary for every multiresolution construction.
- **High-pass proof:** corrected the preceding-proof reference from C₁ (normalization) to C₂ (the even/odd spectral-energy argument). Declared that this proof establishes the forward implication; it does not silently claim to prove the converse.
- **Vanishing moments:** stated the induction hypothesis as vanishing of all lower-order derivatives rather than the ambiguous “assuming this holds for previous derivatives.”
- **Daubechies filters:** stated that the decimal filter values are rounded to four places. Thus the examples are not claimed to satisfy exact filter identities numerically. Identified the compactly supported orthogonal class over which the support length is minimal.

#### Writing and caption polish

Replaced “embedded spaces” by nested spaces and explained their detail components; removed a redundant assertion that the periodized family is a basis; explained the choice of coarsest level without implying it is forced; called forward filter-bank steps **decompositions** rather than refinements; removed “we use here,” “is re-written,” and “One indeed has that”; repaired singular/plural agreement in filter constraints and their caption; made the relation between design requirements and filter identities direct.

### Approximation — `chapters/approximation.tex`

#### Hölder/Sobolev endpoint correction

The chapter explicitly defines noninteger `C^α`, `α=m+β`, `0<β<1`, using β-Hölder continuity of the mth derivatives. Its Sobolev seminorm is the weighted ℓ² sum of Fourier coefficients on the periodic unit cube. Under these exact definitions, the prior unconditional same-order inclusion `C^α⊂H^α` is false for noninteger α.

The text now says:

- For integer α, the periodic C^α ball lies in an H^α ball after adjusting the radius. Bounded derivatives on the finite-volume domain give square-integrable derivatives; Parseval and the finite multinomial expansion of `|k|^(2α)` give the stated Fourier seminorm bound.
- For noninteger α, C^α regularity gives H^s regularity for each `s<α`, but not necessarily for `s=α`. A β-Hölder mth derivative has squared translation increments bounded by `C|h|^(2β)`; the fractional Sobolev integral of order γ<β is bounded near zero by `∫₀¹r^(2β−2γ−1)dr`.

An explicit periodic endpoint counterexample, confirming that the correction is needed, is
`f(x)=Σ_{k≥1}2^(−αk)cos(2π2^k x)` for noninteger `α=m+β`. Derivatives through order m converge uniformly. For the mth derivative, splitting its increment sum at `2^k|h|≈1` bounds
`Σ_k2^(−βk)min(2^k|h|,1)` by `C|h|^β`, proving C^α regularity. In the H^α Fourier seminorm each frequency ±2^k contributes a fixed positive quantity, so the sum diverges. This uses the precise periodic domain and norms already defined in the chapter, not a different Hölder convention.

The first ordinary-derivative display and the later introductory discussion also now explicitly assume integer α where needed. No Fourier decay or Sobolev error-bound formula was changed.

#### Other mathematical repairs

- **BV inclusion:** restricted inclusion of the piecewise-smooth models in bounded variation to `α≥1`, explaining the bounds on derivatives, jumps, and total edge length. Fractional Hölder regularity below one need not imply bounded variation.
- **Indicator total variation:** the equality between total variation and the length of `∂Ω` is now illustrated for a bounded set with piecewise smooth boundary. For an arbitrary finite-perimeter set, perimeter is the measure of its reduced boundary and need not equal the size of its topological boundary. The displayed formula is unchanged and valid under the stated regularity.
- **Indexing conventions:** stated that M is even in both Fourier cutoff constructions; retained and explained the M+1 count for the symmetric cutoff. Corrected the explanation of restriction notation to the actual `f|_(t_i,t_(i+1))` notation. Changed “fewer than K” cartoon edges to “at most K,” matching the union indexed 0 through K−1.
- **Fixed versus selected projection:** distinguished linear projection onto each fixed space from the possibly nonlinear overall map that selects the space from f.
- **Threshold ties:** limited the hard-thresholding rewriting to its threshold-selected coefficient set. Carried the tie convention into the curvelet paragraph.
- **Power laws:** distinguished a power-law upper bound from an exact power law, which alone becomes an exact straight line on logarithmic axes. Removed the redundant reference use that had equated the bound with that equality. The power-law display itself is unchanged.
- **One-dimensional coefficient counting:** restored the omitted **+1** in the second bound of `eq-wavapprox-1d-2`, retaining the coarsest scaling coefficient. Dropping a positive term in an upper-bound chain was unjustified. Its asymptotic estimate and final rate do not change.
- **Fourier/wavelet comparison:** stated the comparison for piecewise Lipschitz signals, the class actually covered by the referenced Fourier bound. The previous “α>1/2” also included nonsmooth fractional pieces outside that bound's stated hypotheses.
- **Curvelet threshold reconstruction:** used the complete index set Λ already introduced for the Parseval frame, rather than summing only the fine-scale `(j,m,k)` atoms. This explicitly retains the low-frequency family needed for a complete reconstruction. The rate claim is unchanged.
- **Support geometry:** restricted the square-support discussion to the isotropic compactly supported wavelets being compared with curvelets.

#### Editorial changes

Removed repeated “coefficients” in the opening approximation definition; repaired a run-on sentence across the Fourier/Sobolev displays; distinguished common decay exponents from signal-dependent constants; explained that the sharper Sobolev estimate uses a better argument under the same derivative assumptions; replaced “coefficient segmentation” and “counting the error” by descriptions of the actual proof steps; used “total length” for singular curves and “piecewise affine” for interpolation on triangles.

The only changed displayed blocks are the restored coarsest-coefficient count and the complete-index curvelet threshold formula. The 74 other displayed blocks are identical to the baseline.

### Compression — `chapters/compression.tex`

- Kept decoding applicable to both signals and images.
- Explained the threshold convention at the reference to nonlinear approximation: compression retains equality at T, whereas the earlier hard-thresholding formula used a strict inequality.
- Introduced the theorem by its rate bound for encoding the support and values. Kept the sufficient upper-approximation-bound qualification concise and linked compression **rates** to approximation **rates** explicitly.
- **Decoder information:** included the integer amplitude bound A in the header, alongside M and N. The value-code length uses `ceil(log₂(2A+1))`, so the decoder must know A. Under the proof's hypotheses A, M, and N grow at most polynomially in K, and their self-delimiting integer descriptions are logarithmic; the rate proof is unchanged.
- Described entropy coding as reducing **expected** file size, not guaranteeing that every individual sequence becomes shorter. Removed a redundant probability-model assumption and clarified the assignment of Huffman codeword lengths.
- Rewrote the three bit-plane steps with **the encoder** as subject; a coefficient does not itself encode bits. Identified when significance first occurs, when its sign is emitted, and when the next magnitude digit is emitted.
- Distinguished a coefficient from its position index in the context-model explanation.

All 24 display blocks, labels, reference uses, and JPEG-2000 assets are identical to the baseline.

### denoising

Clarified the distinction between a random observation and its realization, the centering and linearity assumptions, the meaning of filter normalization, frequency ordering, finite projection sizes, threshold ties, frame multiplicity, wavelet blocks, and variance-stabilizing transforms. Corrected the pixel-shift spacing and clean/noisy coefficient labels.

Coverage: the entire active chapter, including captions and proofs. Disabled blocks and substantive full-line comments match the baseline.

#### 1. Source line 57

- Before: `where $w$ is the noise residual. In the analysis, $w$ is random and so is the noisy image $f$; an acquired measurement is one realization.`
- After: `where $w$ is the noise residual. Both $w$ and $Y$ are random; an acquired image is one realization of $Y$, denoted by $y$ or $f$ below.`
- Rationale: Clarifies the distinction between the random observation Y and its realizations y or f; formulas unchanged.

#### 2. Source line 80

- Before: `The simplest model takes each $w_n$ to be a centered Gaussian with variance $\si^2$, with independent entries $w_n$: $w\sim\Nn(0,\sigma^2\Id_N)$. This is Gaussian white noise.`
- After: `The simplest model assumes independent centered Gaussian entries $w_n$ with variance $\si^2$: $w\sim\Nn(0,\sigma^2\Id_N)$. This is Gaussian white noise.`

#### 3. Source line 93

- Before: `In the additive model, $Y$ has mean $f_0$, so denoising estimates that mean from a single realization.`
- After: `With centered additive noise, $Y$ has mean $f_0$, so denoising estimates that mean from a single realization.`
- Rationale: Makes explicit the centering assumption required for the mean identity.

#### 4. Source line 125

- Before: `A linear estimator $\Dd(Y) = \tilde f$ of $f_0$ depends linearly on $Y$, so that $\Dd(f+g) = \Dd(f)+\Dd(g)$.`
- After: `A linear estimator $\Dd(Y) = \tilde f$ of $f_0$ preserves addition, $\Dd(f+g) = \Dd(f)+\Dd(g)$, and scalar multiplication.`
- Rationale: States both defining properties of linearity; additivity alone is insufficient without further assumptions.

#### 5. Source line 160

- Before: `where $Z_s$ ensures that $\sum_i h_{s,i}=1$ (low pass).`
- After: `where $Z_s$ normalizes the filter so that $\sum_i h_{s,i}=1$, preserving constant signals.`
- Rationale: Distinguishes unit DC gain from the separate low-pass property.

#### 6. Source line 192

- Before: `The denoising error can be decomposed as`
- After: `The triangle inequality separates the smoothing error from the remaining noise:`

#### 7. Source line 196

- Before: `The filter width $s$ should be chosen to balance noise removal ($\norm{h_s \star w}$ decreases with $s$) against excessive smoothing of singularities ($\norm{h_s\star f_0-f_0}$ increases with $s$).`
- After: `The filter width $s$ balances noise removal against excessive smoothing of singularities. A wider filter typically reduces $\norm{h_s \star w}$ while increasing $\norm{h_s\star f_0-f_0}$.`
- Rationale: Avoids asserting monotonicity for every finite-grid filter family and noise realization.

#### 8. Source line 277

- Before: `Since $W$ and $F$ are independent, the same computation as above carries over so that`
- After: `Independence of $W$ and $F$ allows the same risk calculation:`

#### 9. Source line 313

- Before: `This denoising scheme is thus parameterized by an integer $M>0$.`
- After: `This denoising scheme is parameterized by the number $1\leq M\leq N$ of retained coefficients.`
- Rationale: States the finite-dimensional range of M already used in the theorem.

#### 10. Source line 315

- Before: `For instance, when $\Bb$ is the discrete Fourier basis, this corresponds to an ideal low-pass filter against a (discretized) Dirichlet kernel.`
- After: `For the discrete Fourier basis ordered by increasing frequency, this is ideal low-pass filtering by convolution with a discrete Dirichlet kernel.`
- Rationale: The ordering assumption is necessary for the retained coefficients to be low frequencies.

#### 11. Source line 354

- Before: `\begin{prop} Assuming that \eql{\label{eq-discr-sobol}`
- After: `\begin{prop} Let $\alpha>0$. If \eql{\label{eq-discr-sobol}`
- Rationale: Makes the positive regularity exponent explicit in the weighted coefficient bound.

#### 12. Source line 359

- Before: `then \eq{ \foralls M, \quad \norm{f_0 - f_{0,M}^{\text{lin}}}^2 \leq C M^{-2\al}`
- After: `then, for $1\leq M\leq N$, \eq{ \foralls M, \quad \norm{f_0 - f_{0,M}^{\text{lin}}}^2 \leq C M^{-2\al}`
- Rationale: Specifies the admissible finite-dimensional range of the projection size.

#### 13. Source line 406

- Before: `Hard thresholding discards coefficients below a chosen magnitude:`
- After: `Hard thresholding retains only coefficients whose magnitude exceeds the threshold $T$:`
- Rationale: Matches the strict inequality in the displayed definition, including the threshold equality case.

#### 14. Source line 508

- Before: `Left: sparse signal $a$, right: noisy signal.`
- After: `Left: clean sparse coefficients $a_0$; right: noisy coefficients $a$.`
- Rationale: Corrects caption notation: a0 is clean, while a=a0+z is noisy.

#### 15. Source line 610

- Before: `$\De = \{ 0, 1/N, 2/N, 3/N \}^2$, which corresponds to discrete translation by $\{ 0, 1, 2, 3 \}^2$ on the discretized image.`
- After: `$\De = \{ 0, 1/N_0, 2/N_0, 3/N_0 \}^2$ on an $N_0\times N_0$ image, corresponding to shifts by $\{ 0, 1, 2, 3 \}^2$ pixels.`
- Rationale: Corrects the physical shift spacing to 1/N0, not 1/N: N=N0^2 is the total pixel count.

#### 16. Source line 643

- Before: `Although $\Bb_{\text{inv}}$ is no longer an orthogonal basis, it still satisfies a conservation of energy formula`
- After: `Although $\Bb_{\text{inv}}$ is no longer an orthogonal basis, orthonormality of each translated copy gives the energy and reconstruction identities`

#### 17. Source line 657

- Before: `The frame $\Bb_{\text{inv}}$ might contain up to $|\De| |\Bb|$ frame elements.`
- After: `Counted with multiplicity, the translated frame has $|\De| |\Bb|$ elements.`
- Rationale: Clarifies the distinction between the indexed frame and its distinct atoms; multiplicities affect tight-frame weights.

#### 18. Source line 795

- Before: `where each $B_k$ is a square of $s \times s$ coefficients; the block size $s$ is a parameter of the method.`
- After: `where each $B_k$ contains $s\times s$ neighboring coefficients at a fixed scale and orientation. The block size $s$ is a parameter of the method.`
- Rationale: Clarifies which neighboring coefficients belong to one wavelet block.

#### 19. Source line 802

- Before: `and the block thresholding \eq{`
- After: `The block estimator \eq{`

#### 20. Source line 866

- Before: `Many imaging devices sample an image through a photon-counting operation.`
- After: `Many imaging devices measure intensity by counting photons.`

#### 21. Source line 917

- Before: `A variance-stabilizing transform $\phi : \RR \rightarrow \RR$ first maps the image to $\phi(f)$, for which an additive Gaussian white-noise model is approximately appropriate:`
- After: `A variance-stabilizing transform $\phi : [0,+\infty) \rightarrow \RR$, applied pixelwise, maps the counts to $\phi(f)$. An additive Gaussian white-noise model can then be a useful approximation:`
- Rationale: Corrects the domain for the square-root transforms that follow; neither is defined on all real inputs.

#### 22. Source line 924

- Before: `Two popular variance stabilization functions for Poisson noise are the Anscombe mapping`
- After: `Two common variance-stabilizing transforms for Poisson noise are the Anscombe transform`

#### 23. Source line 928

- Before: `and the mapping of Freeman and Tukey`
- After: `and the Freeman--Tukey transform`

#### 24. Source line 1018

- Before: `where $\Ga$ denotes the Gamma function, which extends the factorial to noninteger arguments. One thus has`
- After: `Here $\Ga$ is the Gamma function and $\psi$ is its logarithmic derivative, the digamma function. The transformed observation satisfies`

#### 25. Source line 1022

- Before: `for strictly positive intensities, where $z_n=\log(w_n)-c$ is a zero-mean additive noise with variance $\psi^{(1)}(K)$.`
- After: `for strictly positive intensities, where $z_n=\log(w_n)-c$ is centered additive noise with variance $\psi^{(1)}(K)$; $\psi^{(1)}$ is the derivative of the digamma function.`
- Rationale: Explains the previously unexplained trigamma notation; the variance formula is unchanged.

#### 26. Source line 359

- Before: `then, for $1\leq M\leq N$, \eq{ \foralls M, \quad`
- After: `then \eq{ \foralls M\in\{1,\ldots,N\}, \quad`
- Rationale: Put the admissible finite-dimensional range directly in the displayed quantifier, avoiding a duplicate quantifier; the approximation bound is unchanged.

#### 27. Source line 800

- Before: `E_k=\frac{1}{|B_k|} \sum_{ m \in B_k } |\dotp{f}{\psi_m}|^2, } The block estimator`
- After: `E_k=\frac{1}{|B_k|} \sum_{ m \in B_k } |\dotp{f}{\psi_m}|^2. } The block estimator`
- Rationale: End the sentence defining block energy before introducing the estimator.

#### 28. Source line 808

- Before: `S_T^{\text{block},q} (\dotp{f}{\psi_m}) = A_T^q(\sqrt{E_k}) \dotp{f}{\psi_m}. } for $q`
- After: `S_T^{\text{block},q} (\dotp{f}{\psi_m}) = A_T^q(\sqrt{E_k}) \dotp{f}{\psi_m}, } for $q`
- Rationale: Keep the parameter range in the sentence containing the block-threshold definition.

### variational-priors

Clarified the discrete-gradient scaling and the isotropic TV norm. Qualified the general gradient-descent discussion, explained continuous-time scaling, distinguished spatial gradients from energy gradients, corrected the sign in the TV illustration description, and improved the transitions between flows, fidelity regularization, smoothing, and nonsmooth algorithms.

Coverage: the entire active chapter, including captions and proofs. Disabled blocks and substantive full-line comments match the baseline.

#### 1. Source line 23

- Before: `When finite, the prior energy $J(f) \in \RR$ should be small for images in the target class $f \in \Theta$.`
- After: `A prior energy assigns small values to images in the target class $\Theta$ and larger values, possibly $+\infty$, to images that depart from that model.`

#### 2. Source line 60

- Before: `The discrete gradient at each pixel is`
- After: `With the grid-spacing factor absorbed into the energy scaling, the discrete gradient at each pixel is`
- Rationale: Explains why finite differences omit division by the pixel spacing.

#### 3. Source line 68

- Before: `Figure \ref{fig-grad-discrete} shows gradient vectors pointing in the direction of steepest increase of the function represented by $f$.`
- After: `Figure \ref{fig-grad-discrete} shows discrete gradient vectors approximating the local direction of steepest intensity increase.`
- Rationale: Avoids equating finite differences with exact continuous derivatives.

#### 4. Source line 169

- Before: `Similarly, a discrete TV energy is defined as the $\lun$ norm of the gradient field`
- After: `The isotropic discrete TV energy instead sums the Euclidean magnitudes of the gradient vectors:`
- Rationale: Clarifies that this is a mixed l1-l2 norm, not entrywise anisotropic l1.

#### 5. Source line 174

- Before: `where the $\lun$ norm of a vector field $v \in \RR^{N \times 2}$ is`
- After: `Here the notation $\normu{v}$ for a vector field $v \in \RR^{N \times 2}$ means`

#### 6. Source line 195

- Before: `If $J$ is a smooth function of the image $f$, the discrete energy can be minimized by gradient descent`
- After: `Gradient descent seeks to decrease a smooth discrete energy by iterating`
- Rationale: Avoids an unconditional claim of global minimization for a general smooth, possibly nonconvex prior.

#### 7. Source line 199

- Before: `where the step size $\tau$ must be small enough to guarantee convergence.`
- After: `where the step size $\tau$ is chosen to ensure descent. The quadratic and smoothed TV energies below admit explicit sufficient bounds on $\tau$.`
- Rationale: Replaces an unsupported general convergence assertion with the specific convergence bounds established below.

#### 8. Source line 201

- Before: `As the step size $\tau$ tends to zero, the iteration index $k$ is replaced by continuous time.`
- After: `To pass to continuous time, associate iteration $k$ with time $t=k\tau$ and let $\tau$ tend to zero.`

#### 9. Source line 257

- Before: `Discrete Laplacian and discrete TV gradient.`
- After: `Discrete Laplacian (left) and negative discrete TV gradient (right).`
- Rationale: Corrects the sign in the caption: the displayed divergence is minus the TV energy gradient.

#### 10. Source line 303

- Before: `This formula is singular at points where the gradient vanishes. Figure \ref{fig-discrete-operators-lapl} shows a TV gradient, which appears noisy in smooth areas because $\norm{\nabla f_n}$ is small in such regions.`
- After: `The displayed quotient is undefined where the spatial gradient vanishes. Figure \ref{fig-discrete-operators-lapl}, right, shows the negative TV gradient; normalization by small gradient magnitudes makes it sensitive to perturbations in smooth regions.`
- Rationale: Matches the figure sign and identifies which gradient occurs in the denominator.

#### 11. Source line 335

- Before: `Figure \ref{fig-discrete-operators-tveps} displays the regularized gradient magnitude for several smoothing parameters.`
- After: `Figure \ref{fig-discrete-operators-tveps} displays the regularized spatial-gradient magnitude for several smoothing parameters.`

#### 12. Source line 359

- Before: `In practice, the flow is computed with a discrete gradient descent \eqref{eq-gradmin-flow}.`
- After: `In practice, gradient descent~\eqref{eq-gradmin-flow} discretizes this flow in time.`

#### 13. Source line 361

- Before: `Figure \ref{fig-flow} shows a comparison between the heat flow and the total variation flow for a small value of $\epsilon$.`
- After: `Figure \ref{fig-flow} compares heat flow with smoothed TV flow for a small $\epsilon$.`

#### 14. Source line 399

- Before: `Figure \ref{fig-denoising-cameraman}, top row, shows an example of this oracle estimation of the best stopping time.`
- After: `Figure \ref{fig-denoising-cameraman}, top row, shows reconstructions obtained with such an oracle choice of stopping time.`

#### 15. Source line 461

- Before: `where the descent step size $\tau>0$ should be small enough. This gradient descent is similar to the time-discretized minimization flow \eqref{eq-gradmin-flow}, except that the data fidelity term generally prevents the flow from converging to a constant image.`
- After: `where $\tau>0$ is chosen according to the smoothness of the objective. The iteration adds a data fidelity force to the time-discretized flow~\eqref{eq-gradmin-flow}, counteracting the prior's tendency to smooth the image toward a constant.`

#### 16. Source line 463

- Before: `The following sections explain how to compute $f^\star_\la$ for these priors.`
- After: `Computing $f^\star_\la$ then requires smoothing the prior or using algorithms adapted to nonsmooth objectives, as developed later in the book.`

#### 17. Source line 505

- Before: `Section~\ref{sec-fb-dual} presents an alternative that avoids this $\epsilon$ smoothing.`
- After: `Section~\ref{sec-fb-dual} presents an algorithm for the original nonsmooth TV penalty.`

### inverse-problems

Clarified the operator and noise assumptions, what an oracle knows, the zero-noise limit, the capped-inverse/Tikhonov comparison, and the identifiability needed to recover kernel components. Improved the compactness argument, Sobolev-source interpretation, BV definition, and tomography geometry.

Coverage: the entire active chapter, including captions and proofs. Disabled blocks and substantive full-line comments match the baseline.

#### 1. Source line 13

- Before: `Increasing the resolution of signals and images requires solving an ill-posed inverse problem: inverting a linear measurement operator that reduces resolution. This chapter uses the convex regularization introduced in Chapter \ref{chap-variational} to stabilize this inverse problem.`
- After: `Recovering fine detail from coarse measurements leads to an inverse problem: invert a measurement operator that suppresses or discards information. Such problems are often ill-posed. This chapter uses the convex regularization introduced in Chapter \ref{chap-variational} to stabilize recovery.`
- Rationale: Avoids suggesting that every resolution problem is necessarily ill-posed or linear.

#### 2. Source line 15

- Before: `We consider a linear map $\Phi : \Ss \rightarrow \Hh$, usually assumed continuous, where $\Ss$ is a Hilbert space or, more generally, a Banach space.`
- After: `We consider a bounded linear map $\Phi : \Ss \rightarrow \Hh$, where the signal space $\Ss$ is a Hilbert space or, more generally, a Banach space. The data space $\Hh$ is a Hilbert space.`
- Rationale: Makes the boundedness and data-space structure used in the subsequent squared-norm analysis explicit.

#### 3. Source line 23

- Before: `In most applications, $\Hh=\RR^P$ is finite-dimensional because an acquisition device records only finitely many observations, whose number $P$ is often small.`
- After: `An acquisition device records finitely many observations, so applications usually take $\Hh=\RR^P$. The number $P$ may be small relative to the desired reconstruction dimension.`

#### 4. Source line 47

- Before: `For example, convolution can be defined on $\Ss=\Hh=L^2(\TT^d)$; see Proposition~\ref{prop-ti-convol-l2}.`
- After: `A periodic model takes $\Ss=\Hh=L^2(\TT^d)$; see Proposition~\ref{prop-ti-convol-l2}.`

#### 5. Source line 118

- Before: `We first seek an optimal estimator among inversion methods that act diagonally in a fixed basis.`
- After: `We first seek an oracle estimator among inversion methods that act diagonally in a fixed basis, allowing the gains to depend on the unknown clean signal.`
- Rationale: Explains the information available to the oracle before optimization uses the clean coefficients.

#### 6. Source line 130

- Before: `where the noise coefficients are centered Gaussian variables of variance $\sigma^2$. We work in finite dimension; in infinite dimension, white noise must be interpreted as a generalized random process.`
- After: `where $w$ is Gaussian white noise with variance $\sigma^2$. In a real orthonormal basis, its coefficients are independent centered Gaussians. In a complex Fourier basis they still have mean zero and second moment $\sigma^2$, which suffices for the risk calculation below. We work in finite dimension; infinite-dimensional white noise requires a generalized random-process interpretation.`
- Rationale: Clarifies the second-moment requirement in the complex Fourier case without incorrectly asserting independence of conjugate-frequency coefficients.

#### 7. Source line 157

- Before: `As $\sigma \to 0$, assuming $\phi_k \neq 0$ and $c_k\neq0$, the oracle estimator approaches direct inversion, independently of the signal amplitude: $$ \lambda_k = \frac{1}{\phi_k}. $$`
- After: `For $\phi_k \neq 0$ and $c_k\neq0$, the oracle gains approach direct inversion as the noise vanishes: $$ \lambda_k \longrightarrow \frac{1}{\phi_k}\qquad\text{as }\sigma\to0. $$`
- Rationale: Corrects an equality for positive noise to the limit actually justified by the preceding formula.

#### 8. Source line 228

- Before: `Computing the SVD of a full matrix $\Phi \in \RR^{N \times N}$ has complexity $N^3$.`
- After: `Computing the SVD of a dense matrix $\Phi \in \RR^{N \times N}$ typically costs $O(N^3)$ operations.`

#### 9. Source line 245

- Before: `For an operator $\Phi$ with the expansion~\eqref{eq-svd-operators}, $\norm{\Phi}_{\Ll(\Ss,\Hh)}=\si_1$.`
- After: `For a nonzero operator $\Phi$ with the expansion~\eqref{eq-svd-operators}, $\norm{\Phi}_{\Ll(\Ss,\Hh)}=\si_1$.`
- Rationale: The zero operator has no positive singular value sigma1 in the reduced decomposition.

#### 10. Source line 338

- Before: `Figure~\ref{fig-bound-regul-variance}, left, illustrates a regularized inverse obtained by thresholding singular values.`
- After: `Figure~\ref{fig-bound-regul-variance}, left, illustrates a capped reciprocal filter, which limits the amplification of small singular values.`
- Rationale: Uses the same capped-inverse interpretation as the existing original/TikZ comparison notes.

#### 11. Source line 385

- Before: ``The aim is to choose $\lambda$ as a function of the noise level and to quantify the rate at which $f_\lambda$ approaches $f_0$. One first needs to ensure at least $f_0=f_0^+$, which in turn requires that $f_0\in\overline{\Im(\Phi^*)}=\ker(\Phi)^\bot$. Indeed, a limitation of the squared-norm regularization studied here is that necessarily $f_\la\in\Im(\Phi^*)\subset\ker(\Phi)^\bot$ so that no information can be recovered inside $\ker(\Phi)$. Nonlinear priors can achieve a ``super-resolution'' effect and recover this missing information.``
- After: ``We choose $\lambda$ according to the noise level and seek a rate for $f_\lambda\to f_0$. This requires $f_0=f_0^+$, equivalently $f_0\in\overline{\Im(\Phi^*)}=\ker(\Phi)^\bot$: squared-norm regularization always produces $f_\la\in\Im(\Phi^*)\subset\ker(\Phi)^\bot$, so it cannot recover a kernel component of $f_0$. Nonlinear priors can recover such components when additional structure makes the signal identifiable.``
- Rationale: Qualifies recovery of unobserved components by the necessary identifiability assumption.

#### 12. Source line 401

- Before: `Since the vectors $v_m$ are Fourier modes, \eqref{eq-source-cond-quad} becomes`
- After: `Since the vectors $v_m$ are Fourier modes, \eqref{eq-source-cond-quad} gives, up to constants and the Fourier normalization, a bound of the form`
- Rationale: Polynomial asymptotics of the multiplier imply an equivalent weighted bound, not exact equality of the ball radius.

#### 13. Source line 405

- Before: `which defines a Sobolev ball with radius $\rho$ and smoothness order $\al\be$.`
- After: `This is a Sobolev-type coefficient bound of smoothness order $\al\be$.`

#### 14. Source line 416

- Before: `\begin{thm}\label{thm-sublin-quad} Under the source condition`
- After: `\begin{thm}\label{thm-sublin-quad} Let $\delta>0$. Under the source condition`
- Rationale: The displayed positive regularization-parameter choice uses a positive noise bound; delta=0 is a separate limiting case.

#### 15. Source line 491

- Before: `Bounding $\mu_{\la}(\si) \leq C_\la = \frac{1}{2\sqrt{\la}}$.`
- After: `Left: capping reciprocal singular values. Right: the Tikhonov bound $\mu_{\la}(\si) \leq C_\la = \frac{1}{2\sqrt{\la}}$.`
- Rationale: The bound 1/(2 sqrt(lambda)) belongs to the right Tikhonov panel, not to the left capped reciprocal curve.

#### 16. Source line 507

- Before: `The rate~\eqref{eq-rate-tikhon} saturates:`
- After: `The rate in Theorem~\ref{thm-sublin-quad} saturates:`
- Rationale: References the rate theorem rather than the triangle-inequality decomposition.

#### 17. Source line 513

- Before: `Figure~\ref{fig-bound-regul-variance} compares the curve for quadratic regularization~\eqref{eq-sol-quad-reg-svd} (right) with the simpler thresholding curve (left) which does not suffer from saturation.`
- After: `The capped inverse in Figure~\ref{fig-bound-regul-variance}, left, avoids this finite-order saturation, whereas Tikhonov regularization is limited to source orders $\beta\leq2$.`
- Rationale: Names the illustrated filter precisely and explains the comparison.

#### 18. Source line 533

- Before: `The Lagrange multiplier $\la$ balances these two terms and can be difficult to choose in practice.`
- After: `The regularization parameter $\la>0$ balances these two terms and can be difficult to choose in practice.`
- Rationale: Avoids identifying a penalty parameter with a Lagrange multiplier before a constrained formulation is specified.

#### 19. Source line 539

- Before: `We may assume $y \in \Im(\Phi)$ without loss of generality because the orthogonal decomposition`
- After: `We may assume $y \in \Im(\Phi)$ without changing the minimizers, since orthogonality gives`

#### 20. Source line 543

- Before: `so that one can replace $y$ by $\Proj_{\Im(\Phi)}(y)$ in~\eqref{eq-ip-regul}.`
- After: `The first term is independent of $f$, so replacing $y$ by $\Proj_{\Im(\Phi)}(y)$ in~\eqref{eq-ip-regul} changes only an additive constant.`

#### 21. Source line 553

- Before: `Equivalently, every nonempty sublevel set $\enscond{f}{J(f) \leq c}$ is bounded, and is compact if $J$ is lower semicontinuous, for all $c$.`
- After: `Equivalently, every sublevel set $\enscond{f}{J(f) \leq c}$ is bounded. In finite dimension, lower semicontinuity of $J$ then makes these sets compact.`

#### 22. Source line 563

- Before: `Let $h$ be any solution of~\eqref{eq-ip-regul-noiseless}, so that $\Phi h = y$.`
- After: `Lower semicontinuity and coercivity ensure that $J$ attains its minimum on the nonempty closed constraint set $\{f:\Phi f=y\}$. Let $h$ be such a minimizer.`
- Rationale: Justifies existence of the constrained minimizer invoked by the compactness proof.

#### 23. Source line 694

- Before: `More generally, a function $f$ is of bounded variation when its distributional gradient $Df$ is a finite vector-valued Radon measure.`
- After: `More generally, an integrable function $f$ is of bounded variation when its distributional gradient $Df$ is a finite vector-valued Radon measure.`
- Rationale: Includes the L1 requirement in the standard BV(R^d) definition.

#### 24. Source line 733

- Before: `The two regimes $\epsilon \rightarrow \{0,+\infty\}$ have the asymptotic descriptions`
- After: `For fixed $u$, the limits as $\epsilon\downarrow0$ and $\epsilon\to+\infty$ are described by`
- Rationale: States that the large-epsilon remainder is for fixed u, not uniform over all vectors.

#### 25. Source line 931

- Before: `In medical imaging, a scanner computes projections of the human body along rays $\De_{t,\th}$ defined by`
- After: `In an idealized two-dimensional tomography model, each projection integrates the unknown image along a line $\De_{t,\th}$ defined by`
- Rationale: The set defined by one affine equation in R2 is a line, not a ray; clarifies the idealized acquisition model.

#### 26. Source line 935

- Before: `where we restrict ourselves to two-dimensional projections to simplify the exposition.`
- After: `Here $\tau_\theta=(\cos\theta,\sin\theta)$ is the unit normal to the line.`
- Rationale: Defines the normal vector used throughout the Radon and backprojection formulas.

#### 27. Source line 937

- Before: `The scanning process computes a Radon transform, which integrates the function to be reconstructed along rays`
- After: `The resulting Radon transform is`

#### 28. Source line 982

- Before: `Relation \eqref{eq-fourier-slice} shows that knowing $R f$ is equivalent to knowing the Fourier transform of $f$ along rays,`
- After: `By~\eqref{eq-fourier-slice}, knowing $R f$ is equivalent to knowing the Fourier transform of $f$ along the corresponding radial lines:`
- Rationale: The frequency parameter xi ranges over all real values, giving full lines.

#### 29. Source line 403

- Before: `\sum_m \norm{m}^{2\al\be} |\hat f_m|^2 \leq \rho^2 < +\infty, } This is a Sobolev-type`
- After: `\sum_m \norm{m}^{2\al\be} |\hat f_m|^2 \leq \rho^2 < +\infty. } This is a Sobolev-type`
- Rationale: Match the displayed formula punctuation to the following complete sentence; inequality unchanged.

#### 30. Source line 383

- Before: `Figure~\ref{fig-bound-regul-variance} compares the curve for quadratic regularization~\eqref{eq-sol-quad-reg-svd} (right) with the simpler thresholding curve (left).`
- After: `Figure~\ref{fig-bound-regul-variance} compares quadratic regularization~\eqref{eq-sol-quad-reg-svd} (right) with a capped reciprocal filter (left).`
- Rationale: Use the same precise name for the left-hand filter throughout its discussion; it caps reciprocal singular values rather than discarding them.

### sparse-regularization

Improved the sparsity and convex-relaxation discussion, the analysis/synthesis terminology, and the algorithmic transitions. Clarified threshold tie conventions, repaired the scalar-objective caption, stated consistency for the equality-constrained limit, and matched the forward–backward split to the coefficient-space objective.

Coverage: the entire active chapter, including captions and proofs. Disabled blocks and substantive full-line comments match the baseline.

#### 1. Source line 17

- Before: `Chapter \ref{chap-approximation} shows how an orthonormal basis $\Bb = \{ \psi_m \}_m$ adapted to an image $f$ in a class $f \in \Th$ can approximate it using only a few atoms from $\Bb$.`
- After: `Chapter \ref{chap-approximation} shows that an orthonormal basis $\Bb = \{ \psi_m \}_m$ adapted to an image class $\Th$ can approximate each $f\in\Th$ using relatively few atoms.`

#### 2. Source line 19

- Before: `We measure representation complexity in $\Bb$ by the $\lzero$ prior, which counts nonzero coefficients in $\Bb$:`
- After: `The $\lzero$ prior measures representation complexity by counting the nonzero coefficients:`

#### 3. Source line 29

- Before: `It provides an ideal sparsity measure for the coefficient vector $x$ of $f$ in $\Bb$.`
- After: `It measures sparsity directly, without accounting for the magnitudes of the nonzero coefficients.`

#### 4. Source line 37

- Before: `For suitable wavelet bases, bounded-variation image models give decay estimates for the approximation error $\norm{f-f_M}$; see Section \ref{sec-signal-models}.`
- After: `For suitable wavelet bases and the bounded-variation image classes specified in Section \ref{sec-signal-models}, the error $\norm{f-f_M}$ satisfies quantitative decay estimates.`
- Rationale: Ties the approximation claim to the specific bounded BV assumptions already stated in the referenced section.

#### 5. Source line 55

- Before: `The ideal sparsity prior $\Jz$ is difficult to optimize because $\Jz(f)$ is nonconvex as a function of $f$.`
- After: `The ideal sparsity prior $\Jz$ is nonconvex, which makes optimization difficult.`

#### 6. Source line 57

- Before: `Minimizing $\Jz$ in a general inverse problem leads to a difficult combinatorial optimization problem; see Chapter~\ref{chap-invpbm}.`
- After: `In a general inverse problem, minimizing $\Jz$ entails a combinatorial search over possible coefficient supports.`

#### 7. Source line 73

- Before: `The boundary of this convex range gives the $\lun$ prior $\Ju$:`
- After: `The smallest convex exponent, $q=1$, gives the $\lun$ prior $\Ju$:`

#### 8. Source line 87

- Before: `It can be rewritten in the orthonormal basis as`
- After: `Orthonormality separates the objective into coefficientwise terms:`

#### 9. Source line 128

- Before: `is the hard thresholding introduced in~\eqref{eq-hard-thresh-denoise}, and`
- After: `is hard thresholding, with the alternative tie convention discussed below, and`
- Rationale: Avoids calling two threshold rules identical when one keeps and the other discards coefficients exactly at T.

#### 10. Source line 153

- Before: `Remaining panels: dependence on $\la$ of $F(x) \eqdef \frac{1}{2}\norm{\cdot-y}^2+\la |\cdot|$.`
- After: `Remaining panels: the scalar objective $F(x) \eqdef \frac{1}{2}|x-y|^2+\la |x|$ for increasing $\la$.`
- Rationale: Repairs the argument mismatch in the caption definition of F(x); the scalar objective and normalization are unchanged.

#### 11. Source line 196

- Before: `The analysis prior measures the sparsity of correlations with the dictionary atoms`
- After: `The analysis prior penalizes the sum of the magnitudes of correlations with the dictionary atoms:`

#### 12. Source line 214

- Before: `As $\la \rightarrow 0$, we consider the limiting constrained problem`
- After: `For consistent data $y\in\Im(A)$, the limit $\la\downarrow0$ leads to the constrained problem`
- Rationale: The stated equality-constrained limit is feasible only for data in the range of A.

#### 13. Source line 339

- Before: `The upper bound and tangency imply`
- After: `The surrogate upper bound and equality at the current iterate imply`

#### 14. Source line 369

- Before: `For~\eqref{eq-bpdn}, use the splitting`
- After: `For the coefficient-space objective~\eqref{eq-lasso-lagr-ip}, rescaled by $\lambda$, use the splitting`
- Rationale: Matches the displayed split involving A and ||x||1 to its coefficient-space objective, rather than the signal-space objective involving Phi.

#### 15. Source line 393

- Before: `The limited frequency band illustrates the information lost from $f$ during acquisition.`
- After: `The narrow transmitted band illustrates how acquisition suppresses information about $f_0$.`

#### 16. Source line 480

- Before: `Sparse 1D deconvolution using orthogonal wavelets.`
- After: `One-dimensional deconvolution with sparsity in an orthogonal wavelet basis.`

#### 17. Source line 532

- Before: `For noiseless inpainting, a small $\lambda>0$ approximates the constrained minimum-sparsity problem; exact interpolation is recovered only in the limit, under the assumptions of Proposition~\ref{prop-gamma-conv-regul}.`
- After: `For noiseless inpainting, a small $\lambda>0$ approximates constrained $\ell^1$ minimization. Exact interpolation is obtained in the limit under the assumptions of Proposition~\ref{prop-gamma-conv-regul}.`
- Rationale: Distinguishes the l1 relaxation from exact minimization of the number of nonzero coefficients.

#### 18. Source line 73

- Before: `The smallest convex exponent, $q=1$, gives the $\lun$ prior $\Ju$:`
- After: `The smallest exponent yielding a convex prior, $q=1$, gives the $\lun$ prior $\Ju$:`
- Rationale: Attribute convexity to the prior rather than the numerical exponent.

### sparse-theory

Clarified the polytope and sign-cone geometry, the distinction between existence of one independent-support minimizer and properties of all minimizers, and the uniqueness of fitted data. Refined the Bregman, entropy, certificate, support-stability, ERC, and off-grid explanations, with explicit assumptions and corrected matrix-norm interpretations.

Coverage: the entire active chapter, including captions and proofs. Disabled blocks and substantive full-line comments match the baseline.

#### 1. Source line 21

- Before: `We study the penalized problem~\eqref{eq-lasso-lagr-ip} and its constrained counterpart~\eqref{eq-lasso-constr-ip}, written here as`
- After: `We study the penalized problem~\eqref{eq-lasso-lagr-ip}, with $\lambda>0$, and its constrained counterpart~\eqref{eq-lasso-constr-ip}:`
- Rationale: States the positive parameter required by the division by lambda.

#### 2. Source line 63

- Before: `Figure~\ref{fig-homotopy} illustrates the exact or approximate recovery of sparse vectors $x_0$ and the importance of choosing the parameter $\la$ carefully.`
- After: `Figure~\ref{fig-homotopy} shows how the recovered coefficients change with $\lambda$, illustrating the effect of regularization on sparsity and reconstruction accuracy.`

#### 3. Source line 100

- Before: `The vectors recovered by $\ell^1$ minimization are mapped by $A$ to the boundary of the scaled polytope $\norm{x_0}_1 A B$.`
- After: `For a nonzero recovered vector $x_0$, its measurement $Ax_0$ lies on the boundary of the scaled polytope $\norm{x_0}_1 A B$.`
- Rationale: Keeps the nonzero assumption required by the polytope characterization explicit.

#### 4. Source line 106

- Before: `The cones $A\Cc_s$ partition $\RR^P$.`
- After: `Including the lower-dimensional cones and the origin, their images $A\Cc_s$ partition $\RR^P$.`
- Rationale: Clarifies boundary strata in the sign-pattern partition under the stated uniqueness assumption.

#### 5. Source line 154

- Before: `In particular, $\supp(x_\la) \subset \sat(\eta_\la)$.`
- After: `In particular, $\supp(x_\la) \subset \sat(\eta_\la)$, where $\sat(\eta)=\{i:|\eta_i|=1\}$ is the saturation set.`
- Rationale: Defines the saturation notation at its first use.

#### 6. Source line 186

- Before: `The next proposition guarantees a Lasso solution supported on linearly independent regressors, that is, columns of $A$. This property motivates the term basis pursuit.`
- After: `The next proposition shows that at least one minimizer uses linearly independent columns of $A$. It does not assert that every minimizer has this property.`

#### 7. Source line 196

- Before: `\wrapf{sparse-theory/proofs/injectivity}{Trajectory $(x_t)_t$.}`
- After: `\wrapf{sparse-theory/proofs/injectivity}{Coefficient trajectories along a support-reducing direction.}`

#### 8. Source line 206

- Before: `Although a coefficient vector $x_\la$ solving~\eqref{eq-lasso-lagr} may not be unique, its fitted data, or predictor, are unique: the denoised version of $y$ is $Ax_\la$.`
- After: `Although the coefficient vector $x_\la$ may not be unique, its fitted data $Ax_\la$, also called the predictor, are unique.`

#### 9. Source line 320

- Before: `For $\la>0$, the objective maximized in~\eqref{eq-dual-lasso} is strongly concave.`
- After: `For $\la>0$, the dual objective~\eqref{eq-dual-lasso} is strongly concave.`

#### 10. Source line 339

- Before: `for a convex regularizer $J$.`
- After: `for a proper convex regularizer $J$ and $\lambda>0$. We state estimates for any minimizer, without assuming uniqueness.`
- Rationale: Makes the domain and parameter assumptions explicit; the later theorem is conditional on minimizer existence.

#### 11. Source line 351

- Before: `For any $\eta \in \partial J(x_0)$, we define the associated Bregman divergence as`
- After: `For $x_0$ with nonempty subdifferential and a chosen $\eta \in \partial J(x_0)$, define the Bregman divergence`

#### 12. Source line 365

- Before: `If $J(x)=\sum_i x_i(\log x_i-1)+\iota_{\RR_+^N}(x)$ is the entropy, then`
- After: `For the negative-entropy functional $J(x)=\sum_i x_i(\log x_i-1)+\iota_{\RR_+^N}(x)$, one obtains`
- Rationale: Corrects the sign convention: the displayed convex functional is negative entropy, up to an affine term.

#### 13. Source line 369

- Before: `is the so-called Kullback--Leibler divergence for $x_0\in\RR_{++}^N$ and $x\in\RR_+^N$, with the convention $0\log0=0$.`
- After: `This is the generalized Kullback--Leibler divergence for $x_0\in\RR_{++}^N$ and $x\in\RR_+^N$, with the convention $0\log0=0$. For probability vectors, the linear terms sum to zero.`
- Rationale: Distinguishes the generalized divergence for arbitrary nonnegative vectors from the probability-vector case.

#### 14. Source line 400

- Before: `Since $D_\eta(x_\lambda|x_0)\geq0$, the completed-square identity gives`
- After: `Since $D_\eta(x_\lambda|x_0)\geq0$, the completed-square bound gives`
- Rationale: The displayed derivation is an inequality obtained from optimality, not an exact identity for the divergence.

#### 15. Source line 435

- Before: `The next lemma shows that the $\ell^1$ Bregman divergence controls the $\ell^1$ error on coordinates where $\eta$ is strictly below the saturation threshold.`
- After: `The next lemma quantifies this control through the gap between $|\eta_i|$ and one.`

#### 16. Source line 497

- Before: `where we used $x_{0,J^c}=0$ and~\eqref{eq-predict-source}.`
- After: `where we used $x_{0,J^c}=0$ and~\eqref{eq-predict-source}. The mixed operator norm $\norm{B}_{p,q}$ maps $\ell^p$ to $\ell^q$, as defined in Remark~\ref{rem-operator-norm}.`
- Rationale: Explains the mixed norms at their first substantive use; their definition otherwise appears later.

#### 17. Source line 543

- Before: `Applying this theorem requires constructing a valid certificate~\eqref{eq-sourcecond-l1}. % We now give a construction that, when valid, guarantees both a linear rate and support stability.`
- After: `We now construct a candidate certificate explicitly. When it satisfies a strict off-support bound, it guarantees both a linear rate and support stability.`

#### 18. Source line 568

- Before: `The main obstacle in computing~\eqref{eq-minnorm-certif} is the inequality constraint $\norm{\eta_0}_\infty\leq1$.`
- After: `The inequality constraint $\norm{\eta_0}_\infty\leq1$ makes the minimum-norm problem~\eqref{eq-minnorm-certif} difficult to solve explicitly.`

#### 19. Source line 634

- Before: `Choosing $\lambda$ proportional to $\norm{A^*w}_\infty$ with a sufficiently large constant gives an $O(\norm{w})$ error as the noise tends to zero.`
- After: `When $\norm{A^*w}_\infty>0$, choosing $\lambda$ proportional to this quantity with a sufficiently large constant gives an $O(\norm{w})$ error as the noise tends to zero.`
- Rationale: Excludes the zero-correlation case from a rule that would otherwise set lambda=0 despite the theorem requiring lambda>0.

#### 20. Source line 695

- Before: `Theorem~\ref{thm-support-stable} requires small $\norm{w}$ and $\la$ so that regularization does not eliminate small nonzero coefficients of $x_0$.`
- After: `Theorem~\ref{thm-support-stable} requires small correlated noise $\norm{A^*w}_\infty$ and a small $\lambda$, with their ratio controlled, so that regularization does not eliminate the nonzero coefficients of $x_0$.`
- Rationale: States the actual correlated-noise assumption; components of w orthogonal to the range of A do not affect the estimator.

#### 21. Source line 700

- Before: `\item $L$ accounts for the largest correlation between atoms inside and outside the support;`
- After: `\item $L$ is the largest sum of absolute correlations between an off-support atom and all atoms on the support;`
- Rationale: Corrects the interpretation of L=||A_Ic^* A_I||_infinity: it is a maximum row sum, not a largest individual correlation.

#### 22. Source line 703

- Before: `Choosing $\lambda=(1+\gamma)R\delta/S$ for any fixed $\gamma>0$ and sufficiently small $\delta$ gives`
- After: `For fixed $\gamma>0$ and sufficiently small $\delta>0$, choosing $\lambda=(1+\gamma)R\delta/S$ gives`
- Rationale: Keeps lambda positive in the explicit parameter choice.

#### 23. Source line 728

- Before: `Fuchs precertificate $\eta_F$ for a Gaussian matrix $A \in \RR^{p \times n}$ with $n=64$.`
- After: `Fuchs precertificate $\eta_F$ for a Gaussian matrix $A \in \RR^{P \times N}$ with $N=64$ and increasing measurement count $P$.`
- Rationale: Matches the caption dimensions to the P and N used in the chapter and panel labels.

#### 24. Source line 740

- Before: `Assume $A_I$ is injective and $\operatorname{ERC}(I)>0$.`
- After: `The maximum over an empty complement is taken to be zero. Assume $A_I$ is injective and $\operatorname{ERC}(I)>0$.`
- Rationale: Defines the ERC when I contains every column.

#### 25. Source line 771

- Before: `A sparse vector $x$ can be interpreted as the weights $x$ of a discrete measure $m_{x} \eqdef \sum_{i=1}^N x_i \de_{z_i}$, with atoms restricted to the sampling grid.`
- After: `The coefficients of a sparse vector $x$ are the weights of the discrete measure $m_{x} \eqdef \sum_{i=1}^N x_i \de_{z_i}$, whose atoms lie on the reconstruction grid.`

#### 26. Source line 828

- Before: `A single finite grid may miss an overshoot between samples. The major issue is that $\tilde\eta_F$ is only constrained by construction to interpolate $\sign(x_{0,i})$ at points $z_i$ for $i \in I$. Nothing therefore prevents $\tilde\eta_F$ from leaving $[-1,1]$ near an interpolation point.`
- After: `A finite grid may miss an overshoot between samples. The construction constrains $\tilde\eta_F$ only to interpolate $\sign(x_{0,i})$ at support points $z_i$; it does not prevent the function from leaving $[-1,1]$ nearby.`

#### 27. Source line 854

- Before: `The continuous precertificate $\tilde\eta_V=\Aa^*p_V$ is sampled on the grid to obtain $\eta_V$ from $\tilde\eta_V$.`
- After: `Evaluating $\tilde\eta_V=\Aa^*p_V$ on the reconstruction grid gives the discrete precertificate $\eta_V$.`

### compressed-sensing

Tightened the acquisition motivation and the distinction between measurement count and sparsity. Clarified noise bounds, strong and weak recovery thresholds, conditional independence, finite-sample concentration, RIP parameter domains, probability lower bounds, zero-vector cases, and sampling without replacement.

Coverage: the entire active chapter, including captions and proofs. Disabled blocks and substantive full-line comments match the baseline.

#### 1. Source line 6

- Before: `Combined with sparse regularization, such measurements can reconstruct signals from far fewer samples than their ambient dimension. Randomness also makes precise recovery and stability guarantees possible.`
- After: `For signals with sparse representations, suitable random measurements can retain enough information for reconstruction from far fewer samples than the ambient dimension. Sparse regularization then provides recovery and stability guarantees.`

#### 2. Source line 23

- Before: `The aim is to keep the measurement budget $P$ close to $M$ while retaining enough information to reconstruct $f_0$ accurately.`
- After: `The aim is to keep the measurement count $P$ close to the sparsity level $M$, typically with a logarithmic overhead, while retaining enough information for accurate reconstruction.`
- Rationale: Connects the acquisition motivation to the measurement scaling proved later, rather than suggesting P approximately M without qualification.

#### 3. Source line 31

- Before: `A single sensor integrates the reflected light during each exposure, producing the measurement $y_p\in\RR$.`
- After: `A single sensor integrates the light directed toward it during each exposure, producing one measurement $y_p\in\RR$.`

#### 4. Source line 33

- Before: `Let $c_q$ be the surface of mirror $q$, and define the ideal pixel intensity by $f_{0,q}=\int_{c_q}\tilde f_0(s)\,\mathrm{d}s$.`
- After: `Let $c_q$ be the region occupied by mirror $q$ in the focal plane, and define its integrated intensity by $f_{0,q}=\int_{c_q}\tilde f_0(s)\,\mathrm{d}s$.`

#### 5. Source line 79

- Before: `For this discussion, assume that the noise norm is known and set $\epsilon=\norm{w}$.`
- After: `Assume that a noise bound $\norm{w}\leq\epsilon$ is known, so that the true coefficient vector is feasible.`
- Rationale: Only an upper bound on the noise norm is needed, as in the later uniform recovery theorem.

#### 6. Source line 79

- Before: `Convex optimality conditions relate the two formulations through suitable Lagrange multipliers. This relation depends on $y$, may have flat segments, and need not be one-to-one.`
- After: `Convex optimality conditions relate the two formulations through suitable Lagrange multipliers. The parameter correspondence depends on $y$ and need not be one-to-one; a range of parameter values may produce the same solution.`

#### 7. Source line 81

- Before: `We therefore state the following results directly in terms of explicit distributional assumptions on $A$.`
- After: `The results below therefore impose distributional assumptions directly on $A$.`

#### 8. Source line 94

- Before: `\item All $x_0$ so that $\norm{x_0}_0 \leq C_A(P/N)P$ are identifiable.`
- After: `\item Every vector $x_0$ with $\norm{x_0}_0 \leq C_A(P/N)P$ is identifiable.`

#### 9. Source line 95

- Before: `\item Most $x_0$ so that $\norm{x_0}_0 \leq C_M(P/N)P$ are identifiable.`
- After: `\item A typical support/sign pattern with sparsity at most $C_M(P/N)P$ is identifiable.`
- Rationale: States the weak threshold in terms of typical support/sign patterns rather than an unspecified measure on all vectors.

#### 10. Source line 109

- Before: `Recovery guarantees rely on quantitative bounds for the singular values of random matrices.`
- After: `To prove recovery guarantees, we need quantitative control of the singular values of random column submatrices.`

#### 11. Source line 159

- Before: `\paragraph{Superlinear measurement growth.}`
- After: `\paragraph{Finite-sample concentration.}`
- Rationale: The bound is nonasymptotic and does not require a specific superlinear growth regime.

#### 12. Source line 170

- Before: `The singular-value bounds therefore control the Gram matrix as well.`
- After: `Thus the singular-value bounds also control the deviation of the Gram matrix from the identity.`

#### 13. Source line 180

- Before: `Besides an $\ell^2$ error bound, we seek exact recovery of the support when the noise is small.`
- After: `We seek both a small $\ell^2$ reconstruction error and exact support recovery at low noise levels.`

#### 14. Source line 192

- Before: `The columns of $A$, denoted $a_j \in \RR^P$, are normalized by $\norm{a_j}=1$, and the coherence is`
- After: `Assume that each column $a_j\in\RR^P$ has unit norm. The coherence is`

#### 15. Source line 274

- Before: ``We consider here a class of sufficiently ``random'' matrices.``
- After: ``For the probabilistic analysis, we specify independent entries with a common variance and a uniform Gaussian-type tail bound.``

#### 16. Source line 303

- Before: `Condition on $A_I$. The vector $p_F$ is then fixed, while the columns $a_j$ for $j\notin I$ remain independent of it.`
- After: `Conditioning on $A_I$ fixes $p_F$, while the off-support columns $a_j$ remain independent of this vector.`

#### 17. Source line 358

- Before: `The restricted isometry constant $\delta_s$ of a matrix $A\in\RR^{P\times N}$ is the smallest nonnegative number such that`
- After: `For a positive integer $s$, the restricted isometry constant $\delta_s$ of $A\in\RR^{P\times N}$ is the smallest nonnegative number such that`
- Rationale: Makes the sparsity parameter domain explicit.

#### 18. Source line 364

- Before: `for every $s$-sparse vector.`
- After: `Thus $A$ approximately preserves the Euclidean norm on every coordinate subspace of dimension at most $s$.`

#### 19. Source line 417

- Before: `it satisfies $\de_s \leq \de$ with probability $1 - 2e^{ -\de^2 \frac{P}{2C} }$, where $C$ depends only on the sub-Gaussianity parameters appearing in~\eqref{eq-sub-gaussian}.`
- After: `it satisfies $\de_s \leq \de$ with probability at least $1 - 2e^{ -\de^2 \frac{P}{2C} }$, where $C$ depends only on the sub-Gaussian parameters in~\eqref{eq-sub-gaussian}.`
- Rationale: The concentration result gives a lower probability bound, not an exact probability.

#### 20. Source line 421

- Before: `We outline the argument. % The RIP condition`
- After: `The proof combines concentration on each support with a union bound over supports. % The RIP condition`

#### 21. Source line 425

- Before: `For Gaussian matrices, Section~\ref{sec-random-matrix} gives the required estimates through the Gram matrix $B^*B \in \RR^{s \times s}$ of a matrix $B$ of size $(P,s)$, with $B \sim \texttt{randn}(P,s)/\sqrt{P}$.`
- After: `For Gaussian matrices, Section~\ref{sec-random-matrix} provides these estimates for $B^*B\in\RR^{s\times s}$, where $B\in\RR^{P\times s}$ has independent $\Nn(0,1/P)$ entries.`
- Rationale: Replaces programming-specific distribution notation with the mathematical entry distribution already assumed.

#### 22. Source line 466

- Before: `Let $c \in \RR^N$ with $\norm{c}_0 \leq s$ and $\supp(c)=J$.`
- After: `Let $c \in \RR^N$ with $\norm{c}_0 \leq s$ and $\supp(c)=J$, and let $s'\geq1$ be an integer.`
- Rationale: The off-support bound divides by sqrt(s prime), so a positive integer s prime is required.

#### 23. Source line 503

- Before: `Dividing by the restricted coefficient norm when it is nonzero gives`
- After: `If $\eta_L\neq0$, divide by $\norm{\eta_L}$; if $\eta_L=0$, the bound is immediate. In either case,`

#### 24. Source line 517

- Before: `Let $I_0=\supp(x_0)$. Apply Lemma~\ref{lemma-rip-eta}`
- After: `For $x_0=0$, take $p=0$. Otherwise, let $I_0=\supp(x_0)$ and apply Lemma~\ref{lemma-rip-eta}`
- Rationale: Handles the zero signal directly before constructing interpolants on a nonempty support.

#### 25. Source line 588

- Before: `The measurements are complex for a complex basis. Here $\Om\subset\{1,\ldots,N\}$ is drawn uniformly at random among all sets of size $P$.`
- After: `For a complex basis, the measurements are complex. The index set $\Om\subset\{1,\ldots,N\}$ is sampled uniformly without replacement, with $|\Om|=P$.`

#### 26. Source line 600

- Before: `with probability at least $1-N^{-\gamma}$, the normalized operator~\eqref{eq-random-subsample} satisfies $\delta_{2s}\leq c$.`
- After: `the normalized operator~\eqref{eq-random-subsample} satisfies $\delta_{2s}\leq c$ with probability at least $1-N^{-\gamma}$.`

#### 27. Source line 613

- Before: `Good recovery guarantees follow when the sampling and sparsity bases are incoherent.`
- After: `Small mutual coherence therefore yields stronger measurement bounds.`

#### 28. Source line 364

- Before: `Thus $A$ approximately preserves the Euclidean norm on every coordinate subspace of dimension at most $s$.`
- After: `A small $\delta_s$ therefore means that $A$ approximately preserves the Euclidean norm on every coordinate subspace of dimension at most $s$.`
- Rationale: Approximate preservation is informative when the restricted isometry constant is small; its definition alone does not impose this.

### File-by-file actual-edit record

#### `machine-learning.tex` — 42 recorded edits

1. **Editorial:** Name the centered matrix before using its factorization.
2. **Editorial:** Make the caption idiomatic and consistent with the surrounding PCA terminology.
3. **Editorial:** Remove promotional wording and identify what the displayed matrix represents.
4. **Editorial:** State the dataset dimensions and features directly.
5. **Mathematical clarification:** State the moment condition ensuring finite distortion in the continuous formulation.
6. **Mathematical clarification:** Qualify the geometry and handle boundary mass rather than silently integrating overlapping cells.
7. **Mathematical clarification:** Prevent division by zero in the continuous centroid formula.
8. **Mathematical clarification:** State independence before invoking consistency and avoid implying convergence without assumptions.
9. **Mathematical correction:** Match the strict-threshold probability defined immediately above; Markov still applies.
10. **Editorial:** Make the two penalties readable and qualify feature selection by the model structure.
11. **Editorial:** Use the standard name of the method.
12. **Mathematical clarification:** Replace the informal averaging description and state the almost-everywhere meaning of conditional expectation.
13. **Mathematical clarification:** Distinguish the population conditional expectation from its empirical version.
14. **Notation correction:** Use the outer product, rather than an undefined square of a vector, in the moment explanation.
15. **Editorial:** Remove a redundant qualifier and improve agreement.
16. **Editorial:** Use the standard adjectival form.
17. **Editorial:** Replace a vague superlative and an unconditional performance claim with the relevant tradeoff.
18. **Editorial:** Avoid switching the neighbor-count notation and remove subjective wording.
19. **Editorial:** Explain nonparametric structure precisely and specify ties.
20. **Editorial:** Remove the vague antecedent and misleading implication that the query belongs to the training set.
21. **Editorial:** Distinguish holdout from cross-validation and count a fraction rather than a number.
22. **Notation correction:** Reserve k for the number of classes and use R consistently for neighbors.
23. **Notation correction:** Align the caption with the neighbor-count notation.
24. **Notation clarification:** Identify the fitted predictor and remove nested explanatory wording.
25. **Editorial:** Remove repetitive contrast and state the accepted terminology.
26. **Mathematical clarification:** Handle sign at zero explicitly.
27. **Editorial:** Explain the score-to-probability map directly.
28. **Notation clarification:** Define s at its first use and identify the negative class.
29. **Mathematical correction:** Use the linear score in the margin; f currently denotes a probability predictor.
30. **Mathematical clarification:** Specify the class ordering needed for the sign of the parameter difference.
31. **Notation correction:** Replace the undefined probability-vector name h with the established model.
32. **Notation correction:** Align the caption with the probability predictor.
33. **Mathematical clarification:** Specify the real-valued feature setting used by the inner-product formulas.
34. **Editorial:** Use direct and idiomatic feature-map terminology.
35. **Mathematical clarification:** State the positive regularization required by all subsequent inverses.
36. **Mathematical correction:** Correct the claimed substitution: Gram and covariance matrices have different dimensions.
37. **Notation correction:** Use the defined kernel κ, and retain the dependence on responses and regularization.
38. **Editorial:** Remove unnecessary words before the prediction identity.
39. **Editorial:** Give the feature map a clear antecedent and remove a vague qualifier.
40. **Mathematical clarification:** State symmetry explicitly; nonnegative eigenvalues alone are insufficient for a nonsymmetric matrix.
41. **Mathematical clarification:** Exclude zero-norm feature vectors before the normalization division.
42. **Editorial:** Condense the transition to kernel distances.

#### `sec-pca-theory.tex` — 8 recorded edits

1. **Mathematical clarification:** Explain why centering also covers affine reconstruction spaces.
2. **Editorial:** Match the optimization variables to linear maps after centering.
3. **Editorial:** Distinguish the basis matrix from its columns and sharpen the proof transition.
4. **Editorial:** State the reduction without nested filler clauses.
5. **Mathematical correction:** Put every residual-expansion term under the summation; previously i was free in two terms.
6. **Notation clarification:** Use the full property implied by B transpose B equals identity.
7. **Editorial:** Identify the displayed bound without referring to a nonexistent rightmost inequality.
8. **Editorial:** Remove a dangling participle in the concluding proof.

#### `machine-learning-sec-pac.tex` — 10 recorded edits

1. **Editorial:** Clarify that the data form a sample, and state the statistical objective directly.
2. **Editorial:** Remove duplicated qualifications and improve the transition to surrogate losses.
3. **Editorial:** Match the introduction to the subsequent if-and-only-if statement.
4. **Terminology correction:** The displayed identity decomposes excess risk, not squared-error bias and variance.
5. **Editorial:** State the tradeoff concisely and remove an unclear antecedent.
6. **Editorial:** Remove vague lifting language and redundant parentheticals.
7. **Mathematical clarification:** Avoid implying that the target is already inside the radius-R class, and state the integrability needed by the bound.
8. **Mathematical clarification:** Qualify the rate when the radius may grow with n.
9. **Editorial:** Give the representer argument without suggesting differentiability of hinge loss is required.
10. **Editorial:** Replace awkward phrasing by a precise account of the reduction.

#### `optim-ml-smooth.tex` — 0 recorded edits

Read the complete wrapper; chapter title, label, and active inputs are appropriate. No change was needed.

#### `sec-optim-smooth.tex` — 44 recorded edits

1. **Editorial:** Avoid a dangling relative clause around the displayed definition.
2. **Editorial:** Complete the revised least-squares sentence naturally.
3. **Mathematical clarification:** Separate constrained attainment from the coercive reduction and include nonemptiness.
4. **Mathematical clarification:** Give both directions of the equivalence without repetitive eigenvalue phrasing.
5. **Mathematical clarification:** Explain explicitly why the stated counterexample is not differentiable.
6. **Editorial:** Make the two applications in the converse proof readable.
7. **Editorial:** Correct the placement of respectively and identify the weights.
8. **Editorial:** Use an idiomatic definition of local minimizer.
9. **Mathematical clarification:** Specify positive perturbations before dividing by epsilon.
10. **Editorial:** Remove a repeated statement and inconsistent x/x-star antecedent.
11. **Terminology correction:** Match the one-dimensional cubic in the figure and the corrected surrounding text.
12. **Editorial:** Remove repeated vague distance wording from the local-to-global proof.
13. **Editorial:** Remove a duplicated phrase and comma splice, and define C before its first use in this subsection.
14. **Editorial:** Condense the eigendecomposition explanation while preserving normalization and column convention.
15. **Mathematical clarification:** State the positive-definite assumption needed by the following Gaussian density and inverse.
16. **Mathematical clarification:** Cover the zero-eigenvalue case excluded from the inverse/reciprocal formula.
17. **Mathematical correction:** Correct the probability interpretation: sigmoid of a negative margin is not generally the positive-class probability.
18. **Mathematical clarification:** Supply the regular-point hypothesis needed for a tangent-space interpretation.
19. **Mathematical clarification:** State local regularity and explain the meaning of the set-valued argmin limit.
20. **Mathematical clarification:** Replace the unsupported Lagrange-multiplier sign choice by a direct comparison proof.
21. **Mathematical clarification:** Do not claim that stationarity characterizes a global line minimum without convexity; ensure the optimum is interior.
22. **Editorial:** Distinguish the gradients from their negative descent directions.
23. **Mathematical clarification:** Require positive backtracking factors and avoid prescribing an unnecessarily enormous initial step.
24. **Editorial:** Make the quadratic convergence assumptions grammatical and explicit.
25. **Mathematical clarification:** Distinguish a worst-case spectral bound from optimal nonconstant step schedules.
26. **Mathematical correction:** Use one step-size variable and distinguish the actual spectral norm from a bound using possibly nonsharp endpoints.
27. **Mathematical correction:** Supply the compact-interval uniform contraction argument and qualify endpoint equality.
28. **Mathematical correction:** Evaluate the uniform bound, which equals the norm only for sharp eigenvalue bounds.
29. **Editorial:** Identify exactly which quantity is being optimized.
30. **Editorial:** Remove a panel reference from a single-panel figure.
31. **Mathematical clarification:** Remove duplicated conditional clauses and exclude k=0 from the 1/k bound.
32. **Mathematical correction:** Do not equate equally indexed sorted singular values after a nonmonotone spectral transformation.
33. **Mathematical clarification:** Explain the spectral mapping without an incorrect ordering assumption.
34. **Editorial:** Avoid describing a coordinate as a direction and streamline the Hessian definition.
35. **Mathematical correction:** Restore the Taylor remainder; the expansion is not an exact affine identity for logistic loss.
36. **Editorial:** Repair an unfinished sentence and stray capital letter.
37. **Editorial:** Use eigenvalues for the semidefinite-order interpretation.
38. **Editorial:** Use the standard order of matrix adjectives.
39. **Editorial:** Remove speculative cycling language and introduce the completed proof.
40. **Editorial:** Avoid claiming an unspecified extension before supplying it.
41. **Notation correction:** Match the sequence index.
42. **Mathematical proof completion:** Complete the previously omitted convergence-of-iterates argument and variable-step objective/linear bounds using cocoercivity, Fejer monotonicity, and telescoping.
43. **Mathematical clarification:** Fix a positive smoothness bound before inverse-L step formulas.
44. **Layout within proof:** Separate the two new convergence inequalities to preserve readable line lengths.

#### `optim-ml-advanced.tex` — 0 recorded edits

Read the complete wrapper; chapter title, label, and active inputs are appropriate. No change was needed.

#### `sec-regul.tex` — 12 recorded edits

1. **Editorial:** Make explicit that the final comparison applies to every feasible point.
2. **Mathematical clarification:** Spell out the accumulation-point and arbitrary-comparator steps.
3. **Mathematical clarification:** Include the square case and avoid equating full rank with strict overdetermination.
4. **Mathematical clarification:** Give the precise rank condition rather than an imprecise problem-size label.
5. **Editorial:** State the universally valid SVD construction and standard terminology.
6. **Editorial:** Use the objective notation from the definition and specify the differentiability obstruction.
7. **Notation correction:** Replace placeholder dots by the argument in a function-value definition.
8. **Mathematical correction:** Include the entire coordinate penalty under the sum.
9. **Mathematical clarification:** State the response-sign case and justify zero as minimizer rather than drawing a conclusion from the positive half-line alone.
10. **Mathematical clarification:** Specify positivity before using its reciprocal.
11. **Mathematical clarification:** Explain why the argmin update is single-valued; cancellation leaves curvature 1/tau.
12. **Mathematical clarification:** Use the operator norm explicitly, require a positive fixed step, and handle the zero-operator denominator.

#### `sec-stochastic-optim.tex` — 11 recorded edits

1. **Editorial:** Remove repetition.
2. **Editorial:** Distinguish exact finite-sum representation from Monte Carlo approximation.
3. **Editorial:** Use idiomatic word order.
4. **Mathematical clarification:** Include the assumptions omitted from the step-size claim.
5. **Editorial:** Separate scores from labels in the derivative notation.
6. **Mathematical clarification:** Explain both tuning parameters and require valid positive values.
7. **Mathematical clarification:** Define the filtration supporting conditional unbiasedness.
8. **Editorial:** Repair the sentence fragment and identify what is averaged.
9. **Mathematical clarification:** Explain why the SGD unbiasedness proof cannot be reused for SAG.
10. **Editorial:** Remove an unmatched closing parenthesis.
11. **Caption correction:** Correct the left/right panel description after inspecting both actual assets.

#### `sec-autodiff.tex` — 19 recorded edits

1. **Editorial:** Explain the numerical limitation motivating automatic differentiation.
2. **Mathematical clarification:** Define dimensions rather than calling a vector its dimensionality, and distinguish blocks from scalar coordinates.
3. **Mathematical clarification:** Match the vector-block initialization by an identity matrix.
4. **Editorial:** Avoid undefined x_p block notation when s differs from the total scalar dimension p.
5. **Mathematical clarification:** State the invertibility condition for dual numbers.
6. **Editorial:** Correct agreement with the imperative initialize.
7. **Notation correction:** Match the function signature in the worked example.
8. **Mathematical clarification:** Do not infer the comparison from input/output dimensions alone for arbitrary intermediate bottlenecks.
9. **Editorial:** Correct agreement with the imperative initialize.
10. **Mathematical clarification:** State differentiability limitations for common nonsmooth activations.
11. **Mathematical clarification:** Describe a mesh-refinement limit without an ambiguous sequence-to-moving-target convergence statement.
12. **Editorial:** Complete the sentence after the adjoint gradient normalization.
13. **Mathematical clarification:** Avoid conflating one gradient step with an exact argmin layer and avoid overloading the scalar potential by a separable vector function.
14. **Editorial:** Streamline the transition to the reversible update without changing the scheme.
15. **Mathematical clarification:** State the continuity used in passing to the limit.
16. **Mathematical clarification:** Address the single-valued map implicit in the argmin notation.
17. **Mathematical clarification:** Supply the local smoothness required by implicit differentiation.
18. **Editorial:** Explain the essential update order directly without an unnecessary alternative algorithm.
19. **Mathematical clarification:** Avoid claiming that scalar output alone guarantees the smaller dense matrix-chain cost. For widths (2,100,1,1), the displayed forward/reverse counts are 202 and 300.

### Preface

- Reread the complete preface and retained its already concise scope, roadmap, and Numerical Tours reference.

### Shallow Learning (`chapters/perceptrons.tex`)

- Replaced the vague “Multi-layer” subsection title with “Networks of Arbitrary Depth.”
- Defined the network output, its dimension, and the full parameter tuple; explained that the final activation may be omitted for a linear output layer.
- Used `σ` consistently in the activation examples, reserving `ρ` for data distributions; attributed positive homogeneity specifically to ReLU.
- Connected the two-layer formula explicitly to one hidden layer and a linear output layer.
- Defined the joint input–output law in the quadratic loss and explained its empirical-distribution specialization.
- Stated that observations are matrix columns, distinguished the sample count from the neuron count, and explained omission of the constant empirical-loss normalization.
- Added the differentiability assumption needed for the displayed parameter-gradient calculation.
- Repaired the grammar of the universal-approximation proposition and made “nondecreasing” consistent with the allowed discontinuous sigmoids.
- Corrected citation spacing and explained why the proof first uses cosine networks before approximating them with the prescribed sigmoid.
- Made the local Fourier normalization explicit and standardized the capitalization of Barron space.
- Corrected “Sobolev norm” to “squared Sobolev norm” in its comparison with the sum of squared derivative norms.
- Corrected the direction of the coefficient rescaling when rewriting a sum of neurons as an average: the outer coefficients absorb a factor of the width.
- Specified the pointwise supremum norm for bounded functions paired with arbitrary measures; this retains values at atoms and avoids an inappropriate essential-supremum interpretation.
- Added convexity of the domain to the quadratic-upper-bound lemma.
- Stated the Frank–Wolfe step schedule from iteration zero, consistently with the proof's use of the first step equal to one.
- Rewrote the first-variation Lipschitz proposition as complete sentences, restored the established neuron-parameter notation, and stated the bound on its outer coefficient.

### Deep Learning (`chapters/deep-learning.tex`)

- Described the architecture as a composition of parameterized maps.
- Made explicit that pointwise activations preserve dimension.
- Qualified the dense-network operation count by its assumption of pointwise activations.
- Replaced a misleading chronological transition with a direct connection to the earlier logistic models.
- Replaced ambiguous superscript/subscript indexing of convolution filters with explicit input- and output-channel ranges.
- Simplified the translation-equivariance proof by naming each filter directly as the impulse response; connected adjoint filter reversal to the backpropagation formula.
- Made the channel notation consistent between the input representation and the convolution formula.
- Replaced an equality between a bias vector and a scalar with the correct componentwise statement that each channel's bias is spatially constant.
- Replaced an unsupported claim that a residual bottleneck “enforces” regularization with its precise parameter-count and dimensional consequences.
- Stated that attention weights sum to one over keys for each fixed query index.

### Convex Analysis (`chapters/convex-analysis.tex`)

- Reordered the introduction to lead with the chapter's purpose and integrated its notation and references.
- Specified the convention for zero times positive infinity at the endpoints of the extended-real convexity inequality.
- Corrected “denote by” and standardized subdifferential terminology, first-order headings, and paired-name dashes.
- Simplified the explanation of the dual space while preserving the distinction between canonical pairing and inner-product identification.
- Removed repetition in the empty-subdifferential example and defined the extended function by cases in prose, avoiding a square-root expression outside its real domain.
- Wrote the nonzero branches of the absolute-value subdifferential as singleton sets.
- Stated the closed, convex, nonempty set assumptions for the normal cone.
- Gave the precise continuity-at-a-feasible-point qualification for the constrained sum rule.
- Introduced “convex conjugate” as an alternative name for the Legendre–Fenchel transform.
- Separated the Fenchel–Young inequality from its equality characterization.
- Defined the dimensions of the matrix and right-hand side in affine equality duality.
- Used an infimum for the general primal value, since attainment has not yet been assumed.
- Separated the definitions of the Lagrangian, the dual function, and the primal/dual values; repaired chains that incorrectly equated an optimized value to a function of free variables.
- Completed the weak-duality proof by explicitly taking the primal infimum and dual supremum.
- Repaired the grammar of the strong-duality theorem and specified that Slater's strict inequality is componentwise.
- Explained complementary slackness as multiplier maximization at the feasible primal point.
- Distinguished the role of convexity in KKT sufficiency from the role of a constraint qualification in necessity.
- Clarified the auxiliary-variable reformulation and grouped both inner infima inside the outer supremum in the Fenchel–Rockafellar derivation.

### Nonsmooth Convex Optimization (`chapters/optim-nonsmooth.tex`)

- Removed a comma splice in the references and clarified the smooth-objective introduction.
- Distinguished the gradient of a functional from the spatial gradient of a signal/image; removed the incorrect assertion that the latter is the gradient “of f” and belongs to the same scalar-image space.
- Explained the quadratic Taylor remainder and identified the gradient-computation cross-reference as a section.
- Repaired punctuation after gradient and projected-gradient updates.
- Renamed “Sub-gradient Descent” to “Subgradient Method” and explained that a subgradient step need not decrease the objective.
- Made positivity of the bounds explicit in the finite-horizon subgradient step formula.
- Distinguished continuous differentiability from the Lipschitz-gradient assumption in projected-gradient convergence.
- Stated a positive Lasso parameter and supplied the missing zero right-hand side in the nonnegative-variable formulation.
- Distinguished convergence of central-path objective values from convergence of the minimizers themselves.
- Stated the differentiability/invertible-Hessian assumptions for Newton's method and the need for its line search to preserve strict feasibility.
- Evaluated the Hessian and gradient at the same iterate in the Newton stopping criterion.
- Removed a trailing space from the proximal-map heading.
- Reworked the soft-thresholding proof with a scalar input and a distinct proximal output, avoiding collisions with the quadratic-data notation.
- Repaired the sentence introducing nonconvex hard thresholding while preserving its set-valued threshold ties.
- Described the proximal composition formula through its exact orthonormal-row assumption, avoiding ambiguity between a frame's analysis and synthesis operators.
- Corrected the reversed explanation of equality versus membership for a single-valued resolvent.
- Removed a repeated resolvent definition and standardized maximally monotone terminology.
- Specified that the clipping formula acts componentwise.
- Standardized Moreau–Yosida, forward–backward, and Douglas–Rachford headings and names.
- Shortened the surrogate derivation and the explanation recovering projected gradient as a special case.
- Rewrote the dual-splitting introduction, explicitly stated finiteness/continuity, and described when the composed proximal map has an explicit formula.
- Corrected citation spacing and aligned the TV vector-field description, sample count, and mixed-norm definition.
- Stated positivity of the SVM regularization parameter.
- Replaced an insufficiently qualified variable dual step-size prescription by a fixed admissible step, matching the cited convergence theorem; distinguished a dual maximizer from a minimizer of its negative objective.
- Rewrote the constrained-Lasso example and graph-projection explanation into shorter sentences; explained the choice between input-space and output-space linear systems.
- Rewrote the ADMM introduction to identify the relevant dual proximal maps and their relation to alternating primal minimization.
- Separated Lagrangian definitions from saddle values in the ADMM derivation and explained precisely why adding the quadratic penalty preserves the constrained problem.
- Identified the joint augmented-Lagrangian update as the method of multipliers.
- Simplified the multiplier-step explanation and stated its stationarity consequence at the current iterates.
- Rewrote the primal–dual introduction around biconjugacy.
- Corrected the TV example's expression for the function value, replacing a placeholder dot by its argument, and repaired the surrounding punctuation/capitalization.

All existing labels, citations, figure inclusions, and comparison figures were retained. The mathematical adjustments above follow directly from the definitions and derivations in the chapters; no figure data or historical attributions were reconstructed in this editorial pass.

### Final integration corrections

- Aligned the quadratic gradient-descent proof with its iteration index: the update and error recurrences now use `τ_k` beside `x_k`, and the generic spectral bound is explicitly a function of `τ`.
- Completed the sentence introducing the auxiliary-variable reformulation in Fenchel–Rockafellar duality.
- Removed trailing whitespace on the 11 affected newly edited lines, without reformatting untouched source material.
- Independently checked the revised full-range gradient-descent proof, spectral orthogonalization proof, and the convex-duality/proximal/ADMM changes against the surrounding definitions and algebra.

### Verification of this pass

- Rebuilt and published `FundationsDataScience.pdf` (**339 pages**) and all **19 standalone PDFs** in `chapters-pdf/`. All 20 final builds have **zero LaTeX/BibTeX warnings, missing characters, or overfull/underfull boxes**.
- Verified embedded fonts and resolved destinations throughout the 20 published documents: **1,833 internal links**, **150 chapter-to-book links**, and **725 shared label numbers**, with no issues.
- Confirmed that all **91 editable TikZ reconstructions** remain current single-page vector PDFs with embedded fonts and no raster objects or text beyond their page boundaries.
- Matched all **91 comparison pages** in the book and the corresponding **91 comparisons** across the standalone chapters to their manifests. Original/new panel labels and landscape orientation are intact.
- Rendered and visually inspected **all 339 book pages**, with enlarged checks of dense proofs, revised mathematical passages, and comparison pages. No clipping, overlaps, or broken panels were found. Chapter-ending figure pages retain their deliberate whitespace.
- Inspected **38 standalone chapter first/last pages** and **26 additional comparison-boundary pages**; titles, margins, chapter numbering, references, and comparison layouts are clean.
- Preserved existing label definitions and citation commands in all 26 reviewed source files. The final audits and page renders are retained in `build/build-report.json`, `build/pdf-audit.json`, `build/label-number-audit.json`, `build/tikz-pdf-audit.json`, `build/editorial-pass3-reviews/`, and `build/qa/editorial-pass3/`.

## Chapter introductions and publication details

Added a short opening overview to every one of the 19 active chapters. Each overview has three sentences (51–67 words) explaining the central problem, its importance, and the chapter's scope. The shared `\chapteroverview` format places this text beneath the chapter heading in the book's sans-serif typeface and dark ink color.

Overlapping opening prose was absorbed or shortened in Shannon sampling, Fourier analysis, approximation, compression, denoising, inverse problems, sparse-recovery theory, compressed sensing, machine learning, shallow learning, and convex analysis. Technical setup, bibliography references, labels, and formal mathematical content were retained. The two optimization wrappers now contain their own overview before any section input.

### Shared book identity and date

- Added `book-metadata.tex` as the common source for the book title, Gabriel Peyré's name, the affiliation **CNRS & DMA, École Normale Supérieure**, and the edition date. The affiliation is the one already used on the cover.
- The standalone chapter masthead displays the book title, author, and affiliation above the chapter number and title, with distinct type sizes and spacing. It appears on the chapter opening rather than on a separate title page.
- The date appears on the book cover and at the start of every numbered chapter in both editions. It uses the build date through `\today`; it can be fixed for a dated edition in the shared metadata file.
- The cover and PDF title/author metadata use the shared definitions. The cover keeps the book title on three lines, with “Data Sciences” together.

### Chapter opening texts

#### Shannon Sampling Theory

Source: `chapters/shannon.tex`.

> Digital measurements retain only samples of a continuous signal, so reconstruction depends on what information sampling preserves. Understanding these limits guides the choice of sampling rate and the filtering applied before acquisition. Starting from continuous and discrete signal models, we use Fourier analysis to prove a sampling theorem, explain aliasing, and distinguish exact interpolation from the additional error introduced by quantization.

#### Fourier and Convolution

Source: `chapters/fourier.tex`.

> Convolution combines shifted copies of a signal, but evaluating it directly can be costly and obscure the structure of the operation. Fourier coordinates turn convolution into multiplication, providing both efficient algorithms and a way to solve differential equations. We develop this connection on continuous and finite domains, study sampling and the fast Fourier transform, and extend the viewpoint to groups, surfaces, and graphs.

#### Shannon Coding Theory

Source: `chapters/shannon-coding.tex`.

> A sequence of symbols can often be stored using fewer bits when its probabilities and dependencies are taken into account. The challenge is to reduce the expected length while keeping every message unambiguously recoverable. Prefix codes and entropy lead to Shannon's source coding bound; Huffman's construction, block coding, and models of dependent symbols then show how coding methods approach that bound.

#### Wavelets

Source: `chapters/wavelets.tex`.

> Signals contain broad trends and localized details that occur at different scales. Separating these components produces representations suited to compression and denoising, while allowing each scale to be processed efficiently. We build orthonormal wavelet bases from nested approximation spaces, derive fast decomposition and reconstruction algorithms in one and two dimensions, and explain how filter design controls localization, smoothness, and cancellation of polynomials.

#### Linear and Nonlinear Approximation

Source: `chapters/approximation.tex`.

> When only a limited number of coefficients can be retained, the choice of representation determines which features of a signal survive. Approximation rates quantify this tradeoff and help explain the performance of compression and denoising methods. We compare fixed and adaptive coefficient selection, relate Fourier and wavelet errors to smoothness and edges, and examine triangulations and curvelets for images with smooth contours.

#### Compression

Source: `chapters/compression.tex`.

> Compression must turn a small set of transform coefficients into a short binary description without introducing excessive reconstruction error. This requires controlling both the loss from quantization and the cost of recording coefficient values and locations. We derive error and bit-rate bounds for transform coding, explain how probability models reduce coding costs, and describe the wavelet and entropy-coding stages of JPEG-2000.

#### Denoising

Source: `chapters/denoising.tex`.

> Denoising estimates a clean signal from measurements corrupted by noise. Its central difficulty is to suppress random fluctuations without erasing edges, textures, or other features of the underlying signal. We model common acquisition noise, compare linear filters with nonlinear thresholding, analyze how regularity and sparsity govern reconstruction error, and adapt the methods to Poisson and multiplicative noise.

#### Variational Priors and Regularization

Source: `chapters/variational-priors.tex`.

> Variational methods reconstruct images by balancing agreement with measured data against a penalty describing expected image structure. The choice of penalty controls the tradeoff between noise removal and preservation of sharp edges. Starting from Sobolev and total variation energies, we develop their discrete operators, derive the associated diffusion flows, and compare stopping a flow early with minimizing an energy that includes a data fidelity term.

#### Inverse Problems

Source: `chapters/inverse-problems.tex`.

> Inverse problems seek to reconstruct a signal from measurements that blur, omit, or mix its values. Direct inversion can amplify noise or leave several signals indistinguishable, so reconstruction must incorporate assumptions about the unknown signal. We study how spectral and variational regularization stabilize inversion, derive error bounds for quadratic penalties, and develop methods for deconvolution, inpainting, and tomography.

#### Sparse Regularization

Source: `chapters/sparse-regularization.tex`.

> Many signals admit accurate approximations using only a small number of coefficients in a suitable basis or dictionary. Sparse regularization turns this observation into a reconstruction principle for noisy or incomplete measurements. We move from counting nonzero coefficients to a convex penalty, distinguish analysis and synthesis models, derive thresholding and iterative algorithms, and apply them to deconvolution and inpainting.

#### Theory of Sparse Regularization

Source: `chapters/sparse-theory.tex`.

> The theory of sparse recovery asks when a reconstruction is unique, how much noise changes its coefficients, and whether it identifies the correct nonzero locations. A small reconstruction error alone does not guarantee that these locations are recovered. We use convex geometry and dual certificates to establish recovery conditions for the Lasso and its constrained counterpart, then examine their implications for spike deconvolution on increasingly fine grids.

#### Compressed Sensing

Source: `chapters/compressed-sensing.tex`.

> Compressed sensing aims to recover sparse signals from far fewer measurements than their ambient dimension. This can reduce acquisition costs when measurements are expensive, but success depends on how the sensing operator interacts with the signal representation. Motivated by the single-pixel camera, we analyze recovery from random measurements, distinguish guarantees for fixed signals from uniform guarantees, and extend the discussion to structured Fourier sampling.

#### Basics of Machine Learning

Source: `chapters/machine-learning.tex`.

> Machine learning seeks useful structure in data and predictors that remain accurate on new observations. We begin with dimensionality reduction and clustering, then develop regression and classification through empirical risk minimization, regularization, and kernel methods. Connections with inverse problems explain the role of model assumptions, while the final analysis of generalization balances approximation error against the statistical complexity of the chosen model class.

#### Optimization & Machine Learning: Smooth Optimization

Source: `chapters/optim-ml-smooth.tex`.

> Training a model often reduces to minimizing a differentiable objective, so the geometry of that objective determines how optimization algorithms progress. Starting from regression and classification, we develop convexity, derivatives, and optimality conditions, then derive gradient descent and practical step-size rules. Quadratic examples guide the convergence analysis, which extends to smooth convex objectives and explains the effects of conditioning, strong convexity, and acceleration.

#### Optimization & Machine Learning: Advanced Topics

Source: `chapters/optim-ml-advanced.tex`.

> Reliable and efficient learning requires control of model complexity, economical use of large datasets, and accurate derivative calculations. We study regularization through ridge regression and sparsity, followed by stochastic gradient methods that reduce the cost of individual updates. The final section derives automatic differentiation on computational graphs and applies it to feedforward networks, recurrent models, and optimization procedures.

#### Shallow Learning

Source: `chapters/perceptrons.tex`.

> Networks with one hidden layer provide a setting in which to ask which functions neural models can represent and how their size controls approximation error. These questions distinguish expressive power from the difficulty of fitting a model to data. We derive the training gradients, prove universal approximation, and develop quantitative bounds through Barron's theorem, a representation by measures, and the Frank--Wolfe algorithm.

#### Deep Learning

Source: `chapters/deep-learning.tex`.

> Deep networks learn representations through successive transformations of the input. Their architecture controls which structures they can exploit and how efficiently they can be trained, making it essential to understand the operations within each layer. We develop fully connected and convolutional models and their gradients, then introduce residual connections, normalization, attention, and wavelet scattering.

#### Convex Analysis

Source: `chapters/convex-analysis.tex`.

> Convexity turns local optimality conditions into certificates of global solutions, even when objectives have corners or encode constraints. This provides a common language for understanding regularization and designing optimization methods. We develop subgradients and normal cones, study convex conjugates and their relation to smoothness, and derive dual problems and optimality conditions.

#### Nonsmooth Convex Optimization

Source: `chapters/optim-nonsmooth.tex`.

> Many recovery and learning problems combine smooth losses with constraints or penalties that are not differentiable. Solving them efficiently requires understanding both the cost of an iteration and the conditions for convergence. We begin with subgradient, projection, and barrier methods, then develop proximal maps and splitting algorithms, including ADMM and primal--dual schemes.

### Verification

- Rebuilt and published the **341-page** `FundationsDataScience.pdf` and all **19 standalone chapter PDFs**. All 20 final documents compile with **zero warnings, missing characters, or overfull/underfull boxes**. A final cover-only rebuild restored the intended three-line title without changing the chapter editions.
- Verified that every chapter has exactly one three-sentence overview before its technical material. The complete overview is present on each chapter's opening page in both editions.
- Verified the book title, author, affiliation, and **September 5, 2026** date on the cover; checked the title/author/affiliation masthead in all standalone chapter openings and the date before all 38 chapter overviews across the two editions.
- Rendered and visually inspected the cover and all **38 chapter-opening pages**, including the long optimization titles at full resolution. Headings, dates, overviews, and technical text are complete and have no clipping or overlaps.
- Rechecked embedded fonts, **1,833 internal links**, **150 chapter-to-book links**, and **725 shared label numbers** with no issues. All **91 original/TikZ comparisons** remain present in the book and in the corresponding standalone chapters.
- Validation records are in `build/build-report.json`, `build/chapter-overviews-reviews/source-audit.json`, `build/chapter-overviews-reviews/pdf-audit.json`, `build/pdf-audit.json`, `build/label-number-audit.json`, `build/tikz-pdf-audit.json`, and `build/qa/chapter-overviews/`.

## Second TikZ figure pass: mathematical purpose and visual clarity

Completed a contextual review of all **91 reconstructions** on **2026-09-05**. Each drawing was compared with its original scan, its editable source, and the surrounding definitions, formulas, or proof. **51 drawing sources were refined and 40 were retained after review**; one retained drawing also received a correction to its interpretation note. All original assets and the original/new comparison structure are preserved.

The main corrections make sample grids and Fourier transforms consistent, repair projected geometry and kernel scaling, distinguish conditional means from modes, expose the constraints used in the PCA proof, and align optimization arrows, tangents, and network operations with their definitions. The individual decisions are recorded below.

### Sampling and quantization

#### Sampling and spectral overlap (`shannon-sampling-aliasing`)

Drawing refined. Made the lost out-of-band mass and its in-band alias contribution visible with red fills; superimposed the original spectrum for a direct comparison and added exact-recovery/failure conclusions. Linked the supporting sampling theorem and Poisson formula.

#### Uniform scalar quantization (`shannon-quantizer`)

Drawing refined. Added a vertical error segment at a fixed normalized input and the resulting T/2 reconstruction bound. Preserved half-open cells and index/dequantization distinction.

Retained after contextual and visual checks:

- **The cardinal sine** (`shannon-sinc`). Continuous value at zero, all integer zeros, normalization and oscillation signs are correct and readable.
- **Cardinal B-splines of degrees zero, one, and two** (`shannon-splines`). Convolution-defined supports/heights, quadratic pieces and box endpoint conventions are correct; no additional annotation is needed.
- **Two frequencies with identical samples** (`shannon-aliasing`). Both cosines agree at every integer and the spectral folding arrows use the correct signed shifts of 2 pi.


### Fourier analysis

#### Convolution with a box window (`fourier-convolution`)

Drawing refined. Checked convolution equation `eq-convol-1d` and the scanned box-window sketch. The shaded integral correctly covered `[x-epsilon,x+epsilon]`, but the heading immediately above that interval defined the unshifted box as a function of `t`. It now names the actual factor `g(x-t)` and its support; a separate small definition retains the unshifted, unnormalized box. The original signal and shaded geometry are unchanged. The resulting convolution remains an integral, not a local average. Rendered labels, guides, bracket, and integral have clear separation.

#### One radix-two FFT step (`fourier-fft`)

Drawing refined. Checked the decimation-in-frequency recursion immediately before `fig-fft`. The difference branch must be multiplied componentwise by the vector of twiddle factors before its half-size transform. Changed the ambiguous multiplication sign in that block to the chapter's circle-dot notation. Sum/difference definitions, negative exponential sign, transform sizes, and even/odd interleaving are unchanged and correct. The previously repaired odd-frequency label remains clear of its block and connector.

#### Spatial zero padding refines frequency sampling (`fourier-padding-spatial`)

Drawing refined. The previous schematic did not have one consistent sample count: the last spatial marker coincided with the point labeled `T=Q/N`, whereas the grid ends at `T-1/N`, and the frequency marker count did not follow from the physical sample grid. The two continuous profiles were also unrelated illustrative curves.

Replaced these with a fully specified analytic example: `N=12`, `Q=24`, and `f(x)=sin²(pi x)` on `[0,1]`, zero elsewhere. There are now 12 original samples and 12 added zeros, with the final sample at `23/12` and an unoccupied endpoint at `T=2`. The frequency plot has exactly 24 signed bins, `k=-12,...,11`, at angular frequencies `pi k`; the positive Nyquist endpoint is not duplicated. The continuous curve is the magnitude of the true Fourier transform, while distinct red markers show the normalized DFT magnitudes. Added the spatial sample step, concrete sample counts, analytic signal definition, and a short legend. Retained the approximation sign, because a finite Riemann sum is not the continuous integral.

The figure now ties the endpoint, sample counts, frequency step, Nyquist interval, normalization, and actual transform pair together. It introduces no empirical data. The original scan and its intentionally qualified reconstruction status are retained. Detailed derivation and numerical checks appear below.

#### Spectral zero padding interpolates the samples (`fourier-padding-fourier`)

Drawing refined. Checked the chapter signed-frequency padding convention and inverse factor `Q/N`. The existing polynomial has five original samples, 15 fine-grid samples, and exactly five nonzero retained frequency magnitudes. The analytic values and scale factor are internally consistent. Added the already-used values `N=5, Q=15` to the drawing so readers can reconcile the visible sample counts with the generic formula. No curves, coefficients, sample positions, or formulas were changed. The title fits without colliding with the transform arrows or spectral heading.

#### A plane wave and its normal direction (`fourier-wave`)

Drawing refined. Checked the tensor-product plane-wave formula and constant-phase interpretation. The original marked base point was slightly displaced from the phase line that its right-angle marker purported to meet. Moved it onto the actual line, emphasized that line, and rebuilt the marker at the exact incidence. The phase-line tangent is proportional to `(-0.6,1)`; the displayed wave vector is proportional to `(4.2,2.52)`, whose inner product with that tangent is zero. Extended the arrow beyond the plotted stripe field so its omega label cannot be crossed by another phase line. The geometry indicates a normal direction, not an asserted physical propagation direction.

#### A periodic two-dimensional grid (`fourier-torus-discrete`)

Drawing refined. Checked that the displayed four-by-four grid represents a product of cyclic index sets, with the last distinct index followed by zero after one more step. Both wrap arrows preserve the intended coordinate directions and do not identify the first and last displayed vertices. The top coordinate label touched the red wrapping curve. Moved that label above the curve with clear space. No grid geometry, adjacency, modular convention, or formula changed.

#### Colatitude and azimuth on the sphere (`fourier-spherical-coordinates`)

Drawing refined. The old freehand point, equatorial projection, polar arc, and equatorial ellipse could not be obtained from one projected unit sphere. Replaced them with one explicit orthographic construction. The viewing elevation is 30 degrees; the illustrative point has colatitude 50 degrees and azimuth 38 degrees. Its radial vector, equatorial projection, both angle arcs, equator, and poles now use the same projection. The projected poles correctly lie inside the circular silhouette, and the front and back equator use solid and dashed strokes consistently. Added “orthographic view” to prevent reading screen angles as actual spatial angles. Moved angle, point, and pole labels away from curves during the close rendering check.

Theta remains colatitude from the positive north axis, as required by the chapter formula; the ambiguous original angle has not been silently reinterpreted as latitude. The manifest explains the deliberate resolution and the illustrative viewing parameters. The projection also passed an independent check.

Retained after contextual and visual checks:

- **Convolution becomes multiplication** (`fourier-convolution-fourier`). Checked `eq-convol-fourier`, all four arrows, domains, and the chapter Fourier normalization.
- **Cardinal splines by convolution** (`fourier-splines`). Checked the exact centered box, triangular spline, and quadratic spline against their successive convolutions.
- **Translation equivariance** (`fourier-translation-inv`). Checked the definition `T_tau f(x)=f(x-tau)` and the commuting identity `H T_tau=T_tau H`.
- **Fourier modes align at the Dirac comb** (`fourier-poisson-distrib`). Checked the chapter distributional Poisson identity.
- **A square with opposite edges identified** (`fourier-torus`). Checked the product quotient `(R/2pi Z)²`.
- **Heat diffusion is Gaussian convolution** (`fourier-heat`). Checked the Fourier multiplier `exp(-t omega²)`, the one-dimensional Gaussian normalization, and variance `2t`.
- **Continuum and discrete Laplacian spectra** (`fourier-finite-difference-spectrum`). Checked `eq-disc-diff-2` and `eq-fft-lapl`.
- **A small geodesic ball on a surface** (`fourier-laplacian-surface`). Checked the small-geodesic-ball mean-value expansion, including the division by ball volume and coefficient `epsilon²/[2(d+2)]`.
- **A weighted graph and its local Laplacian** (`fourier-weighted-graph`). Checked the local adjacency drawing, symmetric nonnegative edge weights, and chapter sign `Delta=W-D`.


### Source coding

#### Kraft inequality: disjoint descendant blocks (`coding-kraft-necessity`)

Drawing refined. Tinted descendant subtrees, associated codeword rings with their blocks, added edge bits and the general 2^(m-l) descendant count.

#### Kraft inequality: packing prescribed lengths (`coding-kraft-sufficiency`)

Drawing refined. Exposed the dyadic construction with tinted subtrees, edge bits, and a left-to-right largest-block-first packing arrow; retained the unused leaf and strict inequality.

Retained after contextual and visual checks:

- **Entropy on the probability simplex** (`coding-entropy-extrema`). Binary entropy extrema and the numerical ternary entropy contours are consistent with the entropy lemma; the original ambiguity remains explicitly qualified.


### Wavelets

#### Orthogonalizing a cardinal spline (`wavelets-spline-orthogonalization`)

Drawing refined. Rechecked the linear-spline Gramian `A(omega)=(2+cos omega)/3` and the orthogonalization multiplier. The existing numerical profile is consistent with the Fourier coefficients of `A^(-1/2)` and was preserved exactly. Added axis ticks to distinguish the compact generator's interval `[-1,1]` from the displayed orthogonalized interval `[-4,4]`, which uses a compressed horizontal scale. This prevents the narrow-looking central peak from implying a narrower support. Raised the function labels to avoid the orthogonalized peak touching its label, and separated the support descriptions from the new ticks. The four degree examples and all numerical coordinates are unchanged. Manifest notes explicitly state the different horizontal scales.

#### Sharp and smooth frequency partitions (`wavelets-shannon-spline-bands`)

Drawing refined. Checked the chapter dilation convention and the exact Shannon bands. The scale-minus-one band already had the necessary factor two. Added the identity `2^(-j)|psi-hat_(j,0)(omega)|²=|psi-hat(2^j omega)|²` immediately before the all-scale energy partition, explaining why energies at different scales can be drawn with equal height. Retained the exact sharp bands, the qualified schematic smooth-overlap row, and the restriction of the partition statement to almost every nonzero frequency. The new formula is separated from the axis and partition formula, with no curve changes.

#### Wrapping a function around the unit interval (`wavelets-periodization`)

Drawing refined. The previous figure used a Gaussian and a truncated three-term sum while writing the infinite periodization identity. The omitted terms were visually negligible, but the illustrated equality was not literally exact. Replaced the illustrative function by a compact cosine-squared bump centered at `0.87` with support radius `0.25`. Only the original and its left-shifted copy contribute on `[0,1]`, so the visible sum is now the full periodization. Matching endpoint markers have exactly the same value. Retained the original visual lesson of a right-side bump wrapping across the left endpoint. Manifest notes specify the analytic example and continue to avoid claiming that it is a particular orthogonal wavelet.

#### Vanishing moments and flat filter zeros (`wavelets-vanishing-moments`)

Drawing refined. Checked `eq-cond-vm-fourier` and its proof under the stated QMF and differentiability assumptions. The existing magnitude curves correspond to a two-moment filter pair: their zeros have order two, whereas the formulas state the general p-moment equivalence. Added a short sentence distinguishing the drawn magnitudes from the complex filter symbols whose derivatives appear in the formulas. This prevents interpreting a derivative condition as a statement only about magnitude, which would omit phase information. Neither curve nor the zero/moment formulas changed. The footer fits with clear separation from the equivalence display.

Retained after contextual and visual checks:

- **One representative of each periodic translate** (`wavelets-periodic-translates`). Checked the chapter period-one convention and representative indices `0<=n<2^(-j)`.
- **Haar filter coefficients from overlap integrals** (`wavelets-haar-inner-products`). Checked `eq-dfn-h-g` and the Haar sign convention.
- **Complementary low-pass and high-pass energies** (`wavelets-filter-constraints`). Checked the orthogonality and QMF constraints.


### Approximation

#### Approximation rates across signal models (`approximation-rates`)

Drawing refined. Restored the original Sobolev interpretation in both the sketch headings and rate table, and identified the class as `H^alpha` in the explanatory text. Generic “smooth” labels could conflate Sobolev regularity with noninteger Holder regularity, whose endpoint distinction the chapter explicitly discusses. Retained squared-L2 rate conventions, the condition `alpha >= 1` for the displayed piecewise-smooth rates, the wavelet regularity/moment/boundary requirements, and the uniform bound on the BV class. Verified the rate entries: `M^(-2 alpha)` in one dimension, `M^(-alpha)` in two; `M^-1` for the linear/jump-limited one-dimensional cases; `M^-1/2` for the displayed linear BV/cartoon cases and `M^-1` for nonlinear wavelets; `M^-2` for adaptive cartoon triangulations and `M^-2 (log M)^3` for curvelets. These are upper rates under the stated assumptions. Anchored the explanatory paragraph below the actual table boundary and protected short formulas from stretched word spacing. The smooth analytic sketch and smooth closed cartoon contour remain.


### Inverse problems

#### Capped inversion and Tikhonov filtering (`inverse-problems-variance-bound`)

Drawing refined. Kept the explicitly qualified interpretation of the original plateau as the capped inverse `1/max(sigma,tau)`, rather than truncated SVD. Set both drawings to common axis scales and explicitly imposed `tau=sqrt(lambda)>0`. The capped gain is `1/tau`, whereas Tikhonov's maximum is `1/(2 sqrt(lambda))`; their rendered maxima now have the correct factor of two. Both dashed reciprocal curves also have the same scale. This makes their visual comparison mathematically meaningful without changing the chapter's filter formulas.

#### The source-order dependence of the bias bound (`inverse-problems-bias-bound`)

Drawing refined. Verified the corrected maximum at `sigma*=sqrt(beta lambda/(2-beta))` for `0<beta<2`, the limiting supremum `lambda` for `beta=2`, and the divergence for `beta>2` over the entire nonnegative half-line. Added the assumption `lambda>0` and a direct distinction between this uniform half-line supremum and the bounded spectrum of a fixed bounded operator. Thus the last panel cannot be read as claiming that the bias of every fixed bounded operator is infinite. Protected inline spectral relations from stretched spacing.


### Sparse regularization

#### Soft thresholding: 0<lambda<y (`sparse-regularization-soft-objective-2`)

Drawing refined. Restored the one-sided supporting tangent segments present in the original sketch, with slopes `-y-lambda` and `lambda-y` at zero. Added the complete subdifferential `[-y-lambda,lambda-y]` and an unobtrusive tangent legend. For the plotted `y=2`, `lambda=0.7`, the minimizer remains `1.3` and its value `1.155`; both tangent lines are supporting lines below the convex objective. The negative right derivative explains why the minimizer is positive.

#### Soft thresholding: lambda=y (`sparse-regularization-soft-objective-3`)

Drawing refined. Added the same exact one-sided tangent geometry and full subdifferential in the threshold case `lambda=y=2`. The right tangent is horizontal, the left slope is `-4`, and zero is the minimizer with objective value `2`. This now visibly explains the boundary of the zero-minimizer regime.

#### Soft thresholding: lambda>y (`sparse-regularization-soft-objective-4`)

Drawing refined. Added the exact one-sided tangent geometry and full subdifferential for `lambda=2.7>y=2`. The left/right slopes are `-4.7` and `0.7`; the interval contains zero in its interior, and the minimizer remains zero with value `2`. The tangent segments and text do not obscure the objective.

Retained after contextual and visual checks:

- **Hard thresholding as a scalar minimization** (`sparse-regularization-hard-objective`). Checked the caption normalization `F(x)=(x-y)^2+T^2 1_{x!=0}`.
- **Soft thresholding: lambda=0** (`sparse-regularization-soft-objective-1`). Checked `lambda=0`, `y=2`: the graph is the smooth parabola `F(x)=(x-2)^2/2` with minimizer `2` and minimum zero.


### Sparse recovery theory

#### Interior points and failure of exact recovery (`sparse-theory-polytope-proof`)

Drawing refined. Retained the onto-map hypothesis, `r=||x0||1>0`, and normalized point `q=Ax0/r`. Constructed the extended point exactly on the ray through `q=(0.85,-0.72)` and on the polygon edge from `(2,0)` to `(1.2,-1.65)`: its multiplier is `2/(0.85+0.8*0.72/1.65)`. Replaced the second panel's double-headed marker by the directed displacement `h`, marked `q+h`, and moved that label away from the polygon edge. Broke the proof formulas into readable lines. Both panels concern failure of optimality through interior membership; they do not imply uniqueness of a boundary representation.

#### Removing a dependent active coefficient (`sparse-theory-injectivity`)

Drawing refined. Explicitly marked both finite endpoints of the sign-preserving interval: the first coefficient vanishes at `t_-=-1` and the third at `t_+=4` for the displayed trajectories `x=(0.5,3,1)`, `h=(0.5,-0.25,-0.25)`. Added the minimizer hypothesis: `Ah=0` alone preserves fidelity but does not imply a constant l1 norm. At a minimizer, the affine norm has zero slope while signs are fixed; continuity gives constancy at the endpoints too. Distinguished fixed signs in the open interval from endpoint zeros, and moved a trajectory label off its line.

#### Parameter region for sign consistency (`sparse-theory-recovery-region`)

Drawing refined. Verified the two strict conditions `delta+lambda<T/K` and `R delta<S lambda`, and their intersection at `lambda_max=TR/[K(R+S)]`. The plotted values `T/K=3.91` and `S/R=0.7` give the exact intersection `(2.3,1.61)`. Replaced the filled intersection point by an open marker and clarified that the oblique boundaries are excluded. Explicitly retained `delta>=0`: the zero-noise axis is admissible at positive lambda satisfying the displayed strict conditions. The smaller triangle also has the stated strict cutoff.

#### Convolution of a discrete measure (`sparse-theory-convolution-spikes`)

Drawing refined. Corrected a real scale inconsistency: the individual Gaussian kernels had peak `1.8`, while the output used unit-peak Gaussians with the displayed weights. The kernels now have unit peaks, so the spike heights `1.9` and `1.15` and the output `1.9 K(.-z_j)+1.15 K(.-z_i)` agree exactly on a common vertical scale. Aligned the spike positions with the translated-kernel centers at `0.8` and `2.35`, and added corresponding output-location labels. Preserved the distinction between coefficient indices and spatial grid locations, and the equality between the measure operator and `Ax`.

Retained after contextual and visual checks:

- **Bregman divergence as a vertical gap** (`sparse-theory-bregman`). Verified the strictly convex quadratic and its tangent: at the displayed base point the value and slope agree, and the indicated vertical gap is exactly the Bregman divergence.
- **A strict subgradient margin controls the error** (`sparse-theory-bregman-l1`). Checked the required base point `x0=0` and the strict subgradient margin `|eta|<1`.


### Machine learning and PCA

#### PCA directions and projected variance (`ml-pca-variance`)

Drawing refined. Highlighted one actual observation and its orthogonal projection to clarify which coordinates enter the empirical variance. Computed ellipse semiaxes directly from the displayed 36-point cloud, and separated eigenvector labels from the projection marker.

#### Nearest-centroid assignment and Voronoi cells (`ml-voronoi`)

Drawing refined. Added the centroid segment and its perpendicular-bisector right-angle marker, making the construction of a shared boundary explicit.

#### Conditional expectation as the mean of a slice (`ml-conditional-expectation`)

Drawing refined. Replaced the symmetric-noise example by a centered asymmetric Gaussian mixture. The conditional mean is now visibly different from the mode. Color displays the actual joint density for a uniform input, an arrow explains conditioning/normalization, and the normalized slice uses the same noise law and mean function.

#### Four nearest neighbors of a query point (`ml-nearest-neighbors`)

Drawing refined. Shaded the selected neighborhood and labeled its radius as the exact fourth-neighbor distance. Kept distance ordering and all excluded observations unchanged.

#### Logistic probability and transition width (`ml-logistic-one-dimension`)

Drawing refined. Replaced a single generically labeled width arrow with two correctly scaled measurements, each tied to the actual 0.1 and 0.9 crossings for its slope.

#### The same decision boundary at two parameter norms (`ml-logistic-two-dimensions`)

Drawing refined. Labeled the probability contours and added normal dimension arrows, giving the exact width 2 log(9)/t. Moved contour labels away from the clipping frame.

#### The polytope bounding the PCA objective (`ml-pca-linear-program`)

Drawing refined. Changed the 3D projection so the cutting face has visible area (1.4 square cm rather than .02). Distinguished hidden and removed edges, marked the excluded cube vertex, removed an unexplained 2D guide, and exposed feasibility and maximizing vertices. Moved the beta2 label away from the cut-face equation.

Retained after contextual and visual checks:

- **A joint probability model for input and output** (`ml-joint-model`). Gaussian contour geometry and the two coordinate projections illustrate a joint observation clearly; Gaussianity remains an illustrative choice.
- **An affine fit and a nonlinear regression function** (`ml-linear-fit`). The displayed line is the least-squares affine fit to the displayed points, and residuals are vertical; the nonlinear conditional mean is clearly distinguished.
- **Rows are observations; columns are features** (`ml-pca-data-matrix`). Rows, columns, centering, covariance normalization and SVD dimensions agree with the proof.
- **Linear compression followed by reconstruction** (`ml-pca-compression-reconstruction`). Matrix shapes, rank bound and composition order correctly distinguish arbitrary linear maps from the optimal orthogonal projection.
- **Row norms of the rotated basis** (`ml-pca-row-norms`). B has p rows and k orthonormal columns; row vectors and their squared norms correctly produce total mass k.
- **Completing the columns to an orthogonal matrix** (`ml-pca-orthogonal-extension`). The full orthogonal completion gives unit full-row norms and the beta_i <= 1 bound with correct dimensions.


### Smooth optimization

#### Linear regression (`smooth-linear-regression`)

Drawing refined. Replaced the schematic slope 0.45 with the exact least-squares slope of the fifteen displayed observations, approximately 0.4589391976. The line is now labeled by the minimizing parameter rather than an unspecified parameter. Its slope is `sum(a_i y_i)/sum(a_i²)`; the normal-equation residual is approximately −1.3×10⁻¹⁵. Marked and labeled the vertical residual for the selected observation in red, matching the chapter's convention `y_i−〈a_i,x〉`. Retained a homogeneous scalar-feature model with no intercept and stated this interpretation in the notes. Moved the line label above the cloud to keep the right margin clear.

#### A linear classifier (`smooth-linear-classifier`)

Drawing refined. Added the feature-space origin and a right-angle marker between the parameter vector and the separating hyperplane. These annotations make the homogeneous classifier's geometry explicit. Checked that each positive observation has positive score and each negative observation has negative score. Preserved the boundary and normal directions; their dot product is zero to the precision of the stored coordinates. Positioned the origin label below and to the right so that the separating line does not cross the glyph.

#### Existence and uniqueness of minimizers (`smooth-minimizer-examples`)

Drawing refined. Distinguished the two nonattainment mechanisms directly in the drawing: `inf f=−∞` for the downward quadratic and `inf f=0`, not attained, for the exponential. Enlarged the exponential's displayed profile, preserving a positive gap above its horizontal asymptote. Added continuation arrows to the unbounded profiles. Retained the two-point, interval, and singleton minimizer examples and their exact highlighted minimizer sets. Increased the lower canvas margin to accommodate the explanatory labels.

#### Least squares with a nontrivial kernel (`smooth-least-squares-flat`)

Drawing refined. Added a second arrowhead to the highlighted affine minimizer line. The drawing now communicates unboundedness in both kernel directions, as required for `x-star+ker(A)`, rather than resembling a single ray. Retained the trough geometry and the nontrivial-kernel hypothesis.

#### Failure of the secant inequality (`smooth-secant-nonconvex`)

Drawing refined. Added the explicit reversed secant inequality below the picture, together with the definition of the interior point. Verified that the graph and interpolated endpoint value are compared at the same abscissa: the marked graph value is 3 and the secant value is approximately 1.31863. Preserved the original nonconvex geometry. Expanded the lower canvas after rendering exposed a cropped second equation line; the final equation is complete.

#### Convexity and the secant inequality (`smooth-secant-convex`)

Drawing refined. Connected the interpolated endpoint-value label to its red point with a fine leader and placed the label above the graph. Added the convexity inequality and interior-point definition below the axes. Verified the endpoint values, the interior graph value, and the supporting tangent against the representative quadratic. Expanded the lower canvas so both formula lines render completely.

#### Strict convexity (`smooth-strict-convexity`)

Drawing refined. Added a marked interior midpoint on both the graph and secant, joined by a red double arrow showing the positive strictness gap. The graph value is 0.6056, the secant value is 1.4806, and the gap is exactly 0.875. Stated the strict inequality and `0<t<1` explicitly. The supporting tangent, endpoint values, and distinct endpoint positions were retained.

#### Convexity with an affine interval (`smooth-nonstrict-convexity`)

Drawing refined. Added the exact secant-equality identity on the affine interval and defined the interior point. Verified that the two surrounding Bézier pieces are convex and that their slopes join the straight interval without a downward jump. In a cubic parameter, their curvature numerators are respectively `4.464−4.032t−0.108t²` and `19.44+10.08t+62.64t²`, both positive for `0≤t≤1`. Thus the diagram really is convex, while its affine interval prevents strict convexity.

#### The gradient is a vector field (`smooth-gradient-field`)

Drawing refined. Made arrow length agree with the displayed identity `gradient f(x)=x`. Each arrow now starts at a marked point q and ends at 2q, so its displacement is exactly q in the same coordinates. Previously the arrows represented one-half of their base vectors without stating the scale. Added small base-point markers and retained the zero gradient at the minimizer. The accompanying projected bowl remains an uncalibrated illustration of the same radial objective.

#### A stationary point need not be an extremum (`smooth-stationary-inflection`)

Drawing retained; interpretation note corrected. Verified the cubic, its horizontal tangent at the origin, and the precise statement that the stationary inflection point is neither a minimum nor a maximum. Removed the stale manifest claim about an old chapter caption, because the current chapter already uses the corrected terminology. The source drawing and PDF appearance were retained.

#### PCA and least-squares level sets (`smooth-pca-quadratic-geometry`)

Drawing refined. Moved all width annotations out of the ellipses and into separate lines above the two panels. This removes the white label boxes that interrupted several contours and hid the geometry. Preserved all twelve data points, the common eigenvectors, every ellipse, and the exact covariance construction. The empirical mean is zero and the covariance in principal coordinates is `diag(1,0.105625)`. The covariance axis ratio is 0.325; the quadratic level-set ratio is its reciprocal. Retained the explicit common scale factors in the bottom caption. Increased title clearance so the new width lines and the eigenvector labels do not touch.

#### A first-order Taylor approximation (`smooth-taylor-line`)

Drawing refined. Corrected an unused named tangent-point coordinate from 2.6264 to the actual value 2.5466, then used that named point for both the gap arrow and point marker. Added a label for the tangent value at z, complementing the label for the true function value. The exact values are `f(x)=1.4598`, `f(z)=3.1888`, and `T_x(z)=2.5466`; the affine line has derivative 0.418 at x. The remainder is therefore represented by the correct vertical gap.

#### Gradients and nested level sets (`smooth-nested-level-sets`)

Drawing refined. Replaced labels that interrupted the three ellipses with a separate key using the contours' three intensities. All level curves are now complete and the ordering of their values remains explicit. Added the qualification that arrow lengths are schematic; the arrows depict normal directions, which are generally not radial for an ellipse. Removed the redundant teal arrow underneath the highlighted red gradient arrow. Retained the mathematically correct outward normal directions.

#### Exact line search (`smooth-exact-line-search`)

Drawing refined. Added the actual quadratic objective to the drawing, so the level sets and numerical trajectory have an explicit common interpretation. Connected the small final steps that had previously appeared as isolated points. Recomputed the trajectory from `(3,1.1)` using the exact quadratic line-search formula; the stored coordinates agree to their displayed precision. Consecutive gradient inner products vanish to numerical precision, with a maximum residual below 9×10⁻¹⁶. Preserved the right-angle marker and the exact-search formula.

#### The optimal quadratic step size (`smooth-quadratic-contraction`)

Drawing refined. Stated that μ and L are the actual smallest and largest eigenvalues, making the displayed equality for h valid rather than merely an upper bound. Stated the representative ratio `L=3μ` that determines the symbolic drawing. Moved the optimal contraction-value label into open space, removing the white box that had interrupted the smaller dashed branch. Verified the maximum of the two absolute-value branches, the optimal step `2/(L+μ)`, the minimum `(L−μ)/(L+μ)`, and the strict contraction interval `0<τ<2/L`.

Retained after contextual and visual checks:

- **Convexity of functions and sets** (`smooth-convex-examples`). Checked the three function graphs, their epigraphs, and the corresponding set diagrams against the convexity definitions.
- **Classification losses** (`smooth-classification-losses`). Checked the negative-margin variable, the unscaled logistic loss, the hinge breakpoint at −1, and the 0–1 jump at zero.
- **Coercive least squares** (`smooth-least-squares-coercive`). Verified the full-column-rank hypothesis, the unique quadratic minimizer, and the projection of the minimizer onto the domain.
- **Stationarity at local extrema** (`smooth-stationary-extrema`). Verified the two minima and central maximum of the representative double well, their horizontal tangents, and the statement that the derivative vanishes at every marked point.
- **Stationarity for a convex objective** (`smooth-stationary-convex`). Checked that the projected horizontal plane passes through the bowl's minimum and supports the graph from below.
- **A tangent plane** (`smooth-taylor-plane`). Verified the plane against the representative paraboloid: at the selected domain point `(−1,0)`, the surface height is 0.9 and the tangent has equation `z=0.1−0.8u`.
- **The gradient is orthogonal to a level curve** (`smooth-level-set-tangent`). Checked the level-curve membership, the marked curve point, and the perpendicular tangent and gradient at the leftmost point of the ellipse.
- **Step-size effects in gradient descent** (`smooth-gradient-step-size`). Rechecked the two exact trajectories for the displayed quadratic with Hessian `diag(1,5)`.
- **Quadratic upper and lower bounds** (`smooth-curvature-bounds`). Verified the three representative quadratics: the lower curvature is 0.11, the objective curvature is 0.21, and the upper curvature is 0.30.


### Stochastic optimization

#### An unbiased stochastic gradient (`advanced-unbiased-gradient`)

Drawing refined. Added the uniform sampling distribution required for the arithmetic-mean objective. Lightly shaded the triangle joining the three individual gradient endpoints; the mean endpoint is its centroid. The vectors are `(−2,4)`, `(3,3)`, and `(1,2)`, whose mean is `(2/3,3)`. Preserved the exact mean arrow, moved its label away from the triangle boundary, and retained the expectation identity. The new geometry makes unbiasedness visible as an arithmetic average.

#### Stochastic gradient trajectories (`advanced-sgd-trajectory`)

Drawing refined. Replaced the phrase implying that convergent step sizes suffice with a direct description of an illustrative trajectory and shrinking neighborhoods. The assumptions on the objective and noise are also needed for convergence. Replaced the final solid segment by a dotted continuation so that a displayed finite iterate is not identified with the limit minimizer. The ellipses remain qualitative neighborhoods, not confidence regions, and the notes retain that limitation.


### Deep learning

#### A fully connected feedforward network (`deep-learning-fully-connected`)

Drawing refined. Retained the three inputs, two hidden layers of four units, three outputs, and all `12+16+12=40` adjacent-layer connections. Replaced generic input/output component labels by `(x0)_i` and `(x3)_i`, and identified the concrete layer vectors `x0,x1,x2,x3`. The final vector is explicitly a score vector, avoiding confusion with the chapter's training target `y`. The general layer recurrence is preserved. All labels and connections remain legible; crossings are not rendered as nodes.

#### Feature maps in a convolutional network (`deep-learning-convolutional`)

Drawing refined. Corrected activation notation in the two feature-producing arrows: the chapter's `rho_l` already includes downsampling, which this diagram shows as a separate operation. The arrows now use the pointwise activations `tilde rho_0` and `tilde rho_1`, followed by the separately shown downsampling steps. Named the output score vector `x_L`, avoiding confusion with a target `y`. Verified `n_l=bar n_l d_l`, RGB input `d0=3`, and channel counts preserved through each downsampling. The original's unreadable exact spatial/filter dimensions remain explicitly qualified in the manifest rather than invented.


### Convex analysis

Retained after contextual and visual checks:

- **Convexity of functions and sets** (`convex-analysis-convex-examples`). Checked all six panels.
- **Supporting affine functions at a kink** (`convex-analysis-subdifferential`). Checked the displayed convex kink `1.1+|z|+0.12 z^2`: all five drawn slopes from `-1` to `1` are global supporting slopes at zero.
- **Graphs of two subdifferentials** (`convex-analysis-subdiff-l1`). Verified the full closed interval `[-1,1]` at the kink of absolute value.
- **Normal cones at boundary, interior, and exterior points** (`convex-analysis-normal-cones`). Checked the polygon geometry and the outward directions at each selected boundary point, including the negative-quadrant cone at `(0,0)` and the cone bounded by angles `-45` and `45` degrees at `(4.5,1.5)`.


### Nonsmooth optimization

#### Projection and proximal maps (`optim-nonsmooth-prox-projection`)

Drawing refined. Verified every projection pair for the quarter disk: axis projections, the origin, endpoints of the circular arc, and the diagonal arc projection all satisfy the corresponding normal-cone condition. Named the right-hand function explicitly as `f(u,v)=(u^2+2v^2)/2`. For the exact proximal pair `x=(1.4,1.35)` and `z=(0.7,0.45)`, one has `x-z=(0.7,0.9)=grad f(z)`. Replaced separately rounded ellipse radii by radii computed from `f(z)`, so the point lies exactly on the highlighted level set. Shortened the title and retained the essential qualification that the sublevel set depends on the input.

### Verification of the second TikZ pass

- Rebuilt the complete 341-page book and all 19 standalone chapter PDFs. All 20 final documents have zero compilation diagnostics, including warnings, missing characters, and overfull or underfull boxes.
- Preserved all 91 original/TikZ comparisons in the book and the same 91 comparisons across the standalone chapters. The original scanned figures remain beside the reconstructions for review.
- Checked all 91 individual figure PDFs and all 91 integrated book comparison pages visually. Compared the standalone review pages with their book counterparts: 59 render identically above the footer; the other 32 differ by the expected 3 mm mirrored-margin offset and were separately checked visually. No clipping, overlapping labels, or truncated explanatory notes were found in the final comparison pages.
- Verified that all 91 reconstruction PDFs are current, single-page vector figures with embedded fonts, no raster image objects, and no text extending beyond the page boundaries.
- Checked 1,839 internal PDF links, 150 chapter-to-book links, and 725 shared labels. All link destinations and shared equation/figure numbers pass. Adjusted the PCA context sentence so its reference resolves correctly to equation (13.6).
- Independently checked representative quantitative constructions: aliasing mass, PCA covariance and semiaxes, the asymmetric conditional-density mixture, nearest-neighbor distances, FFT padding grids, spherical projection, least-squares regression, convexity gaps, gradient vectors, and proximal geometry. These checks complement the contextual review of every reconstruction; figures that are schematic remain explicitly qualified as such.

Detailed review records, numerical checks, page renders, and final build reports are retained locally under `build/tikz-pass2/`. The PDF and reference audits are recorded in `build/pdf-audit.json`, `build/tikz-pdf-audit.json`, and `build/label-number-audit.json`.

## Separate figure comparison volume

Moved the review material out of the reading editions and into `figure-processing/figure-comparisons.pdf`. The complete book and all standalone chapters now end normally, without trailing original/TikZ comparison sections or review entries in the table of contents.

- Preserved all 91 side-by-side comparisons, their interpretation notes, and their stable drawing identifiers. The separate volume has a cover and 91 comparison pages, ordered by the main book's figure numbers.
- Matched every original asset to its enclosing figure in the chapter source. The comparisons cover 69 numbered book figures; several figures contain multiple separately reconstructed panels. Each comparison prints the exact book figure number and full caption, followed by an explicit panel identifier when needed. The reconstruction's descriptive title remains separate from the book caption.
- Added 25 missing reference labels to existing numbered captions without changing their wording. Added captions and numbers to four previously unnumbered sketches, so every reconstruction has a precise book reference:

| Figure | Caption added to the main book |
| --- | --- |
| 13.4 | Rows of the centered data matrix are observations; columns are features. |
| 13.5 | Linear compression followed by reconstruction. |
| 14.6 | The gradient as a vector field. |
| 14.12 | Quadratic upper and lower bounds for a smooth strongly convex function. |

- Added links from comparison headings and page references to the corresponding figures in `../FundationsDataScience.pdf`. Added 14 chapter bookmarks and 91 individual comparison bookmarks.
- Added `figure-processing/figure-index.json`, recording each stable identifier, current figure number, complete caption in LaTeX, book page, comparison PDF page, panel, and source asset. The manifest fields `book_label` and `book_panel` maintain the association with the original figure.
- Updated the normal build to compile the book first, derive comparison captions and numbers from its stabilized references, then publish the separate volume alongside the book and chapter PDFs. The build rejects missing figure labels and continues to reject any final compilation diagnostics. It uses reference data only from the active chapters.
- Added a dedicated landscape layout that reserves space for the caption and notes before sizing the image panels. Kept the book's fonts and colors, added authorship, affiliation, and date to the comparison cover, and fixed a font-shape substitution when importing the upright definition marker in a mathematical caption.
- Removed the former shared `figure-reviews.tex` integration and documented the new output and reference workflow in the root README, the TikZ README, and `figure-processing/README.md`.

### Verification

The rebuilt main book has 250 pages; the separate comparison volume has 92. All 21 output PDFs (book, 19 chapters, and comparisons) compile with zero warnings, missing characters, or overfull/underfull boxes. Verified the absence of review pages in all 20 reading PDFs, the source-to-caption mapping for all 91 comparisons, and all 754 shared labels. Checked 1,791 internal links and 150 chapter-to-book links in the reading PDFs, plus 203 comparison-to-book links. All destinations resolve.

Visually checked the comparison cover and every comparison page, as well as all pages containing the four newly captioned sketches in the main book and standalone chapters. The comparison PDF has embedded fonts, valid bookmarks, and no text outside page boundaries. Validation records and page renders are retained in `build/figure-processing/`; the main PDF and label audits are in `build/pdf-audit.json` and `build/label-number-audit.json`.


## Publication of the TikZ figures and requested corrections

Replaced all 91 reviewed illustrations with their editable TikZ reconstructions in the main book and standalone chapters. Also replaced four repeated occurrences of the scalar soft-thresholding drawings in the advanced optimization chapter. This updates 95 image occurrences across 70 numbered figures; the separate comparison volume retains all 91 original/new pairs and their interpretation notes. Original assets remain available on disk.

### Figure 1.4: sampling with and without aliasing

Rebuilt the diagram around the four stages in the supplied `figure-processing/shannon-interp/` reference images: original signal, sampling, interpolation kernel, and reconstruction. The two cases appear side by side, each with matched frequency and time views. With unit sample spacing, the examples are `sinc(t/4)^4` (bandlimited to the Nyquist interval) and `sinc(t/2)^4` (overlapping spectral replicas). The signal/spectrum pairs and the aliased reconstruction are computed from exact formulas rather than drawn independently. The reconstruction passes through the sample points in both cases, while only the bandlimited case recovers the original signal. Added consistent axes, the cardinal interpolation formula, spectral-replica guides, and a legend distinguishing the folded spectral mass and the dashed original signal. Adjusted heading and formula spacing after rendering.

### Figure 2.3: convolution and the marked output value

Added a dashed curve representing the convolution with the displayed box kernel and a marked point at `(x,(f*g)(x))`. The curve is the exact moving integral of the illustrative signal; the shaded integral and marked value agree. Kept the chapter's unnormalized box convention: its mass is `2 epsilon`, so the output is an integral, not a unit-mass average. Updated the caption and reconstruction notes accordingly.

### Figure 2.10: FFT connections and interleaving

Separated the two input halves and drew all four connections to the sum and difference blocks. Crossing wires have a small visual bridge so that crossings cannot be mistaken for junctions. Preserved the twiddle multiplication before the lower half-size Fourier transform. Added six explicit arrows showing how the first three entries of the even- and odd-frequency streams alternate in the output. Corrected the accompanying interleaving definition to use zero-based indices consistently with the DFT. Moved the lower formulas to clear the interleaving legend.

### Figure 3.4: entropy contours and the probability constraint

Added the requested left panel in the positive `(p1,p2)` quadrant. It plots actual level sets of `H2(p1,p2)=-p1 log2(p1)-p2 log2(p2)`, the line `p1+p2=1`, and the constrained maximizer `(1/2,1/2)`. The level `H2=1` is tangent to the constraint there. The unconstrained entropy extension has its maximum at `(1/e,1/e)`; the contours correctly reflect this different center. Numeric labels lie on their corresponding contours. Retained the middle binary-entropy restriction and the right three-symbol simplex, and expanded the caption to explain all three views.

### Figure 14.1: normalized logistic upper bound

Divided the logistic curve by `log(2)`, so it meets the binary loss at zero margin with value one and upper bounds it everywhere. Updated the legend, vertical scale, and boundary marker. Applied the same normalization to the chapter's loss definition and derivative, which is the sigmoid divided by `log(2)`. Clarified the convention that zero margin counts as an error and that the bound is tight at that boundary. Retained the hinge loss and the open/closed endpoint markers for the binary loss.

### Layout and references

- Use the vector drawings at their natural size, shrinking only those wider than the available text area. This preserves their lettering instead of inheriting the scans' small layout dimensions.
- Stack separately reconstructed components with explicit panel headings. Four long sequences use continued captions under the same figure number: the scalar penalties in Chapters 10 and 15, the convexity examples, and the Taylor/level-set illustrations. The original figure numbers are preserved.
- Updated directional references and captions where panels moved from horizontal arrangements to stacked layouts. Comparison links point to the page containing the corresponding panel, including continued figures.
- Removed blank paragraph breaks from two captions when transferring them into standard figure environments; this fixes a caption-anchor compilation error.
- Marked every reconstruction as published in the manifests and generated figure index. Updated the root, TikZ, and comparison-volume documentation to describe the published figures and the separate archive of originals.


### Verification of the published replacements

- Rebuilt the 272-page main book, all 19 standalone chapters, and the 92-page comparison PDF. All 21 outputs have zero compilation warnings, missing characters, or overfull/underfull boxes.
- Verified from the actual compilation inputs that every one of the 91 TikZ drawings appears in the book and its owning chapter, and that the four repeated soft-thresholding drawings are also replaced. None of their original assets is included in the reading editions; every original and reconstruction is included in the comparison volume.
- Checked all 91 comparison numbers, complete captions, panel identifiers, and links against the main book. Verified 758 shared labels, 1,791 internal links and 150 chapter-to-book links in the reading PDFs, and 203 comparison-to-book links. All destinations resolve.
- Visually checked all 65 book pages containing replacements, all 91 comparison pages, and 12 standalone pages covering the five requested corrections and the larger panel arrangements. No visible clipping, overlapping labels, or truncated captions remain. An extraction-only bounds flag on the entropy page comes from hidden, clipped text in pre-existing plot assets; the rendered page is clean.
- Confirmed that all 91 individual reconstruction PDFs are single-page vector drawings with embedded fonts, no raster image objects, and text within their page boundaries. Fonts are also embedded in the reading and comparison PDFs.
- Numerically checked the underlying Shannon interpolation formulas (integer-sample error below `1e-15`, comparison with a long cardinal sum below `5e-15`), entropy contour construction (residual below `1e-15`), the convolution against independent quadrature (error below `2e-12`), and the FFT wiring against direct complex DFTs of lengths 8, 16, and 32 (error below `2e-13`). Checked the normalized logistic loss at the boundary and across positive and negative margins, and its derivative against finite differences.

Detailed publication checks, numerical records, and page renders are retained in `build/figure-feedback/`. The separate-volume checks are in `build/figure-processing/`; the final build, PDF-link, and shared-label reports remain in `build/build-report.json`, `build/pdf-audit.json`, and `build/label-number-audit.json`.
