# Major corrections to the English and French editions

Revision of 9 September 2026. This log records substantive mathematical repairs,
important changes of interpretation, and changes to the organization and rendering.
It does not enumerate every corrected spelling, sentence or punctuation mark.
Unless stated otherwise, the changes were made in both languages.

The original files are preserved unchanged in
`legacy/pre-revision-2026-09-08/`, with a SHA-256 manifest. Chapter references below
refer to the revised books, which now contain the same five chapters.

## 1. Information theory and compression

- Corrected the logarithm example: `log₂(1/4) = −2`. The entropy of a uniform
  four-symbol source is **2 bits**, not 1. Corrected the deterministic-source
  example so that its message agrees with the symbol assigned probability one.
- Recomputed the finite coding example. The 25-symbol message has counts
  `(2, 9, 12, 2)`; its fixed-length representation has 50 bits and its displayed
  prefix representation has 42 bits. Its empirical entropy is approximately
  `1.621901` bits per symbol, and the mean prefix length is `42/25 = 1.68`.
  Both complete bit strings were checked against the coding tables.
- Corrected the fixed-length formula to `ceil(log₂ N)` and distinguished a
  codeword from its source symbol. The conventions for nonempty codewords, known
  message length and a one-symbol effective alphabet are now explicit.
- Distinguished distinct codewords, unique decodability and the prefix condition.
  Added the uniquely decodable, non-prefix example `{0, 01}`. It needs lookahead;
  backtracking is a possible cost of other non-prefix codes, not a necessary
  operation for this particular example.
- Corrected the relation between tree depth, word length and number of levels:
  words of length `ℓ` lie at depth `ℓ`, and the complete tree has `ℓ+1` levels
  when the root is counted. A decoder emits the source symbol at a codeword leaf.
- Added a proof of Kraft's inequality using disjoint sets of descendants at a
  common depth. Expanded the entropy lower-bound proof with
  `Z = Σ 2^(−ℓᵥ) ≤ 1`, `qᵥ = 2^(−ℓᵥ)/Z`, and
  `ℓᵥ = −log₂ qᵥ − log₂ Z`. The resulting gap is
  `L − H = KL(p | q) − log₂ Z ≥ KL(p | q) ≥ 0`; the proof also covers zero
  source probabilities.
- Clarified that Shannon's bound concerns **expected length**, not a lower bound
  on every individual codeword. Distinguished empirical frequencies from a
  source probability model and accounted for model and framing overheads.
- Stated the strict Huffman upper bound `H ≤ L < H+1` with its effective-alphabet
  assumptions. Explained the block rate `H ≤ Lₙ/n < H+1/n` for independent,
  identically distributed symbols and the special case of a deterministic source.
- Corrected difference coding to retain the initial pixel. The displayed
  transformation is bijective and its inverse is a cumulative sum. A bijection
  preserves the **joint entropy** of the message; a gain from separate-symbol
  coding is explained through the marginal model, not an unconditional decrease
  of entropy under differentiation.
- Corrected the binary storage units and payload estimates: 16 KiB, 13.36 KiB
  and 9.28 KiB at the displayed rounded rates for a 256 × 256 image. The first
  pixel and the code-table overhead are separate costs.
- Distinguished source coding, channel coding and quantization. Removed claims
  of universal codec optimality; JPEG 2000 supports both lossless and lossy modes.

## 2. Image processing

- Matched the discussion to the actual **256 × 256** flower image. An 8-bit
  grayscale representation contains 65,536 pixels, 524,288 bits and 64 KiB.
  Distinguished decimal kilobytes from binary kibibytes; retained the 240 × 240
  calculation only as a separately identified example. Corrected the zoom to a
  5 × 5 pixel block.
- Defined four-level quantization by `q = floor(a/64)` and the displayed gray
  value by `85q`. Corrected the dark-gray and light-gray indices and distinguished
  a packed two-bit representation from storing those values in eight-bit bytes.
- Explained that reducing both image dimensions by a factor of four reduces the
  pixel count by sixteen. Added the need for prefiltering before subsampling to
  prevent aliasing.
- Recomputed the filter example: the nine values sum to 1273, giving mean
  `1273/9 ≈ 141.44` and median 166. Specified that all output pixels are computed
  from the same input image and that the boundary convention matters.
- Qualified the denoising interpretation. The variance reduction `σ²/m` for an
  average of `m` samples requires independent, zero-mean noise and a locally
  constant signal. A median is not universally better than an average.
- Corrected the centered-difference edge example to `sqrt(19746) ≈ 140.52`.
  The centered gradient at unit spacing has magnitude `ℓ/2`. The response
  vanishes when the two opposite pairs agree; this does not require all five
  stencil values to be equal, and some oscillations can therefore be missed.
- Corrected additive color mixing: red plus green gives yellow; red plus blue
  gives magenta. Distinguished separate color channels from statistically
  independent ones, and ideal complementary CMY values from physical CMYK printing.
- Replaced the identification of the arithmetic RGB mean with luminance by an
  explicitly named mean-channel grayscale conversion. Explained that perceptual
  or relative luminance uses another definition and depends on the color space
  and transfer function; linked the W3C definition of relative luminance.
- Corrected the power-law interpretation: squaring a normalized nonnegative
  intensity darkens it and taking its square root brightens it. The inverse
  relationship holds before rounding and clipping. Verified the existing gamma
  panels against their labels.
- Qualified the color transformation that retains channel residuals: it preserves
  the prescribed channel mean before clipping; clipping can alter the mean and
  color relations. It is not an exact general hue-preserving operation.
- Corrected the index formulas for transposition and clockwise rotation, then
  checked the published arrays against those operations pixel by pixel. Explained
  that crossfading combines values at fixed positions rather than transporting
  objects between positions.

## 3. Sampling, sparsity and inverse problems

- Stated a precise bandlimited sampling theorem, including the frequency bound,
  sampling interval, normalized sinc convention and convergence assumptions.
  Distinguished the ideal infinite-duration statement from finite, noisy,
  quantized measurements. Corrected its Shannon reference to the 1949 paper.
- Clarified that `‖x‖₀` counts nonzero entries and is not a norm. General sparse
  optimization is difficult in the worst case; that does not make every sparse
  reconstruction instance computationally intractable.
- Repaired the best `M`-term approximation formula using Parseval's identity and
  the missing summation index. Addressed ties and zero coefficients. A fraction
  of retained coefficients is not itself a bitrate: locations, amplitudes,
  quantization and coding overhead also matter.
- Corrected sensing-matrix dimensions and the distinction between explicit and
  implicit operators. Invertibility is required for an inverse formula, and
  ill-conditioning is characterized by singular values, not by arbitrary
  eigenvalue magnitudes.
- Corrected the denoising description to threshold the **observed** signal, not
  the unknown clean signal. Best approximation of the observation need not
  minimize the error to the clean signal.
- Clarified the fixed-vector limit as the exponent tends to zero. The limiting
  bounded feasible sets have truncated coordinate axes; they must not be confused
  with the unbounded set defined solely by a constraint on `‖x‖₀`.
- Separated the constrained Lasso radius `τ`, the penalized Lasso weight `λ`,
  and the noise tolerance `ε` in basis pursuit denoising. Removed the false claim
  that taking `τ → ∞` guarantees exact reconstruction: with an underdetermined
  operator, unregularized data fitting can have many solutions.
- Replaced informal recovery claims with a stable, uniform Gaussian compressed
  sensing statement. The matrix entries have variance `1/P`; under
  `P ≥ C M log(eN/M)`, the recovery error is bounded by
  `C₀ σ_M(x₀)₁ / sqrt(M) + C₁ ε`, with stated probability, sparsity range,
  orthonormal-basis and noise assumptions. Exact sparse recovery uses the
  **noiseless formulation `ε = 0`**, not merely a zero realized noise vector
  with an arbitrary positive tolerance.
- Explained the role of centering and normalization when replacing Gaussian
  measurements by physical zero/one masks. Corrected the single-pixel camera
  description: mirrors direct light toward or away from the detector; they do
  not independently record pixel values. Corrected the compression factor to
  `Q/P = 6`, hence `P/Q ≈ 1/6`.
- Qualified claims about security, energy savings and structured MRI acquisition.
  Keeping a sensing matrix secret is not by itself a cryptographic guarantee;
  physical savings depend on acquisition and reconstruction costs.
- Replaced the Lena-based illustrations, including the pixel zoom hidden under a
  generic filename, by reproducible flower examples. Added corresponding red
  boxes in the sampling illustration and embedded TrueType fonts in its plots.
- Computed actual Daubechies wavelet approximations: `db4`, five levels, periodic
  boundaries, and 16,384, 8,192 or 4,096 retained coefficients. Their relative
  Euclidean errors are approximately 2.12%, 3.80% and 5.81% before display clipping.
- Corrected the inverse-problem comparison to display an **actual reconstruction**.
  The deterministic mask retains 6,554 of 65,536 pixels. The stated sparse
  penalized problem is solved by FISTA with `λ = 0.005`; the supplied run reaches
  a relative proximal-gradient residual below `10⁻⁷` in 1,750 iterations.
  The reconstruction PSNR is approximately 21.41 dB before display clipping.

## 4. Optimal transport

- Made both the one-destination-per-source and one-source-per-destination
  constraints explicit for an assignment of equal unit batches. Distinguished
  the sum of journey times from a traveling-salesperson route or a parallel
  completion time.
- Transcribed and verified the Paris cost matrix. The example assignment has
  cost 65; `(5, 6, 2, 1, 3, 4)` has optimal cost **64**. Added a checkable dual
  certificate with potentials `u = (2, 3, 3, 0, 5, 2)` and
  `v = (10, 4, 9, 7, 8, 11)`. All 36 dual constraints hold and the dual value is
  64. Exhaustive enumeration of all 720 permutations confirms the optimum.
- Corrected the worst-case insertion-sort count for 70 entries to
  `70 × 69 / 2 = 2415` comparisons per list; explained the need to sort both
  lists and the availability of `O(n log n)` comparison sorting.
- Added the uncrossing argument for sorted matching on a line. Explained that
  the optimum for absolute distance need not be unique. Qualified the planar
  noncrossing statement to separate transverse crossings from degenerate
  collinear or endpoint configurations.
- Corrected the algorithmic complexity claim: suitable Hungarian assignment
  algorithms have cubic complexity; the general simplex method does not thereby
  acquire a polynomial worst-case guarantee.
- Explained why the doubly stochastic relaxation has a permutation minimizer
  using Birkhoff's theorem. A non-permutation matrix can also be optimal when
  the optimum is not unique.
- Stated nonnegativity and equal total mass for a general transport plan, gave
  the feasible construction `abᵀ/s` for positive common mass `s`, and treated
  zero mass separately. Explained why unequal source and destination masses may
  require splitting a source. Distinguished mass-`n` doubly stochastic matrices
  from probability couplings normalized by `1/n`.
- Defined Wasserstein distance on measures that include their support positions,
  not just arrays of weights. **`W₂` is a metric; `W₂²` is not.** The Dirac-mass
  example at 0, 1 and 2 demonstrates the failure of the triangle inequality for
  the square.
- Stated the essential hypotheses of the quadratic Brenier theorem: finite
  second moments and an absolutely continuous source. The unique optimal plan
  is induced almost everywhere by the gradient of a convex potential under
  those hypotheses, not for every pair of measures and every cost.
- Corrected the interpretation of palette transfer: RGB-space transport changes
  colors without moving spatial pixels. A permutation preserves the target
  empirical color distribution, whereas a barycentric map from a split coupling
  need not reproduce it exactly.
- Defined Wasserstein barycenters using **squared** distances and nonnegative
  weights summing to one. Normalized shape densities to probability measures and
  removed general assertions of a sharp uniform output shape or unconditional
  uniqueness.
- Clarified normalization of word frequencies, the role of the ground cost and
  the loss of word order. Checked the metro illustration against its assignment:
  the retained higher-resolution line-1 image corresponds to `(3, 2, 1, 5, 4)`;
  the previous English caption did not match its lower-resolution illustration.

## 5. Neural networks

- Integrated the previously separate discriminative and generative articles as
  Chapter 5 in **both complete books**, with an independent chapter PDF in each
  language. Removed the repeated transport introduction now covered in Chapter 4.
- Distinguished annotated class labels, one-hot targets, estimated probabilities,
  accuracy and calibration. Single-label softmax classification is not the same
  task as independent multilabel classification.
- Specified the dimensions of weights and biases in a two-layer network, with
  hidden ReLU activations, linear logits and a softmax output. Added the stable
  softmax computation by subtraction of the largest logit.
- Corrected vector-valued squared loss to a squared Euclidean norm with the
  stated averaging convention. Defined one-hot cross-entropy with the natural
  logarithm.
- Separated **backpropagation**, which computes derivatives, from the optimizer
  update. Added explicit two-layer gradient formulas and checked them against
  centered finite differences. Minibatch gradient steps need not decrease the
  full objective at every step, and nonconvex training has no general guarantee
  of reaching a global minimum.
- Explained the derivative convention at the ReLU kink and qualified the
  biological analogy. Stated the relevant universal-approximation theorem with
  its activation, compact domain and linear-output assumptions; representability
  is distinct from training success, sample complexity and generalization.
- Explained convolution through locality and shared weights. Distinguished
  translation equivariance from invariance and noted boundary and subsampling
  effects. Corrected the history to the **2018 Turing Award, announced in 2019**,
  and corrected Yann LeCun's name.
- Defined a generator through a pushforward distribution, rather than as the
  inverse of a classifier. Distinguished conditional and unconditional models;
  removed the suggestion that generative modeling began with GANs in 2014.
- Normalized empirical distributions, including repeated samples, and clarified
  the finite-batch transport objective. An empirical transport value is not
  generally an unbiased estimate of the population transport distance.
- Distinguished the original GAN objective and its ideal Jensen–Shannon
  interpretation from the Wasserstein GAN. A WGAN critic must satisfy a
  **1-Lipschitz constraint for a fixed ground metric**; an arbitrary second neural
  network does not automatically compute or learn a valid transport distance.
- Qualified convergence and mode-coverage claims. Distinguished a synthetic
  interpolation from a transport geodesic or an impersonation deepfake. Identified
  BigGAN as class-conditional and corrected its ICLR publication year to 2019.
- Restored the original English-labeled network diagrams in the English edition.
  Added the missing WGAN coauthor, Soumith Chintala, and corrected its proceedings
  metadata.

## Presentation, encoding and repository organization

- Rewrote and aligned the English and French exposition, with a short motivation
  and roadmap at the start of every chapter. Repaired awkward translations,
  inconsistent terminology, repeated introductions and malformed references.
- Replaced mixed Latin-1/Windows-1252/UTF-8 source encodings with NFC-normalized
  UTF-8. Fixed corrupted French title text, accents, ligatures and punctuation.
  Used T1 encoding, French Babel conventions, embedded fonts and Unicode mappings
  for PDF text extraction.
- Adopted the modern MFDS/DAFT typography: New PX text and mathematics, Source Sans
  headings, a restrained teal palette, consistent captions and running heads.
  Retained the original wave cover illustration. Refined figure placement and
  glossary layouts. Collected the complete-book glossaries in the back matter,
  keeping the chapter glossaries in the standalone PDFs. Consolidated full-book
  acknowledgements in the preface while
  retaining chapter-specific acknowledgements in the standalone files.
- Added dated covers, PDF metadata, bookmarks and ten standalone chapter drivers.
  Standalone chapters carry the book title, author, affiliation and date, and keep
  their chapter and figure numbers from the complete books. The full books show
  the date on the cover rather than repeating it on every chapter opening.
- Merged and cleaned the bibliography, removed duplicate entries, repaired
  malformed records, and retained the references actually cited in the revised
  text. Each standalone chapter has its own bibliography.
- Separated chapter text, shared typography, figures, numerical inputs, programs,
  build tools and generated downloads. Retained only used figure assets in the
  working figure directory. Compacted 15 large figure PDF exports without raster
  downsampling, saving approximately 12.8 MB of figure data.
- Replaced the old presentation page by a responsive bilingual page with both
  complete PDFs and all ten standalone chapters. Download metadata is generated
  from the actual PDFs. Preserved the former full-book URLs and maintained the
  four former transport/neural-network download URLs as updated compatibility
  copies. No interactive edition was added.

## Verification

- Both full books and all ten independent chapters compile with **no final TeX
  or bibliography warnings**, undefined references, missing-glyph diagnostics,
  overfull boxes or underfull boxes.
- All twelve PDF files are checked for author metadata, French text extraction,
  embedded fonts, a single edition date, and consistent chapter/figure numbering.
  Pages and chapter openings were also rendered and visually reviewed.
- **33 independent numerical checks** cover the displayed coding strings, entropy
  identities, image operations, a sparse recovery example, the Paris assignment
  and dual certificate, one-dimensional matching, and the neural-network
  gradients and GAN identity. The figure generator also checks its wavelet
  reconstruction, Parseval errors and numerical stopping criterion.
- Built the presentation page with Jekyll and checked all twelve PDF downloads
  over HTTP, twenty local links/assets/anchors, keyboard access, and layouts at
  1360, 390 and 320 pixels wide. No horizontal overflow or browser script errors
  were found. Verified all 335 archived original files against their SHA-256 hashes.
- Machine-readable results are in `reports/`. These numerical and rendering
  checks support the mathematical and editorial review; they do not constitute
  proofs of the general theorems.
