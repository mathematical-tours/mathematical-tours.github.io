# Figure audit guide

Generated from the current publication captions and experiment data. Links open [the interactive book](https://mathematical-tours.github.io/book/interactive/); use the **Book figure** tab to compare with the approved PDF illustration. Numbers are those of the main PDF, including its continued figures.

Sliders select the computed values listed below. Menus select actual inputs or methods; experiment notes explain numerical conventions. The original cover and explicitly retained illustrations remain unchanged.

| Chapter | Interactive | Retained | Computed states |
| --- | ---: | ---: | ---: |
| [1. Shannon Sampling Theory](#shannon) | 4 | 5 | 50 |
| [2. Fourier and Convolution](#fourier) | 4 | 18 | 87 |
| [3. Shannon Coding Theory](#shannon-coding) | 5 | 2 | 145 |
| [4. Wavelets](#wavelets) | 23 | 7 | 525 |
| [5. Linear and Nonlinear Approximation](#approximation) | 23 | 3 | 992 |
| [6. Compression](#compression) | 4 | 2 | 63 |
| [7. Denoising](#denoising) | 32 | 1 | 666 |
| [8. Variational Priors and Regularization](#variational-priors) | 8 | 1 | 296 |
| [9. Inverse Problems](#inverse-problems) | 8 | 2 | 268 |
| [10. Sparse Regularization](#sparse-regularization) | 9 | 2 | 242 |
| [11. Theory of Sparse Regularization](#sparse-theory) | 3 | 8 | 75 |
| [12. Compressed Sensing](#compressed-sensing) | 6 | 0 | 134 |
| [13. Basics of Machine Learning](#machine-learning) | 11 | 10 | 210 |
| [14. Optimization & Machine Learning: Smooth Optimization](#optim-ml-smooth) | 0 | 13 | 0 |
| [15. Optimization & Machine Learning: Advanced Topics](#optim-ml-advanced) | 5 | 11 | 158 |
| [16. Shallow Learning](#perceptrons) | 0 | 0 | 0 |
| [17. Deep Learning](#deep-learning) | 0 | 1 | 0 |
| [18. Convex Analysis](#convex-analysis) | 0 | 4 | 0 |
| [19. Nonsmooth Convex Optimization](#optim-nonsmooth) | 0 | 1 | 0 |

<a id="shannon"></a>

## 1. Shannon Sampling Theory

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 1.1](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--examples-1.html)**<br>Examples of sound ($d=1$), an image ($d=2$), and a video ($d=3$). | **Zoom factor:** 1, 2, 4, 8, 16<br>**Source:** Flower, Birdsong | 10 |
| **[Figure 1.2](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--examples-2.html)**<br>A color flower image and a synthetic multispectral image with 32 channels, with an example spectrum. | **Spectral channel:** 1, 4, 8, 12, 16, 20, 24, 28, 32 | 9 |
| **[Figure 1.3](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--discretization.html)**<br>Spatial and temporal detail: a flower image with a pixelated crop, and a recorded birdsong waveform with a short-time enlargement. Matching red boxes identify the spatial crop. | **Zoom factor:** 1, 2, 4, 8, 16<br>**Source:** Flower, Birdsong | 10 |
| **[Figure 1.4](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--sampling-aliasing.html)**<br>Sampling and Shannon interpolation in the frequency and time domains, with unit sample spacing. Left: the bandlimiting condition of Theorem [1.2](https://mathematical-tours.github.io/book/interactive/shannon/#thm-shannon-sampling) holds. Right: overlapping spectral replicas cause aliasing; the reconstructed signal still passes through every sample. The rows show the original signal, sampling, the interpolation kernel, and reconstruction. | Retained approved TikZ display. | — |
| **[Figure 1.5](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--shannon-sinc.html)**<br>The sinc interpolation kernel. | Retained approved TikZ display. | — |
| **[Figure 1.6](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--sampling-splines.html)**<br>Cardinal splines as basis functions for interpolation. | Retained approved TikZ display. | — |
| **[Figure 1.7](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--aliasing.html)**<br>Aliasing of a sine wave. This periodic example is interpreted through its discrete spectral lines rather than as an $L^2(\mathbb{R})$ signal. | Retained approved TikZ display. | — |
| **[Figure 1.8](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--shannon-quantizer.html)**<br>A uniform scalar quantizer. | Retained approved TikZ display. | — |
| **[Figure 1.9](https://mathematical-tours.github.io/book/interactive/live/panels/shannon--section3-quantize.html)**<br>Quantizing an image with $K\in\{2,3,4,16\}$ gray levels. The normalized intensities lie in $[0,1)$; finite-level quantizers also specify clipping or endpoint conventions. | **Gray levels K:** 2, 3, 4, 8, 16, 32, 64<br>**Image:** Flower, Cartoon, Mandrill | 21 |

<a id="fourier"></a>

## 2. Fourier and Convolution

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 2.1](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-wav.html)**<br>Fourier modes (left) and several translated Daubechies wavelets at each scale (right). | **Scale level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4 | 12 |
| **[Figure 2.2](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-transforms.html)**<br>Four settings for Fourier analysis, linked by sampling and periodization. | Retained approved TikZ display. | — |
| **[Figure 2.3](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-convolution.html)**<br>Convolution on $\mathbb{R}$. The dashed curve is $f\star g$; the marked value $(f\star g)(x)$ equals the shaded integral. | Retained approved TikZ display. | — |
| **[Figure 2.4](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--filter-box.html)**<br>Signal filtering with a box filter (running average). | Retained approved TikZ display. | — |
| **[Figure 2.5](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--filtering-1d.html)**<br>Filtering an irregular signal with Gaussian filters of increasing width $\sigma$. | **Filter width σ (samples):** 0, 0.5, 1, 2, 4, 8, 16<br>**Signal:** Piecewise smooth, Smooth, Chirp | 21 |
| **[Figure 2.6](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-convolution-fourier.html)**<br>How the Fourier transform converts convolution into multiplication. | Retained approved TikZ display. | — |
| **[Figure 2.7](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--splines.html)**<br>Cardinal splines are defined by successive convolutions. | Retained approved TikZ display. | — |
| **[Figure 2.8](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-translation-inv.html)**<br>Commutative diagram for translation invariance. | Retained approved TikZ display. | — |
| **[Figure 2.9](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-poisson-distrib.html)**<br>Sine waves summed in the Poisson formula. | Retained approved TikZ display. | — |
| **[Figure 2.10](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fft.html)**<br>Diagram of one radix-two FFT step. Both halves of the input enter the sum and difference; the outputs of the two smaller transforms are interleaved into even and odd frequencies. | Retained approved TikZ display. | — |
| **[Figure 2.11](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--padding-spacial.html)**<br>Fourier transform approximation by zero-padding in the spatial domain. | Retained approved TikZ display. | — |
| **[Figure 2.12](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--padding-fourier.html)**<br>Interpolation by zero-padding in the frequency domain. | Retained approved TikZ display. | — |
| **[Figure 2.13](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-wave.html)**<br>2-D sine wave. | Retained approved TikZ display. | — |
| **[Figure 2.14](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--2d-fourier-extension.html)**<br>2D Fourier orthogonal bases. | **Horizontal frequency:** 0, 1, 2, 4, 6, 8<br>**Vertical frequency:** 0, 1, 2, 4, 6, 8 | 36 |
| **[Figure 2.15](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-torus.html)**<br>The two-dimensional torus $\mathbb{T}^2=(\mathbb{R}/2\pi\mathbb{Z})^2$. | Retained approved TikZ display. | — |
| **[Figure 2.16](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-torus-discrete.html)**<br>Discrete 2-D torus. | Retained approved TikZ display. | — |
| **[Figure 2.17](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-2d.html)**<br>2-D Fourier analysis of an image (left), and attenuation of the periodicity artifact using masking (right). | **Boundary taper:** 0, 0.2, 0.4, 0.6, 0.8, 1<br>**Image:** Flower, Cartoon, Mandrill | 18 |
| **[Figure 2.18](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--heat.html)**<br>Heat diffusion as a convolution. | Retained approved TikZ display. | — |
| **[Figure 2.19](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-finite-difference-spectrum.html)**<br>Comparison of the spectra of $\Delta$ and $D_2$. | Retained approved TikZ display. | — |
| **[Figure 2.20](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-laplacian-surface.html)**<br>Computing the Laplacian on a surface. | Retained approved TikZ display. | — |
| **[Figure 2.21](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-spherical-coordinates.html)**<br>Spherical coordinates. | Retained approved TikZ display. | — |
| **[Figure 2.22](https://mathematical-tours.github.io/book/interactive/live/panels/fourier--fourier-weighted-graph.html)**<br>Weighted graph. | Retained approved TikZ display. | — |

<a id="shannon-coding"></a>

## 3. Shannon Coding Theory

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 3.1](https://mathematical-tours.github.io/book/interactive/live/panels/shannon-coding--tree-prefix.html)**<br>Left: the complete tree of binary words of length 3. Right: a prefix code. | **Distribution concentration:** 0, 0.25, 0.5, 1, 2, 3, 4<br>**Alphabet size:** 2, 3, 4, 6, 8 | 35 |
| **[Figure 3.2](https://mathematical-tours.github.io/book/interactive/live/panels/shannon-coding--the-concave-function-underlying-entropy.html)**<br>The concave function $h(u)=-u\log_2u$ underlying entropy. | **Distribution concentration:** 0, 0.25, 0.5, 1, 2, 3, 4<br>**Alphabet size:** 2, 3, 4, 6, 8 | 35 |
| **[Figure 3.3](https://mathematical-tours.github.io/book/interactive/live/panels/shannon-coding--histo-entropy.html)**<br>Three examples of probability distributions with corresponding entropies in bits. | **Distribution concentration:** 0, 0.25, 0.5, 1, 2, 3, 4<br>**Alphabet size:** 2, 3, 4, 6, 8 | 35 |
| **[Figure 3.4](https://mathematical-tours.github.io/book/interactive/live/panels/shannon-coding--entropy-extrema.html)**<br>Entropy is minimized at a point mass and maximized at the uniform distribution. Left: level sets of the entropy formula in the positive quadrant and the constraint $p_1+p_2=1$. Middle: its restriction to that constraint. Right: entropy level sets on the three-symbol probability simplex. | Retained approved TikZ display. | — |
| **[Figure 3.5](https://mathematical-tours.github.io/book/interactive/live/panels/shannon-coding--kraft-trees.html)**<br>Panel 1: full binary tree obtained by completing the tree associated with the code $(c_1=0, c_2=10, c_3=110, c_4=111)$. Panel 2: subtrees with the prescribed codeword lengths packed from left to right. | Retained approved TikZ display. | — |
| **[Figure 3.6](https://mathematical-tours.github.io/book/interactive/live/panels/shannon-coding--huff-code.html)**<br>Successive merges in Huffman’s coding algorithm. | **Distribution concentration:** 0, 0.25, 0.5, 1, 2, 3, 4<br>**Alphabet size:** 2, 3, 4, 6, 8 | 35 |
| **[Figure 3.7](https://mathematical-tours.github.io/book/interactive/live/panels/shannon-coding--codage-differences.html)**<br>Difference coding. From left to right: the signal, its first differences, the corresponding histograms, and a code tree for the differences. | **Alphabet size:** 2, 3, 4, 6, 8 | 5 |

<a id="wavelets"></a>

## 4. Wavelets

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 4.1](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--multires-1d.html)**<br>Translations and dilations of a scaling function generate the approximation spaces. The Shannon scaling function $\varphi(t)=\sin(\pi t)/(\pi t)$ is illustrated, with horizontal coordinate $x/2^j$ and heights relative to $2^{-j/2}$. | Retained approved TikZ display. | — |
| **[Figure 4.2](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--haard-scaling-wav.html)**<br>Haar scaling and wavelet functions. Both have amplitude $2^{-j/2}$ and unit $L^2$ norm. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4 | 12 |
| **[Figure 4.3](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--spline-ortho.html)**<br>Orthogonalization of B-splines to define cardinal orthogonal spline functions. | Retained approved TikZ display. | — |
| **[Figure 4.4](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--haar-1d.html)**<br>1-D Haar multiresolution projection $P_{V_j}f$ of a function $f$. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Figure 4.5](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--haar-appr-details.html)**<br>Three successive Haar approximations (top) and their extracted details (bottom). Coarser approximations use fewer intervals. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Figure 4.6](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--spline-vs-shannon.html)**<br>Spline and Shannon wavelet segmentation of the frequency axis. | Retained approved TikZ display. | — |
| **[Figure 4.7](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--periodize.html)**<br>Periodized basis functions. Only the translations $0\leqslant n<2^{-j}$ are distinct, and $j_0=0$ is a common choice of coarsest scale. | Retained approved TikZ display. | — |
| **[Figure 4.8](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavcoef-1d.html)**<br>Wavelet coefficients. The first column shows the input signal and its packed decomposition; the other panels show individual detail scales. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Figure 4.9](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--haar-inner.html)**<br>Haar filter weights as inner products. | Retained approved TikZ display. | — |
| **[Figure 4.10](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--filterbank-fwd.html)**<br>Forward filter bank decomposition. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Figure 4.11](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--filterbank-pyramidal.html)**<br>Pyramidal computation of wavelet coefficients: approximation vectors $a_j$ on the top row and detail vectors $d_j$ below. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Figure 4.12](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavalgo-1d.html)**<br>Packed coefficient arrays after one, two, and three wavelet decomposition steps. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Figure 4.13](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--inversion-transpose.html)**<br>Analysis and inverse wavelet transforms in block-matrix form. Synthesis uses the transpose of the orthogonal analysis operator. | **Decomposition level:** 1, 2<br>**Wavelet:** haar, db2, db4 | 6 |
| **[Figure 4.14](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--haar-analysis-synthesis-matrix.html)**<br>One-step Haar analysis matrix $W$ and synthesis matrix $W^{\mathsf T}$ for eight samples. | **Decomposition level:** 1, 2<br>**Wavelet:** haar, db2, db4 | 6 |
| **[Figure 4.15](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--filterbank-bwd.html)**<br>Inverse filter-bank reconstruction. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Figure 4.16](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavalgo-2d-aniso.html)**<br>Steps of the anisotropic wavelet transform. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 24 |
| **[Figure 4.17](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavecoefs-iso.html)**<br>Anisotropic (left) versus isotropic (right) wavelet coefficients. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 24 |
| **[Figure 4.18](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavelets-2d.html)**<br>2-D wavelets and a schematic support centered at $(2^j n_1,2^j n_2)$, with width $K2^j$ (right). | **Wavelet:** haar, db2, db4 | 3 |
| **[Figure 4.19](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--multires-2d-haar.html)**<br>2-D Haar approximation $P_{V_j^O}f$ for increasing $j$. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 24 |
| **[Figure 4.20](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavcoefs-2d.html)**<br>2-D wavelet coefficients. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 24 |
| **[Figure 4.21](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavcoefs-example.html)**<br>Images (top) and their wavelet coefficients (bottom). | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 24 |
| **[Figure 4.22](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavalgo-2d-step-fwd.html)**<br>Forward 2-D filterbank step. | **Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 6 |
| **[Figure 4.23](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavalgo-2d.html)**<br>One step of the 2-D wavelet transform algorithm. | **Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 6 |
| **[Figure 4.24](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavalgo-2d-step-bwd.html)**<br>Backward 2-D filterbank step. | **Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 6 |
| **[Figure 4.25](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--filter-constr.html)**<br>Constraints on low-pass and high-pass filters. | Retained approved TikZ display. | — |
| **[Figure 4.26](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--vanmoments.html)**<br>Vanishing moments correspond to vanishing derivatives of the Fourier transform. | Retained approved TikZ display. | — |
| **[Figure 4.27](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--vanishing-moments-1d.html)**<br>Location of large coefficients for Haar (left) and Daubechies-4 (right) wavelets applied to the same piecewise-smooth signal. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Figure 4.28](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--wavelets.html)**<br>Daubechies mother wavelets $\psi$ with increasing numbers $p$ of vanishing moments. | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4 | 12 |
| **[Unnumbered illustration](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--embeded-spaces-1d.html)**<br>Unnumbered illustration: embeded spaces 1d | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 36 |
| **[Unnumbered illustration](https://mathematical-tours.github.io/book/interactive/live/panels/wavelets--embeded-spaces-2d.html)**<br>Unnumbered illustration: embeded spaces 2d | **Decomposition level:** 1, 2, 3, 4<br>**Wavelet:** haar, db2, db4<br>**Image:** Flower, Cartoon | 24 |

<a id="approximation"></a>

## 5. Linear and Nonlinear Approximation

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 5.1](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--approximation.html)**<br>Linear versus nonlinear wavelet approximation. | **Retained coefficients M:** 4, 8, 16, 32, 64, 128, 256<br>**Basis:** Fourier, haar, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 63 |
| **[Figure 5.2](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--decay-coefs.html)**<br>Ordered coefficient magnitudes and the threshold used for nonlinear approximation. | Retained approved TikZ display. | — |
| **[Figure 5.3](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--thresholding-hard.html)**<br>Hard thresholding. | Retained approved TikZ display. | — |
| **[Figure 5.4](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--approximation-class.html)**<br>Examples of image models. The bounded-variation example is generated by ROF total-variation denoising of white noise. | **Smoothing strength:** 0, 0.5, 1, 2, 4, 8 | 6 |
| **[Figure 5.5](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--signal-model-sobolev.html)**<br>Progressive Gaussian filtering of one white-noise image. Smoothing widths increase from left to right; the displayed images have the same mean and variance. | **Smoothing strength:** 0, 0.5, 1, 2, 4, 8 | 6 |
| **[Figure 5.6](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--image-model-cartoon.html)**<br>Progressive dead-leaves construction: later opaque shapes cover earlier ones, producing a family of cartoon images. | **Opaque leaves:** 5, 15, 40, 100, 250, 600 | 6 |
| **[Figure 5.7](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--sample-cartoon.html)**<br>A grayscale frame from *Feline Follies* (1919), a mathematical jump along a smooth graph, and a Gaussian-smoothed cartoon. Film frame: Pat Sullivan, via Wikimedia Commons, public domain. | **Smoothing strength:** 0, 0.5, 1, 2, 4, 8<br>**Image:** Flower, Cartoon, Mandrill, Film cartoon | 24 |
| **[Figure 5.8](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--approx-speed-1.html)**<br>Test images with different structures: smooth, cartoon, flower, and mandrill. | **Smoothing strength:** 0, 0.5, 1, 2, 4, 8<br>**Image:** Flower, Cartoon, Mandrill, Film cartoon | 24 |
| **[Figure 5.9](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--approx-speed.html)**<br>Wavelet approximation error for the images in Figure [5.8](https://mathematical-tours.github.io/book/interactive/approximation/#fig-approx-speed-1). | **Retained coefficients M:** 16, 64, 256, 1024, 4096, 8192, 16384<br>**Basis:** Fourier, DCT, haar, db4<br>**Image:** Flower, Cartoon, Mandrill | 84 |
| **[Figure 5.10](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--approx-effi-1.html)**<br>Fourier, DCT, and wavelet approximations of the same $256\times256$ flower image, each retaining $M=4096$ real orthonormal coefficients. | **Retained coefficients M:** 16, 64, 256, 1024, 4096, 8192, 16384<br>**Basis:** Fourier, DCT, haar, db4<br>**Image:** Flower, Cartoon, Mandrill | 84 |
| **[Figure 5.11](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--approx-effi.html)**<br>Comparison of approximation error decay for different bases. | **Retained coefficients M:** 16, 64, 256, 1024, 4096, 8192, 16384<br>**Basis:** Fourier, DCT, haar, db4<br>**Image:** Flower, Cartoon, Mandrill | 84 |
| **[Figure 5.12](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--approx-comparaison.html)**<br>Summary of linear and nonlinear approximation rates $&#124;\!&#124; f-f_M &#124;\!&#124;^2$ for different classes of 1-D signals and images. | Retained approved TikZ display. | — |
| **[Figure 5.13](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--ourier-approx-1d.html)**<br>Fourier approximation of a signal. | **Retained coefficients M:** 4, 8, 16, 32, 64, 128, 256<br>**Basis:** Fourier, haar, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 63 |
| **[Figure 5.14](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--fourier-approx-2d-linnonlin.html)**<br>Linear (top row) and nonlinear (bottom row) Fourier approximation. | **Retained coefficients M:** 16, 64, 256, 1024, 4096, 8192, 16384<br>**Basis:** Fourier, DCT, haar, db4<br>**Image:** Flower, Cartoon, Mandrill | 84 |
| **[Figure 5.15](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--singular-part.html)**<br>Singular and regular regions of a signal (left) and an image (right). | **Smoothing strength:** 0, 0.5, 1, 2, 4, 8<br>**Image:** Flower, Cartoon, Mandrill, Film cartoon | 24 |
| **[Figure 5.16](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--singular-support-1d.html)**<br>Segmentation of the wavelet coefficients into regular and singular parts. | **Retained coefficients M:** 4, 8, 16, 32, 64, 128, 256<br>**Basis:** Fourier, haar, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 63 |
| **[Figure 5.17](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--wav-approx-1d.html)**<br>1-D wavelet approximation. | **Retained coefficients M:** 4, 8, 16, 32, 64, 128, 256<br>**Basis:** Fourier, haar, db4<br>**Signal:** Piecewise smooth, Smooth, Chirp | 63 |
| **[Figure 5.18](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--approx-2d-wav.html)**<br>2-D wavelet approximation. | **Retained coefficients M:** 16, 64, 256, 1024, 4096, 8192, 16384<br>**Basis:** Fourier, DCT, haar, db4<br>**Image:** Flower, Cartoon, Mandrill | 84 |
| **[Figure 5.19](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--cartoon-wav.html)**<br>Wavelet approximation of a cartoon image. | **Retained coefficients M:** 16, 64, 256, 1024, 4096, 8192, 16384<br>**Basis:** Fourier, DCT, haar, db4<br>**Image:** Flower, Cartoon, Mandrill | 84 |
| **[Figure 5.20](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--triangulation-image.html)**<br>A cartoon with smooth variation on both sides of an elliptical jump, a triangulation elongated along that boundary, and the continuous $P_1$ interpolant. | **Vertices around the boundary:** 32, 40, 48, 56, 64<br>**Smooth variation:** 0.4, 0.7, 1, 1.3 | 20 |
| **[Figure 5.21](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--triangulation-anisoiso.html)**<br>The same mesh over a lightened $P_1$ interpolant (left), with the boxed boundary region enlarged (right). Thin triangles follow the red singularity curve; the arrow indicates its tangent. Smooth regions use nearly equilateral elements. | **Vertices around the boundary:** 32, 40, 48, 56, 64<br>**Smooth variation:** 0.4, 0.7, 1, 1.3 | 20 |
| **[Figure 5.22](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--triangulation-peppers.html)**<br>Compression of a piecewise-smooth cartoon using a contour-aligned mesh and JPEG-2000 at the same transmitted bit budget. From left to right: original, mesh, decoded $P_1$ image, and JPEG-2000 image. Mesh costs include coordinates, nodal values, connectivity, and headers. | **Vertices around the boundary:** 32, 40, 48, 56, 64<br>**Smooth variation:** 0.4, 0.7, 1, 1.3 | 20 |
| **[Figure 5.23](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--curvelet.html)**<br>A curvelet $c_m$ in space (left) and its localization in the Fourier domain (right). | **Radial scale:** 4, 8, 16, 32<br>**Orientation (degrees):** 0, 30, 60, 90, 120, 150 | 24 |
| **[Figure 5.24](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--curvelet-discretization.html)**<br>Sampling pattern for the curvelet positions. | **Radial scale:** 4, 8, 16, 32<br>**Orientation (degrees):** 0, 30, 60, 90, 120, 150 | 24 |
| **[Figure 5.25](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--curvelets-vs-wavelets.html)**<br>Wavelet (left) and curvelet (right) approximation of a cartoon image. | **Radial scale:** 4, 8, 16, 32<br>**Orientation (degrees):** 0, 30, 60, 90, 120, 150 | 24 |
| **[Figure 5.26](https://mathematical-tours.github.io/book/interactive/live/panels/approximation--curvelet-denoising.html)**<br>Comparison of translation-invariant wavelet denoising and curvelet denoising. | **Threshold T / σ:** 1, 1.5, 2, 2.5, 3, 3.5, 4<br>**Image:** Flower, Cartoon | 14 |

<a id="compression"></a>

## 6. Compression

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 6.1](https://mathematical-tours.github.io/book/interactive/live/panels/compression--thresholding-vs-quantizing.html)**<br>Hard thresholding and quantization as functions of the input coefficient. | Retained approved TikZ display. | — |
| **[Figure 6.2](https://mathematical-tours.github.io/book/interactive/live/panels/compression--compression.html)**<br>Image compression using wavelet support coding. | **Quantization step T:** 0.025, 0.05, 0.075, 0.1, 0.15, 0.25, 0.4<br>**Image:** Flower, Mandrill | 14 |
| **[Figure 6.3](https://mathematical-tours.github.io/book/interactive/live/panels/compression--entropy.html)**<br>Three different probability distributions on the same alphabet of $Q=7$ symbols. Entropy is measured in bits. | Retained approved TikZ display. | — |
| **[Figure 6.4](https://mathematical-tours.github.io/book/interactive/live/panels/compression--jpeg-2k-architecture.html)**<br>JPEG-2000 coding architecture. | **Target bits / pixel:** 0.3, 0.5, 0.75, 1, 1.5, 2, 3<br>**Image:** Flower, Mandrill | 14 |
| **[Figure 6.5](https://mathematical-tours.github.io/book/interactive/live/panels/compression--compression-jpeg.html)**<br>JPEG and JPEG-2000 reconstructions of the flower (first pair) and mandrill (second pair), at matched target bit rates. | **Target bits / pixel:** 0.3, 0.5, 0.75, 1, 1.5, 2, 3<br>**Image:** Flower, Mandrill | 14 |
| **[Unnumbered illustration](https://mathematical-tours.github.io/book/interactive/live/panels/compression--quantization-bins.html)**<br>Unnumbered illustration: quantization bins | **Gray levels K:** 2, 3, 4, 8, 16, 32, 64<br>**Image:** Flower, Cartoon, Mandrill | 21 |

<a id="denoising"></a>

## 7. Denoising

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 7.1](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--noise-examples.html)**<br>Gaussian, Poisson, and multiplicative Gamma noise simulated on the same flower image. | **Noise standard deviation σ:** 0, 0.02, 0.04, 0.08, 0.12, 0.2<br>**Image:** Flower, Cartoon | 12 |
| **[Figure 7.2](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--additive-noise-model.html)**<br>Image acquisition, noise, and reconstruction by denoising. | **Noise standard deviation σ:** 0, 0.02, 0.04, 0.08, 0.12, 0.2<br>**Image:** Flower, Cartoon | 12 |
| **[Figure 7.3](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--noise-example-1d.html)**<br>One-dimensional additive Gaussian noise at low, medium, and high levels, with a common clean signal and normalized noise realization. | **Noise standard deviation σ:** 0, 0.02, 0.04, 0.08, 0.12, 0.2 | 6 |
| **[Figure 7.4](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--noise-example-2d.html)**<br>Two-dimensional additive noise: the clean image $f_0$ and the observation $f=f_0+w$, with $w\sim\mathcal N(0,\sigma^2 I)$. | **Noise standard deviation σ:** 0, 0.02, 0.04, 0.08, 0.12, 0.2<br>**Image:** Flower, Cartoon | 12 |
| **[Figure 7.5](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--denoising-clean.html)**<br>Left: clean image, center: noisy image, right: denoised image. | **Gaussian width (pixels):** 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5<br>**Image:** Flower, Cartoon | 18 |
| **[Figure 7.6](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--filtering-denoised.html)**<br>Gaussian filtering: noisy signal, kernel, and filtered signal (top), with their Fourier magnitudes below. | **Gaussian width (pixels):** 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5 | 9 |
| **[Figure 7.7](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--filtering-progression.html)**<br>Gaussian filtering with increasing width: three signals (top) and five images (bottom). | **Gaussian width (pixels):** 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5<br>**Image:** Flower, Cartoon | 18 |
| **[Figure 7.8](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--iltering-optimal-curve.html)**<br>SNR as a function of filter width for a signal (left) and an image (right). | **Gaussian width (pixels):** 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5<br>**Image:** Flower, Cartoon | 18 |
| **[Figure 7.9](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--filtering-optimal-1d.html)**<br>Gaussian filtering at the SNR-optimal widths. From left to right: noisy signal, filtered signal, noisy image, and filtered image. | **Gaussian width (pixels):** 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5<br>**Image:** Flower, Cartoon | 18 |
| **[Figure 7.10](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--wavthresh.html)**<br>Denoising using thresholding of wavelet coefficients. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Threshold:** Soft, Hard, Stein<br>**Image:** Flower, Cartoon | 42 |
| **[Figure 7.11](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--thresholding-hard-vs-soft.html)**<br>Hard and soft thresholding functions. | Retained approved TikZ display. | — |
| **[Figure 7.12](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--wavthresh-2d-hard-vs-soft.html)**<br>SNR as a function of $T/\sigma$ for hard and soft thresholding. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Threshold:** Soft, Hard<br>**Image:** Flower, Cartoon | 28 |
| **[Figure 7.13](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--wavthresh-2d.html)**<br>Comparison of hard (left) and soft (right) thresholding. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Threshold:** Soft, Hard, Stein<br>**Image:** Flower, Cartoon | 42 |
| **[Figure 7.14](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--sparse-signal.html)**<br>Left: clean sparse coefficients $a_0$; right: noisy coefficients $a$. | **Noise standard deviation σ:** 0, 0.02, 0.04, 0.08, 0.12, 0.2 | 6 |
| **[Figure 7.15](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--max-gaussian-mean.html)**<br>Empirical mean (top) and standard deviation (bottom) of the largest absolute Gaussian noise coefficient. | **Number of variables N:** 8, 16, 32, 64, 128, 256, 512, 1024 | 8 |
| **[Figure 7.16](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--wavthresh-2d-soft-ti.html)**<br>Soft thresholding in an orthogonal wavelet basis (left) and translation-invariant hard thresholding (right). | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Threshold:** Soft, Hard, Stein<br>**Image:** Flower, Cartoon | 42 |
| **[Figure 7.17](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--wavthresh-ti.html)**<br>SNR as a function of $T/\sigma$ for translation-invariant thresholding. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Threshold:** Soft, Hard, Stein<br>**Image:** Flower, Cartoon | 42 |
| **[Figure 7.18](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--tiwav.html)**<br>Translation invariant wavelet coefficients. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Threshold:** Soft, Hard, Stein<br>**Image:** Flower, Cartoon | 42 |
| **[Figure 7.19](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--tiwav-coefs-thresh.html)**<br>Left: translation invariant wavelet coefficients, for $j=-8, \omega=H$, right: thresholded coefficients. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Threshold:** Soft, Hard, Stein<br>**Image:** Flower, Cartoon | 42 |
| **[Figure 7.20](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--semisoft.html)**<br>Left: semi-soft thresholder, right: Stein thresholder. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Transition μ:** 1.1, 1.5, 2, 3, 5 | 35 |
| **[Figure 7.21](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--semisoft-optimal.html)**<br>Left: SNR as a function of $\mu$ and $T/\sigma$. Right: SNR versus $\mu$, with $T/\sigma$ optimized separately for each $\mu$. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Transition μ:** 1.1, 1.5, 2, 3, 5 | 35 |
| **[Figure 7.22](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--stein-optimal-curve.html)**<br>SNR as a function of $T/\sigma$ for Stein thresholding. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Image:** Flower, Cartoon | 14 |
| **[Figure 7.23](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--block-wavcoefs.html)**<br>Left: wavelet coefficients, center: block thresholded coefficients, right: denoised image. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Block side:** 1, 2, 4 | 21 |
| **[Figure 7.24](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--block-optimal-curves.html)**<br>SNR as a function of $T/\sigma$ (left) and comparison of different block sizes (right). | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Block side:** 1, 2, 4 | 21 |
| **[Figure 7.25](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--comp-thresh.html)**<br>Left: translation invariant wavelet hard thresholding, center: block orthogonal Stein thresholding, right: block translation invariant Stein thresholding. | **Threshold T / σ:** 0.5, 1, 1.5, 2, 2.5, 3, 4<br>**Block side:** 1, 2, 4 | 21 |
| **[Figure 7.26](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--poisson-distrib.html)**<br>Poisson distributions for various $\lambda$. | **Peak photon count:** 1, 2, 4, 8, 16, 32, 64, 128 | 8 |
| **[Figure 7.27](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--poisson-denoising.html)**<br>Poisson-corrupted flower images at increasing photon-count scales $\lambda_{\max}$. | **Peak photon count:** 1, 2, 4, 8, 16, 32, 64, 128<br>**Image:** Flower, Cartoon | 16 |
| **[Figure 7.28](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--variance-stabilization.html)**<br>Exact Poisson variances of the Anscombe, Freeman–Tukey, and $2\sqrt{x}$ transforms as functions of the mean intensity. | **Peak photon count:** 1, 2, 4, 8, 16, 32, 64, 128<br>**Image:** Flower, Cartoon | 16 |
| **[Figure 7.29](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--poisson-wav.html)**<br>Left: noisy image, center: denoising without variance stabilization, right: denoising after variance stabilization. | **Peak photon count:** 1, 2, 4, 8, 16, 32, 64, 128<br>**Image:** Flower, Cartoon | 16 |
| **[Figure 7.30](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--multnoise.html)**<br>Multiplicative Gamma noise for increasing numbers $K$ of averaged looks, with standard deviation $\sigma=K^{-1/2}$. | **Independent looks K:** 1, 2, 4, 8, 16, 32<br>**Image:** Flower, Cartoon | 12 |
| **[Figure 7.31](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--multnoise-hist-unstab.html)**<br>Histogram of multiplicative noise before (left) and after (right) stabilization. | **Independent looks K:** 1, 2, 4, 8, 16, 32 | 6 |
| **[Figure 7.32](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--multnoise-wav.html)**<br>Left: noisy image, center: denoising after variance stabilization, right: denoising without variance stabilization. | **Independent looks K:** 1, 2, 4, 8, 16, 32<br>**Image:** Flower, Cartoon | 12 |
| **[Unnumbered illustration](https://mathematical-tours.github.io/book/interactive/live/panels/denoising--cycle-spin-principle.html)**<br>Unnumbered illustration: cycle spin principle | **Translation (pixels):** 0, 1, 2, 3, 4, 5, 6, 7<br>**Method:** Orthogonal, Cycle spinning | 16 |

<a id="variational-priors"></a>

## 8. Variational Priors and Regularization

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 8.1](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--grad-discrete.html)**<br>Discrete gradient vectors. | **Display contrast gain:** 0.5, 1, 2, 4, 8, 16<br>**Image:** Flower, Cartoon, Mandrill | 18 |
| **[Figure 8.2](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--discrete-operators.html)**<br>Discrete gradient of the flower image. Red and blue encode the signed horizontal and vertical components with enhanced display contrast; neutral gray represents zero. The last panel shows the Euclidean gradient magnitude. | **Display contrast gain:** 0.5, 1, 2, 4, 8, 16<br>**Image:** Flower, Cartoon, Mandrill | 18 |
| **[Figure 8.3](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--flow.html)**<br>Heat flow (top) and TV flow (bottom): images $f_t$ at increasing times $t$. | **Evolution / regularization strength:** 0, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4<br>**Model:** Heat flow, TV flow, Sobolev regularization, TV regularization<br>**Image:** Flower, Cartoon | 56 |
| **[Figure 8.4](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--discrete-operators-lapl.html)**<br>Discrete Laplacian (left) and negative discrete TV gradient (right). | **Smoothing ε:** 0.001, 0.003, 0.01, 0.03, 0.1, 0.3<br>**Image:** Flower, Cartoon, Mandrill | 18 |
| **[Figure 8.5](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--regul-abs.html)**<br>Regularized absolute value $x \mapsto \sqrt{x^2+\varepsilon^2}$. | Retained approved TikZ display. | — |
| **[Figure 8.6](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--discrete-operators-tveps.html)**<br>Regularized gradient norm $\sqrt{ &#124;\!&#124; \nabla f(x) &#124;\!&#124;^2 + \varepsilon^2 }$. | **Smoothing ε:** 0.001, 0.003, 0.01, 0.03, 0.1, 0.3<br>**Image:** Flower, Cartoon, Mandrill | 18 |
| **[Figure 8.7](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--denoising-flow.html)**<br>Denoised images $f_t$ at several times $t$: Sobolev flow (top) and TV flow (bottom). | **Evolution / regularization strength:** 0, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4<br>**Model:** Heat flow, TV flow, Sobolev regularization, TV regularization<br>**Image:** Flower, Cartoon | 56 |
| **[Figure 8.8](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--denoising-cameraman.html)**<br>Denoising using PDE flows and regularization. | **Evolution / regularization strength:** 0, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4<br>**Model:** Heat flow, TV flow, Sobolev regularization, TV regularization<br>**Image:** Flower, Cartoon | 56 |
| **[Figure 8.9](https://mathematical-tours.github.io/book/interactive/live/panels/variational-priors--denoising-cameraman-snr.html)**<br>SNR as a function of time $t$ for flows (top) and $\lambda$ for regularization (bottom). | **Evolution / regularization strength:** 0, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4<br>**Model:** Heat flow, TV flow, Sobolev regularization, TV regularization<br>**Image:** Flower, Cartoon | 56 |

<a id="inverse-problems"></a>

## 9. Inverse Problems

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 9.1](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--ip.html)**<br>Examples of forward operators in inverse problems. | **Iterations:** 0, 10, 30, 80, 200, 500<br>**Prior:** Sobolev, TV, Wavelet<br>**Missing region:** Square, Stripes, Random | 54 |
| **[Figure 9.2](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--bound-regul-variance.html)**<br>Left: capping reciprocal singular values. Right: the Tikhonov bound $\mu_{\lambda}(\sigma) \leqslant C_\lambda= \frac{1}{2\sqrt{\lambda}}$. | Retained approved TikZ display. | — |
| **[Figure 9.3](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--ip-bound-regul.html)**<br>Bounding $\lambda\frac{\sigma^\beta}{\lambda+\sigma^2} \leqslant D_{\lambda,\beta}$. | Retained approved TikZ display. | — |
| **[Figure 9.4](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--inp-sob-iter.html)**<br>Sobolev projected gradient descent algorithm. | **Iterations:** 0, 10, 30, 80, 200, 500<br>**Prior:** Sobolev, TV, Wavelet<br>**Missing region:** Square, Stripes, Random | 54 |
| **[Figure 9.5](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--inpainting-parrot.html)**<br>Sobolev inpainting of a grid-shaped occlusion on the flower image. | **Iterations:** 0, 10, 30, 80, 200, 500<br>**Prior:** Sobolev, TV, Wavelet<br>**Missing region:** Square, Stripes, Random | 54 |
| **[Figure 9.6](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--inpainting-cameraman.html)**<br>Inpainting with Sobolev and TV regularization. | **Iterations:** 0, 10, 30, 80, 200, 500<br>**Prior:** Sobolev, TV, Wavelet<br>**Missing region:** Square, Stripes, Random | 54 |
| **[Figure 9.7](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--tomo-principle.html)**<br>The modified Shepp–Logan phantom and one Radon projection at angle $\theta=\pi/6$. The marked line and point correspond to the same detector coordinate. | **Projection angle θ (degrees):** 0, 15, 30, 45, 60, 75, 90, 120, 150, 180 | 10 |
| **[Figure 9.8](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--tomo-radon-subsample.html)**<br>Tomography test phantom, discrete Radon projections, and a Cartesian radial Fourier-sampling mask. | **Radial Fourier lines:** 4, 8, 13, 20, 32, 48, 64<br>**Reconstruction:** Pseudoinverse, TV | 14 |
| **[Figure 9.9](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--tomo-pinv.html)**<br>Pseudoinverse reconstruction from partial Cartesian Fourier measurements on 13 and 32 radial lines. | **Radial Fourier lines:** 4, 8, 13, 20, 32, 48, 64<br>**Reconstruction:** Pseudoinverse, TV | 14 |
| **[Figure 9.10](https://mathematical-tours.github.io/book/interactive/live/panels/inverse-problems--tomo-tv.html)**<br>TV-regularized reconstruction and the pseudoinverse for the same partial Fourier measurements. | **Radial Fourier lines:** 4, 8, 13, 20, 32, 48, 64<br>**Reconstruction:** Pseudoinverse, TV | 14 |

<a id="sparse-regularization"></a>

## 10. Sparse Regularization

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 10.1](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--sparsity-wav.html)**<br>Most wavelet coefficients of a natural image have small magnitude. | **Display / selection threshold:** 0, 0.02, 0.05, 0.1, 0.2, 0.4<br>**Image:** Flower, Cartoon, Mandrill | 18 |
| **[Figure 10.2](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--sparsity-lq.html)**<br>$\ell ^q$ balls $\left\{ x \;;\; J_q(x) \leqslant 1 \right\}$ for varying $q$. | **Exponent q:** 0.3, 0.5, 0.75, 1, 1.5, 2, 4, 8<br>**Radius:** 0.5, 0.75, 1, 1.25, 1.5 | 40 |
| **[Figure 10.3](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--varspars.html)**<br>Panel 1: the objective $&#124;\!&#124; \cdot-y &#124;\!&#124;^2+T^2 &#124;\!&#124; \cdot &#124;\!&#124;_0$. Panels 2–5: the scalar objective $F(x) \mathrel{:=}\frac{1}{2}&#124;x-y&#124;^2+\lambda&#124;x&#124;$ for increasing $\lambda$. | Retained approved TikZ display. | — |
| **[Figure 10.3](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--varspars-continued.html)**<br>Panel 1: the objective $&#124;\!&#124; \cdot-y &#124;\!&#124;^2+T^2 &#124;\!&#124; \cdot &#124;\!&#124;_0$. Panels 2–5: the scalar objective $F(x) \mathrel{:=}\frac{1}{2}&#124;x-y&#124;^2+\lambda&#124;x&#124;$ for increasing $\lambda$. | Retained approved TikZ display. | — |
| **[Figure 10.4](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--optim-geom.html)**<br>Geometry of convex optimization: smooth-boundary and corner contact with objective level sets, followed by the smallest scaled $\ell^2$ and $\ell^1$ balls touching an affine constraint. Red dots mark the contact points. | **Exponent q:** 0.3, 0.5, 0.75, 1, 1.5, 2, 4, 8<br>**Affine right-hand side:** 0.5, 0.75, 1, 1.25, 1.5 | 40 |
| **[Figure 10.5](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--spikes-lun.html)**<br>Sparse spike recovery by the pseudoinverse and by $\ell^1$ regularization. | **Regularization λ:** 1e-05, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1<br>**Algorithm:** ISTA, FISTA | 18 |
| **[Figure 10.6](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--spikes-ista-decay.html)**<br>Convergence of the energy and iterates under iterative soft thresholding. | **Regularization λ:** 1e-05, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1<br>**Algorithm:** ISTA, FISTA | 18 |
| **[Figure 10.7](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--deconvolution-1d.html)**<br>One-dimensional deconvolution with sparsity in an orthogonal wavelet basis. | **Regularization λ:** 1e-05, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1<br>**Algorithm:** ISTA, FISTA | 18 |
| **[Figure 10.8](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--deconvolution-2d.html)**<br>Image deconvolution. | **Regularization λ:** 1e-05, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1<br>**Algorithm:** ISTA, FISTA | 18 |
| **[Figure 10.9](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--deconvolution-2d-snr.html)**<br>Reconstruction SNR as a function of the regularization parameter $\lambda$, swept logarithmically from $10^{-5}$ to $1$. Red dots mark the measured maxima. | **Regularization λ:** 1e-05, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1<br>**Algorithm:** ISTA, FISTA | 18 |
| **[Figure 10.10](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-regularization--inpainting-lena.html)**<br>Inpainting with Sobolev and wavelet sparsity priors. | **Iterations:** 0, 10, 30, 80, 200, 500<br>**Prior:** Sobolev, TV, Wavelet<br>**Missing region:** Square, Stripes, Random | 54 |

<a id="sparse-theory"></a>

## 11. Theory of Sparse Regularization

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 11.1](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--homotopy.html)**<br>Lasso solution paths $\lambda\mapsto x_\lambda$ for sparsities $s=3,6,13$, with the same measurement matrix and noise realization. | **Support size s:** 2, 3, 5, 8, 12<br>**Measurements P:** 16, 24, 32, 40, 48 | 25 |
| **[Figure 11.2](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--polytope-proof.html)**<br>Geometry of the polytope characterization of exact recovery. | Retained approved TikZ display. | — |
| **[Figure 11.3](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--polytopes.html)**<br>The linear map $A$ projects the $\ell^1$ ball $B$; the nonlinear recovery map $\mathcal{A}$ selects solutions of [($\mathcal  {P}_0(y)$)](https://mathematical-tours.github.io/book/interactive/sparse-theory/#eq-lasso-constr). | Original retained at the author’s request. | — |
| **[Figure 11.4](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--sparse-theory-injectivity.html)**<br>Coefficient trajectories along a support-reducing direction. | Retained approved TikZ display. | — |
| **[Figure 11.5](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--bregman.html)**<br>Visualization of Bregman divergences. | Retained approved TikZ display. | — |
| **[Figure 11.6](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--sparse-theory-bregman-l1.html)**<br>Bregman divergence controls the $\ell^1$ error on coordinates where $\eta$ does not saturate. | Retained approved TikZ display. | — |
| **[Figure 11.7](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--delaunay.html)**<br>For unit-norm columns, the supporting hyperplane defined by a valid certificate determines an empty spherical cap: its interior contains no signed column $\pm a_i$. | Original retained at the author’s request. | — |
| **[Figure 11.8](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--sparse-theory-recovery-region.html)**<br>Region in the $(\lambda,\delta)$ plane where sign consistency holds. | Retained approved TikZ display. | — |
| **[Figure 11.9](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--certif-cs.html)**<br>Fuchs precertificate $\eta_F$ for a Gaussian matrix $A \in \mathbb{R}^{P \times N}$ with $N=64$ and increasing measurement count $P$. | **Support size s:** 2, 3, 5, 8, 12<br>**Measurements P:** 16, 24, 32, 40, 48 | 25 |
| **[Figure 11.10](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--spikes.html)**<br>Convolution operator. | Retained approved TikZ display. | — |
| **[Figure 11.11](https://mathematical-tours.github.io/book/interactive/live/panels/sparse-theory--certif-convol.html)**<br>Continuous representations of the precertificates $\eta_F$ and $\eta_V$ for a convolution operator $A$. | **Support size s:** 2, 3, 5, 8, 12<br>**Measurements P:** 16, 24, 32, 40, 48 | 25 |

<a id="compressed-sensing"></a>

## 12. Compressed Sensing

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 12.1](https://mathematical-tours.github.io/book/interactive/live/panels/compressed-sensing--single-pixel.html)**<br>Simulated compressed sensing of the flower image, with $Q=4096$ pixels. The three wavelet-$\ell^1$ reconstructions use $P=512$, $1024$, and $2048$ orthonormal signed Walsh–Hadamard measurements, respectively. | **Measurements P:** 256, 512, 768, 1024, 1536, 2048, 3072, 4096<br>**Image:** Flower, Cartoon | 16 |
| **[Figure 12.2](https://mathematical-tours.github.io/book/interactive/live/panels/compressed-sensing--marcenko-pastur.html)**<br>Empirical eigenvalue distributions approaching the Marchenko–Pastur law. | **Aspect ratio β = s / P:** 0.1, 0.2, 0.4, 0.6, 0.8<br>**Rows P:** 40, 80, 160, 320 | 20 |
| **[Figure 12.3](https://mathematical-tours.github.io/book/interactive/live/panels/compressed-sensing--marcenko-pastur-laws.html)**<br>Marchenko–Pastur densities $f_\beta$ for several aspect ratios $\beta$. | **Aspect ratio β = s / P:** 0.1, 0.2, 0.4, 0.6, 0.8<br>**Rows P:** 40, 80, 160, 320 | 20 |
| **[Figure 12.4](https://mathematical-tours.github.io/book/interactive/live/panels/compressed-sensing--max-gaussian.html)**<br>Maximum of $N$ independent standard Gaussian variables, compared with the dashed curve $\sqrt{2\log(N)}$. | **Number of variables N:** 8, 16, 32, 64, 128, 256, 512, 1024 | 8 |
| **[Figure 12.5](https://mathematical-tours.github.io/book/interactive/live/panels/compressed-sensing--phase-trans.html)**<br>Phase transitions. Right: empirical probabilities that recovery criteria hold as sparsity varies. Blue: weak ERC; black: ERC; green: $&#124;\!&#124; \eta_F &#124;\!&#124;_\infty\leqslant 1$; red: exact recovery. | **Measurements P:** 12, 20, 28, 36, 44<br>**Support size s:** 1, 2, 4, 6, 8, 10, 12 | 35 |
| **[Figure 12.6](https://mathematical-tours.github.io/book/interactive/live/panels/compressed-sensing--rip-const.html)**<br>Evolution of lower bounds $\hat\delta_k$ on the RIP constant. | **Measurements P:** 12, 20, 28, 36, 44<br>**Support size s:** 1, 2, 4, 6, 8, 10, 12 | 35 |

<a id="machine-learning"></a>

## 13. Basics of Machine Learning

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 13.1](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--cov.html)**<br>Empirical covariance of the Iris measurements and its eigenvalues. | **Second principal component:** 2, 3, 4<br>**Dataset:** Iris, Digits | 6 |
| **[Figure 13.2](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--ml-pca-variance.html)**<br>Principal directions of variation. | Retained approved TikZ display. | — |
| **[Figure 13.3](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--pca.html)**<br>PCA projections of the data into two and three dimensions. | **Second principal component:** 2, 3, 4<br>**Dataset:** Iris, Digits | 6 |
| **[Figure 13.4](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--ml-pca-data-matrix.html)**<br>Rows of the centered data matrix are observations; columns are features. | Retained approved TikZ display. | — |
| **[Figure 13.5](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--ml-pca-compression-reconstruction.html)**<br>Linear compression followed by reconstruction. | Retained approved TikZ display. | — |
| **[Figure 13.6](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--pca-var-proof.html)**<br>Panel 1: the linear-program bound in [(13.6)](https://mathematical-tours.github.io/book/interactive/machine-learning/#eq-variational-pca). Panel 2: matrix $B$. Panel 3: its orthogonal extension $\tilde B$. | Retained approved TikZ display. | — |
| **[Figure 13.7](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--ml-voronoi.html)**<br>$k$-means clusters according to Voronoi cells. | Retained approved TikZ display. | — |
| **[Figure 13.8](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--kmeans.html)**<br>Left: iterations of the $k$-means algorithm. Right: class histograms after the $k$-means optimization. | **Lloyd iteration:** 0, 1, 2, 3, 5, 10, 20, 40<br>**Clusters K:** 2, 3, 4, 6, 8 | 40 |
| **[Figure 13.9](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--lloyd.html)**<br>Continuous $k$-means (Lloyd) iterations $0,2,3,5,30$. The top row uses a uniform density; the bottom row uses the nonuniform density shown in grayscale. | **Lloyd iteration:** 0, 1, 2, 3, 5, 10, 20, 40<br>**Clusters K:** 2, 3, 4, 6, 8 | 40 |
| **[Figure 13.10](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--ml-joint-model.html)**<br>Probabilistic modeling. | Retained approved TikZ display. | — |
| **[Figure 13.11](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--bound-regul.html)**<br>Conditional expectation. | Retained approved TikZ display. | — |
| **[Figure 13.12](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--ml-linear-fit.html)**<br>Linear regression. | Retained approved TikZ display. | — |
| **[Figure 13.13](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--ml-nearest-neighbors.html)**<br>Nearest neighbors. | Retained approved TikZ display. | — |
| **[Figure 13.14](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--hist-classif.html)**<br>Classification boundaries produced by $R$-nearest neighbors. | **Neighbors R:** 1, 3, 5, 9, 15, 31, 61<br>**Geometry:** Moons, Circles | 14 |
| **[Figure 13.15](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--logistic-sharpness.html)**<br>Logistic classification in one and two dimensions, showing how $&#124;\!&#124; \beta &#124;\!&#124;$ controls the width of the probability transition. | Retained approved TikZ display. | — |
| **[Figure 13.16](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--losses.html)**<br>Classification losses as functions of the negative margin $s=-y\langle x,\,\beta\rangle$, with zero margin counted as an error. The unscaled logistic loss need not upper-bound the 0–1 loss; division by $\log 2$ gives the upper-bound normalization used above. | **Negative margin s:** -3, -2, -1, 0, 0.5, 1, 2, 3 | 8 |
| **[Figure 13.17](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--separation-influ.html)**<br>Effect of class separation on predicted classification probabilities. | **Class separation:** 0.2, 0.5, 1, 1.5, 2, 3<br>**Inverse regularization C:** 0.1, 1, 10, 100 | 24 |
| **[Figure 13.18](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--digits.html)**<br>Digit images and their two- and three-dimensional PCA projections. | **Second principal component:** 2, 3, 4<br>**Dataset:** Iris, Digits | 6 |
| **[Figure 13.19](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--digits-classes.html)**<br>Digit classification. Left: the first nine class probabilities on an affine two-dimensional PCA slice through the training mean. Right: the mixture of all ten class colors, with projected observations overlaid. | **Second principal component:** 2, 3, 4<br>**Dataset:** Iris, Digits | 6 |
| **[Figure 13.20](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--kernel.html)**<br>Regression using a Gaussian kernel. | **Kernel width σ:** 0.05, 0.1, 0.2, 0.4, 0.8, 1.6<br>**Regularization λ:** 0.0001, 0.001, 0.01, 0.1, 1 | 30 |
| **[Figure 13.21](https://mathematical-tours.github.io/book/interactive/live/panels/machine-learning--classes-kernel.html)**<br>Nonlinear classification using a Gaussian kernel. | **Kernel width σ:** 0.05, 0.1, 0.2, 0.4, 0.8, 1.6<br>**Regularization λ:** 0.0001, 0.001, 0.01, 0.1, 1 | 30 |

<a id="optim-ml-smooth"></a>

## 14. Optimization & Machine Learning: Smooth Optimization

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 14.1](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--ml-ex.html)**<br>Panels 1 and 2: linear regression and a linear classifier. Panel 3: the 0-1 loss and its normalized logistic and hinge upper bounds, both tight at zero margin. | Retained approved TikZ display. | — |
| **[Figure 14.2](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--minimizer-exists.html)**<br>Left: no minimizer. Middle: multiple minimizers. Right: a unique minimizer. | Retained approved TikZ display. | — |
| **[Figure 14.3](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--least-square.html)**<br>Coercivity condition for least squares. | Retained approved TikZ display. | — |
| **[Figure 14.4](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--cvx-vs-noncvx.html)**<br>Convex and nonconvex functions; strictly convex and convex but not strictly convex functions. Here $x_t=(1-t)x+ty$ and $\ell_t=(1-t)f(x)+tf(y)$, with $0<t<1$. | Retained approved TikZ display. | — |
| **[Figure 14.5](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--cvx-set.html)**<br>Comparison of convex functions $f : \mathbb{R}^p \rightarrow \mathbb{R}$ (for $p=1$) and convex sets $C \subset \mathbb{R}^p$ (for $p=2$). | Retained approved TikZ display. | — |
| **[Figure 14.6](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--smooth-gradient-field.html)**<br>The gradient as a vector field. | Retained approved TikZ display. | — |
| **[Figure 14.7](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--first-order.html)**<br>Local extrema (panel 1), a stationary inflection point (panel 2), and a global minimizer (panel 3). | Retained approved TikZ display. | — |
| **[Figure 14.8](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--link-pca.html)**<br>Left: point clouds $(a_i)_i$ and their principal directions. Right: the quadratic part of $f(x)$. | Retained approved TikZ display. | — |
| **[Figure 14.9](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--expansion-taylor.html)**<br>Panels 1 and 2: first-order Taylor approximations in one and two dimensions. Panels 3 and 4: the gradient is orthogonal to the level set; the diagram illustrates the proof. | Retained approved TikZ display. | — |
| **[Figure 14.9](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--expansion-taylor-continued.html)**<br>Panels 1 and 2: first-order Taylor approximations in one and two dimensions. Panels 3 and 4: the gradient is orthogonal to the level set; the diagram illustrates the proof. | Retained approved TikZ display. | — |
| **[Figure 14.10](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--gradesc.html)**<br>Effect of the step size $\tau$ on gradient descent (panel 1) and exact line search (panel 2). | Retained approved TikZ display. | — |
| **[Figure 14.11](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--grad-desc-contract.html)**<br>Contraction constant $h(\tau)$ for a quadratic function. | Retained approved TikZ display. | — |
| **[Figure 14.12](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-smooth--smooth-curvature-bounds.html)**<br>Quadratic upper and lower bounds for a smooth strongly convex function. | Retained approved TikZ display. | — |

<a id="optim-ml-advanced"></a>

## 15. Optimization & Machine Learning: Advanced Topics

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 15.1](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--ml-sparsity-lq.html)**<br>$\ell ^q$ balls $\left\{ x \;;\; \sum _k &#124;x_k&#124;^q \leqslant 1 \right\}$ for varying $q$. | **Exponent q:** 0.3, 0.5, 0.75, 1, 1.5, 2, 4, 8<br>**Radius:** 0.5, 0.75, 1, 1.25, 1.5 | 40 |
| **[Figure 15.2](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--ml-varspars.html)**<br>Effect of $\lambda$ on the penalized scalar objective $F(x)\mathrel{:=}\frac12(x-y)^2+\lambda&#124;x&#124;$. | Retained approved TikZ display. | — |
| **[Figure 15.2](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--ml-varspars-continued.html)**<br>Effect of $\lambda$ on the penalized scalar objective $F(x)\mathrel{:=}\frac12(x-y)^2+\lambda&#124;x&#124;$. | Retained approved TikZ display. | — |
| **[Figure 15.3](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--bgd.html)**<br>Objective-error evolution for batch gradient descent in logistic classification. | **Step-size multiplier:** 0.1, 0.3, 0.6, 1, 1.5, 2<br>**Data passes:** 5, 10, 20, 40, 60 | 30 |
| **[Figure 15.4](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--advanced-unbiased-gradient.html)**<br>Unbiased stochastic gradient estimate (panel 1) and a schematic SGD trajectory (panel 2). | Retained approved TikZ display. | — |
| **[Figure 15.5](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--sgd-traject.html)**<br>Repeated SGD runs for $k\mapsto x_k\in\mathbb{R}$. Left: the distribution of iterates at each iteration. Right: individual trajectories. | **Step size τ:** 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.5<br>**Gradient noise σ:** 0.1, 0.3, 0.5, 1 | 28 |
| **[Figure 15.6](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--sgd.html)**<br>Objective-error evolution for SGD in logistic classification; the dashed curve shows batch gradient descent. | **Step-size multiplier:** 0.1, 0.3, 0.6, 1, 1.5, 2<br>**Data passes:** 5, 10, 20, 40, 60 | 30 |
| **[Figure 15.7](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--compariso-sgd.html)**<br>Objective gaps for SGD, averaged SGD (SGA), and SAG, displayed on logarithmic axes. | **Step-size multiplier:** 0.1, 0.3, 0.6, 1, 1.5, 2<br>**Data passes:** 5, 10, 20, 40, 60 | 30 |
| **[Figure 15.8](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--compgraph.html)**<br>A computational graph. | Original retained at the author’s request. | — |
| **[Figure 15.9](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--forward-backward.html)**<br>Derivative dependencies in forward mode (left) and reverse mode (right). | Original retained at the author’s request. | — |
| **[Figure 15.10](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--dag-example-simple.html)**<br>Example of a simple computational graph. | Original retained at the author’s request. | — |
| **[Figure 15.11](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--matrix-mult.html)**<br>Matrix multiplication costs for forward mode (left) and reverse mode (right) along a chain. | Original retained at the author’s request. | — |
| **[Figure 15.12](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--mlp.html)**<br>Computational graph for a feedforward architecture. | Original retained at the author’s request. | — |
| **[Figure 15.13](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--mlp-param.html)**<br>Multilayer-perceptron parameterization. | Original retained at the author’s request. | — |
| **[Figure 15.14](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--recur.html)**<br>Computational graph for a recurrent architecture. | Original retained at the author’s request. | — |
| **[Figure 15.15](https://mathematical-tours.github.io/book/interactive/live/panels/optim-ml-advanced--recurrent-param.html)**<br>Recurrent residual perceptron parameterization. | Original retained at the author’s request. | — |

<a id="perceptrons"></a>

## 16. Shallow Learning

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |

<a id="deep-learning"></a>

## 17. Deep Learning

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 17.1](https://mathematical-tours.github.io/book/interactive/live/panels/deep-learning--deep-fc-cnn.html)**<br>Panel 1: fully connected network. Panel 2: convolutional neural network. | Retained approved TikZ display. | — |

<a id="convex-analysis"></a>

## 18. Convex Analysis

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 18.1](https://mathematical-tours.github.io/book/interactive/live/panels/convex-analysis--convex-examples.html)**<br>Convexity and strict convexity for functions and sets. | Retained approved TikZ display. | — |
| **[Figure 18.2](https://mathematical-tours.github.io/book/interactive/live/panels/convex-analysis--convex-analysis-subdifferential.html)**<br>The subdifferential | Retained approved TikZ display. | — |
| **[Figure 18.3](https://mathematical-tours.github.io/book/interactive/live/panels/convex-analysis--subdiff-l1.html)**<br>Subdifferential of the absolute value and a piecewise affine convex function. | Retained approved TikZ display. | — |
| **[Figure 18.4](https://mathematical-tours.github.io/book/interactive/live/panels/convex-analysis--convex-analysis-normal-cones.html)**<br>Normal cones | Retained approved TikZ display. | — |

<a id="optim-nonsmooth"></a>

## 19. Nonsmooth Convex Optimization

| Publication figure and exact caption | Controls and available values | States |
| --- | --- | ---: |
| **[Figure 19.1](https://mathematical-tours.github.io/book/interactive/live/panels/optim-nonsmooth--prox-proj.html)**<br>Proximal map and projection map. | Retained approved TikZ display. | — |
