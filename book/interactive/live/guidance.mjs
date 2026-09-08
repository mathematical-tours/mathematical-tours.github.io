/** Reading prompts follow the mathematics of each experiment, not its filename. */
const G=(prompt,observe,focus)=>({prompt,observe,focus});
const families={
  acquisition:G('Zoom into the same observation, then switch between the flower and the birdsong.', 'A closer view reveals the discrete samples behind an apparently continuous image or waveform.', 'zoom'),
  quantization:G('Increase the number of gray levels and watch the quantization steps shrink.', 'The midpoint rule has a maximum error of 1/(2K). More levels reduce the visible intensity bands.', 'K'),
  fourier:G('Change the dyadic scale and compare global Fourier modes with translated wavelets.', 'Wavelets concentrate information in space; Fourier modes extend across the whole domain.', 'scale'),
  coding:G('Compare a uniform distribution with a concentrated one, keeping the alphabet size fixed.', 'Entropy is an ideal average bit count. Huffman lengths must be integers, so the average length usually lies above entropy.', 'skew'),
  wavelet1d:G('Increase the decomposition depth and follow where the signal’s detail coefficients appear.', 'Coarse coefficients describe broad variations; detail coefficients record what each coarser approximation loses.', 'level'),
  wavelet2d:G('Explore another scale or wavelet and follow the colored subband boundaries.', 'The three detail subbands respond to different edge orientations. The coarse block occupies the upper-left corner.', 'level'),
  models:G('Increase the smoothing strength and compare broad variations with sharp transitions.', 'Smoothness within a region and a jump across its boundary are different kinds of image structure.', 'strength'),
  approximation1d:G('Compare linear and largest-coefficient approximations at the same coefficient budget.', 'Keeping the M largest magnitudes minimizes squared error in the chosen orthonormal basis. The best basis still depends on the signal.', 'M'),
  approximation2d:G('Increase the coefficient budget, then compare bases without changing the image.', 'Compare reconstruction quality at equal M. A larger SNR means a smaller squared reconstruction error.', 'M'),
  mesh:G('Add vertices along the curved boundary and inspect how the thin triangles follow its tangent.', 'A narrow direction resolves the jump; a long tangential direction avoids spending as many vertices along the edge.', 'angular'),
  directional:G('Rotate the selected frequency window and examine the corresponding spatial atom.', 'The frequency direction is normal to the atom’s long axis. Finer scales make the atom narrower and more elongated.', 'angle'),
  compression:G('Change the target rate and compare the decoded images at their reported actual rates.', 'Small files spend part of their budget on headers. The achieved bit rate is more informative than the requested target alone.', 'rate'),
  noise:G('Increase the noise level while keeping the clean input fixed.', 'The same normalized noise sample is scaled by σ. Additive Gaussian noise has constant variance across pixels, before display clipping.', 'sigma'),
  gaussian:G('Move from no filtering to stronger smoothing and look for the best compromise.', 'Smoothing suppresses noise and also removes detail. The clean input lets us measure this tradeoff using SNR.', 'width'),
  threshold:G('Raise the threshold, then compare shrinkage rules at the same noise level.', 'A low threshold leaves noise; a high threshold removes meaningful coefficients. Soft shrinkage also biases surviving coefficients.', 'factor'),
  poisson:G('Increase the photon count while keeping the underlying image fixed.', 'Relative shot noise decreases as the photon budget grows. At low counts, a Gaussian noise model is a poor approximation.', 'photons'),
  multiplicative:G('Increase the number of independent looks and compare the fluctuations.', 'The multiplier has mean 1 and variance 1/K. Averaging more looks reduces relative noise while preserving the expected intensity.', 'looks'),
  gradient:G('Increase the display gain to reveal small derivatives, then inspect another image.', 'Red and blue encode signed components. Gain changes visibility, not the finite-difference operator or the computed gradient.', 'gain'),
  pde:G('Compare the smoothing models on the same input, then vary their strength.', 'Diffusion evolves an initial image; regularization balances fidelity to that image against a prior. Their strength parameters have different meanings.', 'strength'),
  inpainting:G('Change the missing region and compare priors at the same iteration count.', 'Each prior fills gaps using a different model of structure. Check the known samples as well as the appearance inside the gap.', 'steps'),
  tomography:G('Increase the number of sampled Fourier lines and compare reconstruction methods.', 'Undersampling leaves unobserved directions. A prior can reduce artifacts, but it does not turn missing measurements into observed data.', 'lines'),
  sparse:G('Vary λ and compare the recovered signal with the reference.', 'Too little regularization retains noise; too much removes real features. Iteration speed and reconstruction quality answer different questions.', 'lam'),
  compressed:G('Increase the measurement count and watch the image reconstruction become better constrained.', 'The adjoint is a backprojection. Sparse reconstruction uses the wavelet prior to estimate information missing from the measurements.', 'P'),
  lq:G('Compare q below, at, and above 1, then change the scale of the set.', 'At q = 1 the boundary has corners on the axes. For q < 1 the set is nonconvex; q = 2 gives the Euclidean ball.', 'q'),
  certificates:G('Change the support size or measurement count and check the off-support certificate values.', 'Interpolating the support signs is only part of the condition. A strict off-support bound below 1 certifies the relevant nondegeneracy.', 's'),
  random:G('Vary the problem dimensions and compare a finite random experiment with its theoretical reference.', 'A finite random sample fluctuates. An asymptotic law or a sampled bound is not an exact description of every realization.', 'P'),
  pca:G('Keep PC1 fixed and choose another component for the second displayed coordinate.', 'PCA maximizes captured variance without using class labels. Good visual separation in a projection does not itself establish classification accuracy.', 'component'),
  clustering:G('Start at iteration 0 and follow the assignments, centroids, and energy together.', 'Each Lloyd update lowers or preserves the quantization energy. It can converge to a local minimum, and changing K changes the problem.', 'iteration'),
  learning:G('Change the model’s locality or regularization and inspect the resulting predictions.', 'Flexible decision boundaries can follow small features in the training sample. A visually convincing fit is not a validation score.', 'sigma'),
  stochastic:G('Compare optimization methods at the same number of equivalent data passes.', 'A batch gradient uses the full dataset; a stochastic update uses one observation. Compare computational budgets as well as convergence curves.', 'factor'),
};

const views={
  'examples-2':G('Move through the spectral channels and locate the selected wavelength on the spectrum.', 'These channels are synthesized from the RGB flower. They illustrate spectral sampling rather than a measured hyperspectral acquisition.', 'channel'),
  'filtering-1d':G('Begin with width zero, then increase the Gaussian width and inspect the Fourier magnitudes.', 'The filter smooths rapid variations by attenuating high frequencies. The periodic boundary convention joins the ends of the signal.', 'sigma'),
  '2d-fourier-extension':G('Change one frequency at a time and observe how the stripes rotate and become denser.', 'The wave vector (k, ℓ) is perpendicular to the constant-phase stripes. Real and imaginary parts differ by a quarter period.', 'k'),
  'fourier-2d':G('Increase the boundary taper and compare both spectra on their shared display scale.', 'Tapering suppresses the jump at the periodic boundary, while also modifying the image near its edges.', 'taper'),
  'the-concave-function-underlying-entropy':G('Move the probability dots along the entropy-contribution curve.', 'Each probability contributes −p log₂ p bits. The total entropy is the sum over all symbols, not the value at one dot.', 'skew'),
  'huff-code':G('Change the probabilities and follow the successive merges of the two lightest clusters.', 'The merge tree produces integer code lengths. Highly likely symbols should receive shorter words.', 'skew'),
  'codage-differences':G('Change the quantization alphabet and compare the signal histogram with its first differences.', 'Small differences can have a concentrated distribution even when the original intensities do not. Decoding also needs the first sample.', 'K'),
  'haard-scaling-wav':G('Change the dyadic level and compare the scaling function with the wavelet.', 'Dilation changes both the support width and amplitude. The normalization preserves the L² norm; here j is minus the displayed level.', 'level'),
  'wavelets':G('Compare the support and oscillations of Haar and Daubechies wavelets.', 'Increasing the level contracts these basis functions. This is a dilation parameter, rather than a count of image pixels.', 'level'),
  'inversion-transpose':G('Change the wavelet and compare the analysis matrix with its transpose.', 'The product WᵀW stays equal to the identity because this finite periodic transform is orthogonal.', 'wavelet'),
  'haar-analysis-synthesis-matrix':G('Inspect the analysis and synthesis matrices for different orthogonal wavelets.', 'Columns of W are transforms of coordinate vectors; Wᵀ performs the inverse transform.', 'wavelet'),
  'wavelets-2d':G('Compare the three tensor-product detail atoms for Haar and Daubechies wavelets.', 'Combining a scaling function with a wavelet gives one directional detail; combining two wavelets gives the diagonal detail.', 'wavelet'),
  'wavalgo-2d-step-fwd':G('Change the image or wavelet and follow a single two-dimensional analysis step.', 'Filtering and decimation along each coordinate separate the approximation from three directional detail bands.', 'wavelet'),
  'wavalgo-2d-step-bwd':G('Change the wavelet and inspect the bands used for reconstruction.', 'Synthesis recombines approximation and directional details. Perfect reconstruction requires the matching analysis and synthesis filters.', 'wavelet'),
  'image-model-cartoon':G('Add opaque leaves and follow how occlusions create new boundaries.', 'Only the visible parts of earlier leaves remain. This produces sharp edges between simple regions.', 'count'),
  'approximation-class':G('Increase the denoising strength in the bounded-variation example.', 'The BV image is obtained by TV denoising white noise. TV encourages small total variation while allowing jump discontinuities.', 'strength'),
  'signal-model-sobolev':G('Smooth the same white noise progressively and watch its spatial correlation grow.', 'Display variance is normalized: compare correlation length and structure rather than overall contrast.', 'strength'),
  'triangulation-anisoiso':G('Add boundary vertices and inspect the enlarged view of the same mesh.', 'Check the direction of the long edges near the jump. The enlargement uses the same P1 interpolant and triangles as the overview.', 'angular'),
  'triangulation-peppers':G('Change the mesh resolution and compare P1 coding with JPEG-2000 at equal transmitted size.', 'The mesh budget includes positions, nodal values, triangle indices, and a header. This is an explicit codec comparison, not an optimal rate theorem.', 'angular'),
  'curvelet-discretization':G('Rotate the frequency direction and compare the translation grid with the spatial atom.', 'The larger grid spacing lies along the atom’s tangent. The denser normal spacing resolves its rapid oscillations.', 'angle'),
  'curvelets-vs-wavelets':G('Rotate the elongated atoms away from the tangent and inspect the mismatch.', 'The ellipses represent localization. Tangential alignment lets elongated atoms cover a smooth edge efficiently.', 'angle'),
  'curvelet-denoising':G('Compare directional and wavelet denoising at the same threshold multiplier.', 'Different transforms capture different structures. A transform that is well aligned with edges need not win on every image.', 'factor'),
  'compression':G('Increase the quantization step and follow the loss of small wavelet coefficients.', 'The dead zone removes small coefficients; midpoint reconstruction quantizes those that remain. Count surviving coefficients separately from encoded bit rate.', 'rate'),
  'noise-example-1d':G('Increase σ and compare the noisy signal with the same clean reference.', 'The normalized noise sample is held fixed, so the control changes its amplitude rather than drawing a different disturbance.', 'sigma'),
  'noise-examples':G('Increase the noise level and compare three observation models on the same image.', 'The three variances match at clean intensity one. Poisson and multiplicative noise remain signal-dependent at other intensities.', 'sigma'),
  'sparse-signal':G('Increase the noise and look for the few nonzero signal coefficients.', 'A sparse signal concentrates its energy on a small support. Thresholding will try to distinguish those spikes from noise.', 'sigma'),
  'cycle-spin-principle':G('Translate the input and compare orthogonal denoising with cycle spinning.', 'Decimation makes an orthogonal wavelet denoiser sensitive to shifts. Averaging shifted estimates reduces that sensitivity.', 'shift'),
  'semisoft':G('Change the threshold, then adjust the width of the transition region.', 'Semisoft shrinkage interpolates between zeroing small values and keeping large values. The two thresholds are T and μT.', 'mu'),
  'semisoft-optimal':G('Compare semisoft transition widths at the same threshold.', 'An optimum in measured SNR depends on this particular signal and noise realization.', 'mu'),
  'poisson-distrib':G('Increase the Poisson mean and inspect the shift and spread of the mass function.', 'The mean and variance both equal the photon count. Relative fluctuations therefore scale as the inverse square root of the count.', 'photons'),
  'grad-discrete':G('Increase the arrow gain and inspect forward differences on the pixel grid.', 'Arrow direction follows (dx, dy). The image coordinate y increases downward; gain changes the drawing, not the derivative.', 'gain'),
  'discrete-operators-lapl':G('Vary ε and compare the TV force with the unchanged Laplacian.', 'The TV force normalizes the gradient by √(|∇f|² + ε²). Smaller ε increases its sensitivity near small gradients; the Laplacian does not depend on ε.', 'eps'),
  'discrete-operators-tveps':G('Change ε and inspect how the regularized TV operator responds near small gradients.', 'A positive ε smooths the singularity at zero gradient. It affects the operator itself, unlike a display-contrast control.', 'eps'),
  'tomo-principle':G('Rotate the projection direction and follow the corresponding one-dimensional profile.', 'Each profile value integrates the phantom along a line. It is a Radon projection, not a row of image pixels.', 'angle'),
  'sparsity-wav':G('Increase the selection threshold and examine which wavelet coefficients survive.', 'A coefficient map can contain many small values and relatively few large ones. Sparsity depends on both the image and its representation.', 'factor'),
  'spikes-ista-decay':G('Compare ISTA and FISTA at the same λ, using the convergence curve.', 'Both target the same regularized objective. Faster objective decay does not imply that a different λ would recover the signal better.', 'lam'),
  'deconvolution-2d-snr':G('Sweep λ and locate the largest reconstruction SNR.', 'The best SNR here uses the known clean image. It is an oracle comparison, not a parameter-selection rule for unknown data.', 'lam'),
  'optim-geom':G('Compare the contact points of scaled ℓ¹, ℓ², and selected-q balls with the affine line.', 'Minimizing the gauge amounts to finding the smallest ball that touches the constraint. The ℓ¹ contact lies on an axis in this example.', 'q'),
  'homotopy':G('Change the support size and follow how the coefficients enter along the λ path.', 'For large λ all coefficients vanish. Reading the path toward smaller λ reveals how the active set grows or changes.', 's'),
  'marcenko-pastur':G('Increase the matrix size at fixed aspect ratio and compare the histogram with the limiting density.', 'Finite-sample eigenvalues can lie outside the limiting support; all of them remain in the histogram.', 'P'),
  'marcenko-pastur-laws':G('Vary s/P and watch the limiting spectrum spread toward zero.', 'The lower edge is (1 − √β)² and the upper edge is (1 + √β)². These are asymptotic support edges.', 'beta'),
  'max-gaussian':G('Increase N and compare the sampled maxima with the √(2 log N) scale.', 'This experiment takes the maximum of the Gaussian variables without absolute values. The reference describes leading-order growth.', 'N'),
  'max-gaussian-mean':G('Increase N and compare the maximum absolute noise coefficient with its leading scale.', 'Even when every coefficient has unit variance, the largest absolute value grows with the number of coefficients.', 'N'),
  'phase-trans':G('Change the measurement count and compare recovery frequency with the sufficient criteria.', 'A criterion may fail even when recovery succeeds. Only 32 random trials are used at each setting.', 'P'),
  'rip-const':G('Increase the selected support size and inspect the sampled RIP bound.', 'The plotted values are lower bounds obtained from sampled supports. Small plotted values do not certify the worst-case RIP constant.', 's'),
  'digits-classes':G('Choose a different PCA plane and inspect the mixture of class probabilities.', 'The classifier acts in the full feature space. Colors show predictions on an affine PCA slice; the reported held-out accuracy is measured in the full space.', 'component'),
  'hist-classif':G('Increase the number of neighbors and inspect the smoother voting regions.', 'The displayed value is the fraction of neighbors in class 1. It is an empirical vote, not a calibrated probability guarantee.', 'k'),
  'losses':G('Move the negative margin through zero and compare both logistic normalizations.', 'Dividing logistic loss by log 2 makes it touch the binary loss at zero and upper-bound it everywhere. The unscaled logistic is the likelihood convention.', 'margin'),
  'separation-influ':G('Separate the class centers while keeping their internal displacements fixed.', 'C is an inverse penalty: larger C weakens quadratic regularization. Training accuracy alone does not measure generalization.', 'separation'),
  'kernel':G('Compare narrow and wide kernels at fixed λ, then vary regularization.', 'A narrow kernel can follow individual observations; a wide kernel couples more distant observations. The noise sample stays fixed.', 'sigma'),
  'classes-kernel':G('Change the Gaussian-feature width and inspect the nonlinear class boundary.', 'The model penalizes its finite feature coefficients. Its flexibility depends on both the feature width and λ.', 'sigma'),
  'sgd-traject':G('Increase the constant step size and compare the long-run fluctuations.', 'The stationary variance is τσ²/(2 − τ), for 0 < τ < 2. A larger stable step can converge faster initially while leaving more variance.', 'tau'),
};

const help={
  image:'Use the same image when comparing parameter settings.',
  signal:'Use the same signal when comparing methods or scales.',
  wavelet:'Daubechies N has N vanishing moments and 2N filter taps; Haar has one vanishing moment.',
  basis:'The budget counts real orthonormal coefficients in every basis.',
  level:'Number of dyadic analysis steps. More steps give a coarser approximation.',
  M:'Number of retained real coefficients, not the number of bytes in a file.',
  K:'',
  lam:'Larger λ gives the regularization term more weight.',
  C:'Larger C weakens the quadratic penalty.',
  sigma:'',
  photons:'Expected photon count at unit clean intensity.',
  looks:'The multiplicative noise variance is 1/K.',
  factor:'The threshold is T = (selected multiplier) × σ.',
  gain:'Display amplification only; the numerical derivative is unchanged.',
  beta:'The number of columns divided by the number of rows.',
  P:'Number of measurements or matrix rows.',
  component:'The second coordinate of the plot; the first coordinate remains PC1.',
  iteration:'0 is the initialization, before any centroid update.',
  steps:'More iterations need not remove the bias introduced by the chosen prior.',
  mu:'The upper semisoft threshold is μ times the lower threshold T.',
  q:'For q < 1 this is a nonconvex gauge, not a norm.',
  epochs:'One full data pass costs n individual gradients or one batch gradient.',
};

export function displayValue(value){return ({haar:'Haar',db2:'Daubechies 2',db4:'Daubechies 4'})[value]??String(value);}
export function controlHelp(experiment,control){
  const {family,view}=experiment;
  if(control.id==='level'&&['haard-scaling-wav','wavelets'].includes(view))return 'Dilation level L, with j = −L: increasing L contracts the support and multiplies amplitudes by 2^(L/2).';
  if(control.id==='factor'&&family==='stochastic')return 'Multiplier of the step schedule stated below; larger steps are not always faster.';
  if(control.id==='factor'&&view==='sparsity-wav')return 'Threshold for selecting the displayed wavelet coefficients.';
  if(control.id==='rate')return view==='compression'?'Larger T coarsens the quantization and removes more coefficients.':'Requested encoded bits per image pixel; the actual rate includes file headers.';
  if(control.id==='strength'&&family==='pde')return 'Heat time and Sobolev λ equal eight times this value; TV time and TV λ equal this value.';
  if(control.id==='radius'&&view==='optim-geom')return 'The value b in x₁ + 2x₂ = b; each ball is scaled to its first contact.';
  if(control.id==='K')return family==='coding'?'Number of source symbols.':family==='quantization'?'Number of uniformly spaced reconstruction levels.':'Number of clusters; changing K uses a new fixed initialization.';
  return help[control.id]||'';
}

export function getGuide(experiment){
  const guide={...families[experiment.family],...views[experiment.view]};
  if(!guide.prompt)throw new Error(`Missing exploration guide for ${experiment.id}`);
  const control=experiment.controls.find(c=>c.id===guide.focus)||experiment.controls[0];
  let presets=[
    {label:displayValue(control.values[0]),values:{[control.id]:control.values[0]}},
    {label:displayValue(control.values.at(-1)),values:{[control.id]:control.values.at(-1)}},
  ];
  if(control.id==='q')presets=[.5,1,2].map(q=>({label:`q = ${q}`,values:{q}}));
  else if(control.id==='margin')presets=[-1,0,1].map(margin=>({label:margin<0?'Correct side':margin===0?'At the boundary':'Wrong side',values:{margin}}));
  else if(control.id==='iteration')presets=[0,3,40].map(iteration=>({label:iteration===0?'Initialization':`Iteration ${iteration}`,values:{iteration}}));
  else if(control.id==='skew')presets=[{label:'Uniform',values:{skew:0}},{label:'Concentrated',values:{skew:4}}];
  else if(control.type==='range')presets.forEach(p=>p.label=`${control.label.replace(/\([^)]*\)/g,'').trim()}: ${p.label}`);
  return {...guide,presets};
}
