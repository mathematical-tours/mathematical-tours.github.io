# Inputs for regenerated figures

- `iris.npz`: the unmodified feature and label arrays returned by scikit-learn's bundled `load_iris()`. Measurements are in centimeters; figure code states when centering is applied. `iris-source.txt` preserves the dataset description and provenance supplied by the package.
- `digits.npz`: the unmodified 1797 by 64 optical handwritten-digit array and labels from scikit-learn's bundled `load_digits()`. Values range from 0 to 16; generators explicitly normalize when needed. `digits-source.txt` preserves the package's source and citation information.
- `parrot-cage.png`, `parrot-observation.png` and `parrot-known-mask.png`: historical inputs retained from the first reconstruction pass. They are no longer used by the generators. The cage-removal proposal now uses the supplied flower and an explicitly generated occlusion mask.

- `mandrill.png`: the 256 × 256 mandrill image extracted without resampling from the embedded raster in the historical mandrill wavelet PDF, at the author’s explicit request for Figures 5.8 and 6.5. Its companion approximation-error curve also uses this input. `mandrill-source.json` records the extraction and hash.
- `felix-1919.jpg`: a frame from *Feline Follies* (1919), credited to Pat Sullivan on [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Felix_1919.jpg), where it is identified as public domain. It supplies the author-requested real grayscale cartoon in Figure 5.7; `felix-1919-source.json` records the URL and hash.

Every generated numerical figure records the SHA-256 of each input it uses. Other photographic experiments use `data/flower.png`, with crop and resampling information recorded where relevant. Audio examples use `data/bird.wav`. No new reconstruction uses Lena. Historical comparison panels remain available solely to document the previous book figures.
