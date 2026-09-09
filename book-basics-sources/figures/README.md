# Illustration notes

The revised books retain the original wave cover and the original information,
image-processing, transport, and neural-network illustrations where their content
is appropriate. Unused assets have been removed from this working directory;
the unchanged original collection is preserved in `../legacy/`.

The larger original transport exports are used for legibility. Selected PDF
assets were compacted with Ghostscript without raster downsampling; vector
annotations were retained, page dimensions checked, and rendered output compared.
The file-size and rendering checks are recorded in
`../reports/figure-compaction.json`. This optimization is already applied to the
supplied figure files and is not required for a normal book build.

The updated sampling and inverse-problem figures use Maïtine Bergounioux's flower
photograph, inherited from the original image-processing chapter. The input is
`../data/flower.png`; `../figures-code/chapter-sampling/generate.py` computes the
wavelet approximations, blur, deterministic mask, and actual sparse reconstruction.
The PDF plots have embedded TrueType fonts and boxed axes. Error measurements
are computed before clipping intensities for display.

The English neural-network chapter uses the original English-labeled diagrams;
the French chapter uses the French-labeled versions. Transport and information
diagrams retain their original embedded labels, with revised captions in each
language. The metro illustration and its displayed assignment have been checked
against one another. Citations for the single-pixel camera and the application
figures remain in the corresponding chapter text and bibliography.
