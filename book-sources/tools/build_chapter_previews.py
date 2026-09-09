#!/usr/bin/env python3
"""Build compact illustrations for the Foundations book's chapter catalogue."""
import argparse
import hashlib
import html
import json
from pathlib import Path

import pypdfium2 as pdfium

ROOT = Path(__file__).resolve().parents[1]

# Current figure IDs come from the published reader's figure catalogue. The
# asset index selects a complete panel when the book uses separate panel PDFs.
SELECTIONS = [
    ('shannon', 'shannon--aliasing', 0, 'Sampling and aliasing',
     'Two sine waves passing through the same sampled points'),
    ('fourier', 'fourier--fourier-2d', 0, 'Image spectra',
     'Images and their two-dimensional Fourier spectra'),
    ('shannon-coding', 'shannon-coding--tree-prefix', 0, 'Prefix codes',
     'Binary trees illustrating the structure of a prefix code'),
    ('wavelets', 'wavelets--wavcoefs-2d', 0, 'Wavelet coefficients',
     'Image wavelet coefficients arranged by scale and orientation'),
    ('approximation', 'approximation--triangulation-image', 0, 'Adaptive approximation',
     'A piecewise-smooth image approximated on triangles aligned with its boundary'),
    ('compression', 'compression--compression-jpeg', 0, 'Image compression',
     'JPEG and JPEG-2000 reconstructions of flower and mandrill images'),
    ('denoising', 'denoising--denoising-clean', 0, 'Removing noise',
     'A clean flower image, its noisy observation, and a denoised reconstruction'),
    ('variational-priors', 'variational-priors--discrete-operators', 0, 'Image gradients',
     'Horizontal and vertical image derivatives encoded in red and blue'),
    ('inverse-problems', 'inverse-problems--tomo-principle', 0, 'Tomographic projections',
     'A Shepp-Logan phantom and one of its one-dimensional Radon projections'),
    ('sparse-regularization', 'sparse-regularization--sparsity-lq', 0, 'The geometry of sparsity',
     'Norm balls illustrating the geometry of sparse regularization'),
    ('sparse-theory', 'sparse-theory--certif-convol', 0, 'Recovery certificates',
     'Dual certificate curves for recovering spikes from convolution measurements'),
    ('compressed-sensing', 'compressed-sensing--single-pixel', 0, 'Fewer measurements',
     'Compressed sensing reconstructions with increasing numbers of measurements'),
    ('machine-learning', 'machine-learning--pca', 0, 'Principal components',
     'A point cloud, its principal directions, and a lower-dimensional representation'),
    ('optim-ml-smooth', 'optim-ml-smooth--gradesc', 0, 'Gradient descent',
     'Gradient descent trajectories with different step sizes'),
    ('optim-ml-advanced', 'optim-ml-advanced--compariso-sgd', 0, 'Stochastic optimization',
     'Convergence curves for stochastic and variance-reduced gradient methods'),
    ('perceptrons', None, 0, 'One hidden layer',
     'A shallow network with two inputs, four nonlinear hidden neurons, and a linear output'),
    ('deep-learning', 'deep-learning--deep-fc-cnn', 1, 'Convolutional networks',
     'Local connections between successive layers of a convolutional neural network'),
    ('convex-analysis', 'convex-analysis--convex-analysis-subdifferential', 0, 'Supporting hyperplanes',
     'Subgradients and supporting lines to a convex function'),
    ('optim-nonsmooth', 'optim-nonsmooth--prox-proj', 0, 'Proximity and projection',
     'The geometry of a proximal map and projection onto a convex set'),
]


def shallow_network():
    """Illustrate the chapter’s one-hidden-layer model, sum_k u_k sigma(<v_k,x> + b_k)."""
    elements = []
    inputs, hidden = [120, 240], [78, 146, 214, 282]
    for y in inputs:
        for z in hidden:
            elements.append(f'<path d="M 85 {y} L 278 {z}"/>')
    for y in hidden:
        elements.append(f'<path d="M 322 {y} L 538 180"/>')
    edges = '\n'.join(elements)
    nodes = ''.join(f'<circle cx="300" cy="{y}" r="21" fill="#e1eeeb"/>'
                    f'<text x="300" y="{y+7}">σ</text>' for y in hidden)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="640" height="340" viewBox="0 0 640 340" role="img" aria-labelledby="title desc">
<title id="title">A neural network with one hidden layer</title>
<desc id="desc">Two input coordinates feed four nonlinear neurons. A linear output combines their activations with weights u_k.</desc>
<rect width="640" height="340" fill="white"/>
<defs><marker id="arrow" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M 0 0 L 6 3 L 0 6" fill="#8ca7a1"/></marker></defs>
<text x="320" y="29" text-anchor="middle" fill="#20464e" font-family="Georgia,serif" font-size="22">f<tspan baseline-shift="sub" font-size="15">θ</tspan>(x) = Σₖ uₖ σ(vₖ · x + bₖ)</text>
<g fill="none" stroke="#8ca7a1" stroke-width="1.7" marker-end="url(#arrow)">{edges}</g>
<g stroke="#24737a" stroke-width="2" font-family="Georgia,serif" font-size="24" text-anchor="middle">
<circle cx="64" cy="120" r="21" fill="#f0f5f3"/><circle cx="64" cy="240" r="21" fill="#f0f5f3"/>
<circle cx="560" cy="180" r="21" fill="#f5e8ce"/>
<g fill="#20464e" stroke="none">{nodes}<text x="64" y="127">x₁</text><text x="64" y="247">x₂</text><text x="560" y="187">Σ</text></g>
</g>
<g fill="#5f7778" font-family="Arial,sans-serif" font-size="17" text-anchor="middle"><text x="64" y="330">Input</text><text x="300" y="330">Hidden layer</text><text x="560" y="330">Output</text></g>
</svg>'''


def build(root, site):
    catalogue = json.loads((root.parent/'book/interactive/live/figures.json').read_text())
    figures = {figure['id']: figure for figure in catalogue}
    output = site/'img/previews'
    output.mkdir(parents=True, exist_ok=True)
    metadata = {}
    for number, (slug, figure_id, panel, label, alt) in enumerate(SELECTIONS, 1):
        item = {'number': number, 'label': label, 'alt': alt, 'figure_id': figure_id}
        if figure_id:
            figure = figures[figure_id]
            assert figure['chapter'] == slug
            source = root/figure['current_assets'][panel]
            with pdfium.PdfDocument(str(source)) as document:
                page = document[0]
                bitmap = page.render(scale=min(1000/page.get_width(), 700/page.get_height()))
                preview = bitmap.to_pil().convert('RGB')
                bitmap.close()
                page.close()
            preview.thumbnail((760, 480))
            name = f'foundations-{slug}.webp'
            preview.save(output/name, 'WEBP', quality=90, method=6)
            anchor = figure['book_label']
            item.update(source=source.relative_to(root).as_posix(),
                        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                        figure_number=figure['figure_number'], width=preview.width, height=preview.height)
        else:
            name = f'foundations-{slug}.svg'
            (output/name).write_text(shallow_network())
            anchor = 'eq-2layers-nn'
            item.update(source='chapters-tex/perceptrons.tex', width=640, height=340,
                        source_sha256=hashlib.sha256((root/'chapters-tex/perceptrons.tex').read_bytes()).hexdigest(),
                        note='Website illustration of the chapter’s one-hidden-layer model; no numbered book figure.')
        reader = root.parent/'book/interactive'/slug/'index.html'
        assert f'id="{html.escape(anchor, quote=True)}"' in reader.read_text(), f'Missing reader anchor: {slug}#{anchor}'
        item.update(image='/img/previews/'+name, online=f'/book/interactive/{slug}/#{anchor}')
        metadata[slug] = item
    target = site/'_data/foundations_chapter_previews.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(metadata, ensure_ascii=False, indent=2)+'\n')
    total = sum((site/item['image'].lstrip('/')).stat().st_size for item in metadata.values())
    print(f'Built {len(metadata)} chapter previews ({total/1024:.0f} KiB total).')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--book-root', type=Path, default=ROOT)
    parser.add_argument('--site-root', type=Path, help='Defaults to the parent website repository.')
    args = parser.parse_args()
    root = args.book_root.resolve()
    build(root, args.site_root.resolve() if args.site_root else root.parent)


if __name__ == '__main__':
    main()
