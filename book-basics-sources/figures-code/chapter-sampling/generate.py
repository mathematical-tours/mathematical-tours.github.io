"""Reproduce the flower sampling, approximation and inverse-problem figures.

The photograph is inherited from the original book (Maïtine Bergounioux's
image-processing material). Input intensities are normalized to [0, 1];
reconstructed intensities are clipped only when saving display images.
Run from any directory with Python, NumPy, SciPy, Pillow, Matplotlib, PyWavelets.
"""
from pathlib import Path
import json
import os
ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'build/matplotlib'))
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
import pywt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch

OUT = ROOT / 'figures/chapter-sampling'
OUT.mkdir(parents=True, exist_ok=True)
(ROOT / 'reports').mkdir(exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'pdf.fonttype': 42,
                     'axes.spines.top': True, 'axes.spines.right': True,
                     'axes.linewidth': .7, 'savefig.bbox': 'tight'})
flower = np.asarray(Image.open(ROOT / 'data/flower.png').convert('RGB'), float)
f = flower.mean(axis=2) / 255
assert f.shape == (256, 256)

def save_image(name, a):
    Image.fromarray(np.rint(255*np.clip(a, 0, 1)).astype('uint8')).save(OUT / name)

def analysis(a):
    return pywt.coeffs_to_array(pywt.wavedec2(a, 'db4', mode='periodization', level=5))[0]

_, slices = pywt.coeffs_to_array(pywt.wavedec2(f, 'db4', mode='periodization', level=5))
def synthesis(c):
    return pywt.waverec2(pywt.array_to_coeffs(c, slices, output_format='wavedec2'),
                        'db4', mode='periodization')

save_image('original.png', f)
coeff = analysis(f)
assert np.max(np.abs(synthesis(coeff)-f)) < 1e-12
assert abs(np.linalg.norm(coeff)-np.linalg.norm(f)) < 1e-10
order = np.argsort(-np.abs(coeff.ravel()), kind='stable')
approx_report = []
for ratio in (4, 8, 16):
    m = f.size // ratio
    c = np.zeros(f.size)
    c[order[:m]] = coeff.ravel()[order[:m]]
    a = synthesis(c.reshape(coeff.shape))
    save_image(f'approx-{ratio}.png', a)
    err = np.linalg.norm(a-f)
    assert abs(err-np.linalg.norm(c.reshape(coeff.shape)-coeff)) < 1e-10
    approx_report.append({'M': m, 'relative_l2_error': float(err/np.linalg.norm(f))})
assert all(approx_report[i]['relative_l2_error'] < approx_report[i+1]['relative_l2_error'] for i in range(2))

save_image('blurred.png', gaussian_filter(f, 2, mode='reflect'))
rng = np.random.default_rng(20260908)
mask = np.zeros(f.size)
p = round(.1*f.size)
mask[rng.choice(f.size, p, replace=False)] = 1
mask = mask.reshape(f.shape)
y = mask*f
save_image('masked.png', y)

# FISTA for 1/2 ||mask*(Psi*x-f)||^2 + lambda ||x||_1.
# The orthonormal wavelet transform and mask give gradient Lipschitz constant 1.
lam = .005
c = np.zeros_like(f)
z = c.copy()
t = 1.
for iteration in range(1, 8001):
    grad = analysis(mask*(synthesis(z)-f))
    u = z-grad
    next_c = np.sign(u)*np.maximum(np.abs(u)-lam, 0)
    next_t = (1+np.sqrt(1+4*t*t))/2
    z = next_c+(t-1)/next_t*(next_c-c)
    c, t = next_c, next_t
    if iteration % 50 == 0:
        u = c-analysis(mask*(synthesis(c)-f))
        prox = np.sign(u)*np.maximum(np.abs(u)-lam, 0)
        mapping = np.linalg.norm(c-prox)/max(1., np.linalg.norm(c))
        if mapping < 1e-7:
            break
restored = synthesis(c)
assert mapping < 2e-7
save_image('reconstructed.png', restored)
report = {'seed': 20260908, 'pixels': int(f.size), 'observed_pixels': p,
          'wavelet': 'db4 (four vanishing moments)', 'boundary': 'periodization',
          'levels': 5, 'approximations': approx_report, 'lambda': lam,
          'iterations': iteration, 'relative_proximal_gradient_residual': float(mapping),
          'objective': float(.5*np.linalg.norm(mask*(restored-f))**2+lam*np.abs(c).sum()),
          'psnr_db_before_display_clipping': float(-10*np.log10(np.mean((f-restored)**2)))}
(ROOT / 'reports/sampling-figures.json').write_text(json.dumps(report, indent=2)+'\n')

# Pixel sampling illustration with corresponding boxes and a nearest-neighbor zoom.
fig, axes = plt.subplots(1, 2, figsize=(6.1, 2.5), gridspec_kw={'width_ratios':[1, .85]})
for ax in axes:
    ax.set_xticks([]); ax.set_yticks([])
x0, y0, width = 124, 154, 16
axes[0].imshow(f, cmap='gray', vmin=0, vmax=1)
axes[0].add_patch(Rectangle((x0-.5, y0-.5), width, width, fill=False, ec='#c9333b', lw=1.6))
axes[1].imshow(f[y0:y0+width, x0:x0+width], cmap='gray', vmin=0, vmax=1, interpolation='nearest')
for spine in axes[1].spines.values():
    spine.set_edgecolor('#c9333b'); spine.set_linewidth(1.6)
fig.add_artist(ConnectionPatch((x0+width-.5, y0-.5), (0, 1), 'data', 'axes fraction',
                              axesA=axes[0], axesB=axes[1], color='#c9333b', lw=.6))
fig.tight_layout(w_pad=2)
fig.savefig(OUT/'discrete-image.pdf')
plt.close(fig)

# Continuous signal, discrete observations, and a magnified sampling interval.
def signal(t):
    return .65*np.sin(2*np.pi*3*t)+.3*np.cos(2*np.pi*7*t)
fig, axes = plt.subplots(1, 2, figsize=(6.1, 2.5), sharey=True)
for ax, interval in zip(axes, [(0,1),(.3,.5)]):
    tt = np.linspace(*interval, 1200)
    ts = np.arange(0, 1.0001, 1/40)
    ts = ts[(ts>=interval[0]) & (ts<=interval[1])]
    ax.plot(tt, signal(tt), color='#146c78', lw=1.2)
    ax.plot(ts, signal(ts), 'o', color='#cb3e42', ms=3)
    ax.set_xlim(interval); ax.set_ylim(-1.1,1.1); ax.set_xlabel('$t$')
    ax.grid(alpha=.15)
axes[0].set_ylabel('$f(t)$')
fig.tight_layout(w_pad=1)
fig.savefig(OUT/'discrete-signal.pdf')
plt.close(fig)
print(json.dumps(report, indent=2))
