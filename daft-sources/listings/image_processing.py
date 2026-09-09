"""Image filtering, valid template matching, and Fourier resampling."""
import numpy as np
from scipy.signal import fftconvolve
from PIL import Image


def gaussian_kernel(sigma, radius=None):
    if sigma <= 0:
        raise ValueError('Gaussian width must be positive.')
    radius = int(np.ceil(4*sigma)) if radius is None else int(radius)
    if radius < 0:
        raise ValueError('Kernel radius must be nonnegative.')
    k = np.arange(-radius, radius+1)
    g = np.exp(-(k[:, None]**2+k[None, :]**2)/(2*sigma**2))
    return g/g.sum()


def smooth_image(image, sigma):
    kernel = gaussian_kernel(sigma)
    radius = kernel.shape[0]//2
    padded = np.pad(np.asarray(image, dtype=float), radius, mode='reflect')
    return fftconvolve(padded, kernel, mode='valid')


def apply_gaussian(path, sigma=1):
    with Image.open(path) as source:
        image = np.asarray(source.convert('L'), dtype=float)/255
    return image, smooth_image(image, sigma)


def sliding_sum(image, width, power=1):
    image = np.asarray(image, dtype=float)
    if image.ndim != 2 or not 1 <= width <= min(image.shape):
        raise ValueError('Window must fit within a two-dimensional image.')
    table = np.pad(np.cumsum(np.cumsum(image**power, axis=0), axis=1), ((1, 0), (1, 0)))
    return (table[width:, width:]-table[:-width, width:]
            - table[width:, :-width]+table[:-width, :-width])


def normalized_correlation(image, template):
    f, g = np.asarray(image, float), np.asarray(template, float)
    if g.ndim != 2 or g.shape[0] != g.shape[1] or g.shape[0] > min(f.shape):
        raise ValueError('Use a square template that fits inside the image.')
    f, g = f-f.mean(), g-g.mean()
    energy = np.sum(g*g)
    if energy == 0:
        raise ValueError('Normalized correlation is undefined for a constant template.')
    width = g.shape[0]
    numerator = fftconvolve(f, g[::-1, ::-1], mode='valid')
    s1, s2 = sliding_sum(f, width), sliding_sum(f, width, 2)
    variance = np.maximum(s2-s1*s1/g.size, 0)
    valid = variance > 32*np.finfo(float).eps*np.maximum(s2, 0)
    result = np.full(variance.shape, np.nan)
    np.divide(numerator, np.sqrt(variance*energy), out=result, where=valid)
    return np.clip(result, -1, 1)


def fft_translate(image, shift):
    image = np.asarray(image)
    if image.ndim != 2 or any(n % 2 == 0 for n in image.shape):
        raise ValueError('Use odd sizes to avoid a Nyquist interpolation ambiguity.')
    k = np.fft.fftfreq(image.shape[0])[:, None]
    l = np.fft.fftfreq(image.shape[1])[None, :]
    phase = np.exp(-2j*np.pi*(k*shift[0]+l*shift[1]))
    y = np.fft.ifft2(np.fft.fft2(image)*phase)
    return y.real if np.isrealobj(image) else y


def fft_shear(image, amount, axis=0):
    image = np.asarray(image)
    if image.ndim != 2 or any(n % 2 == 0 for n in image.shape) or axis not in (0, 1):
        raise ValueError('Use odd image sizes and axis 0 or 1.')
    n, m = image.shape[axis], image.shape[1-axis]
    phase = np.exp(-2j*np.pi*np.fft.fftfreq(n)[:, None]*amount*(np.arange(m)-(m-1)/2))
    if axis == 1:
        phase = phase.T
    y = np.fft.ifft(np.fft.fft(image, axis=axis)*phase, axis=axis)
    return y.real if np.isrealobj(image) else y


def fft_rotate(image, angle):
    if abs(angle) >= np.pi:
        raise ValueError('Split large rotations into angles of modulus less than pi.')
    a = -np.tan(angle/2)
    return fft_shear(fft_shear(fft_shear(image, a, 0), np.sin(angle), 1), a, 0)


def rotation_examples(image, angles=(0, np.pi/12, np.pi/6, np.pi/4)):
    image = np.asarray(image)
    padding = max(image.shape)//2+1
    padded = np.pad(image, padding, mode='constant', constant_values=1)
    return np.stack([fft_rotate(padded, angle) for angle in angles])


def lowpass_kernel(n):
    if n < 3 or n % 2 == 0:
        raise ValueError('Use an odd length >= 3.')
    frequency = abs(np.fft.fftfreq(n))
    response = (frequency < .25).astype(float)
    response[frequency == .25] = .5
    return np.fft.ifft(response).real


def smooth_lowpass_kernel(n, epsilon):
    if not 0 <= epsilon <= 1:
        raise ValueError('Transition parameter must lie in [0, 1].')
    if epsilon == 0:
        return lowpass_kernel(n)
    if n < 3 or n % 2 == 0:
        raise ValueError('Use an odd length >= 3.')
    frequency = abs(np.fft.fftfreq(n))
    t = np.clip((frequency-(1-epsilon)/4)/(epsilon/2), 0, 1)
    return np.fft.ifft((1+np.cos(np.pi*t))/2).real


def filter_response(kernel, length=2048):
    kernel = np.asarray(kernel)
    n = kernel.size
    if n % 2 == 0 or length < n:
        raise ValueError('Use an odd kernel length and a larger plotting grid.')
    index = np.rint(np.fft.fftfreq(n)*n).astype(int)
    padded = np.zeros(length, dtype=complex)
    padded[index % length] = kernel
    omega = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(length))
    return omega, np.fft.fftshift(np.fft.fft(padded))
