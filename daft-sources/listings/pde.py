"""Fourier solvers with explicit time, grid, and boundary conventions."""
import numpy as np


def heat_explicit(x, s, steps):
    x = np.asarray(x)
    y = np.array(x, dtype=np.result_type(x, float), copy=True)
    if y.ndim != 1 or y.size < 3 or s < 0 or steps < 0:
        raise ValueError('Use a vector of length >= 3 and nonnegative parameters.')
    for _ in range(steps):
        y += s*(np.roll(y, 1)+np.roll(y, -1)-2*y)
    return y


def heat_theta(x, s, theta, steps):
    x = np.asarray(x)
    if x.ndim != 1 or x.size < 3 or s < 0 or not 0 <= theta <= 1 or steps < 0:
        raise ValueError('Invalid grid, time step, or theta.')
    a = -4*s*np.sin(np.pi*np.fft.fftfreq(x.size))**2
    multiplier = (1+(1-theta)*a)/(1-theta*a)
    y = np.fft.ifft(multiplier**steps*np.fft.fft(x))
    return y.real if np.isrealobj(x) else y


def heat_theta_2d(x, s, theta, steps):
    x = np.asarray(x)
    if x.ndim != 2 or min(x.shape) < 3 or s < 0 or not 0 <= theta <= 1 or steps < 0:
        raise ValueError('Invalid grid, time step, or theta.')
    a = -4*s*(np.sin(np.pi*np.fft.fftfreq(x.shape[0]))[:, None]**2
              + np.sin(np.pi*np.fft.fftfreq(x.shape[1]))[None, :]**2)
    multiplier = (1+(1-theta)*a)/(1-theta*a)
    y = np.fft.ifft2(multiplier**steps*np.fft.fft2(x))
    return y.real if np.isrealobj(x) else y


def odd_extend(interior):
    y = np.asarray(interior)
    for axis in range(y.ndim):
        shape = list(y.shape); shape[axis] = 1
        zero = np.zeros(shape, dtype=y.dtype)
        y = np.concatenate([zero, y, zero, -np.flip(y, axis)], axis=axis)
    return y


def sine_transform(x, axis=0):
    x = np.asarray(x)
    n = x.shape[axis]+1
    shape = list(x.shape); shape[axis] = 1
    zero = np.zeros(shape, dtype=x.dtype)
    odd = np.concatenate([zero, x, zero, -np.flip(x, axis)], axis=axis)
    y = np.take(np.fft.fft(odd, axis=axis), np.arange(1, n), axis=axis)/(-2j)
    return y.real if np.isrealobj(x) else y


def sine_transform_2d(x):
    return sine_transform(sine_transform(x, axis=0), axis=1)


def poisson_zero_boundary(rhs, h):
    rhs = np.asarray(rhs)
    if h <= 0 or rhs.ndim < 1 or min(rhs.shape) == 0:
        raise ValueError('Use a nonempty interior grid and positive spacing.')
    extended = odd_extend(rhs)
    symbol = np.zeros(extended.shape)
    for axis, length in enumerate(extended.shape):
        shape = [1]*rhs.ndim; shape[axis] = length
        symbol -= 4*np.sin(np.pi*np.arange(length)/length).reshape(shape)**2/h**2
    zero = (0,)*rhs.ndim
    symbol[zero] = 1
    spectrum = np.fft.fftn(extended)/symbol
    spectrum[zero] = 0
    solution = np.fft.ifftn(spectrum)
    interior = tuple(slice(1, n+1) for n in rhs.shape)
    return solution[interior].real if np.isrealobj(rhs) else solution[interior]


def poisson_dirichlet(source, boundary, h):
    source, boundary = np.asarray(source), np.asarray(boundary)
    if source.shape != boundary.shape or min(source.shape) < 3:
        raise ValueError('Source and boundary arrays must share a grid with an interior.')
    center = (slice(1, -1),)*source.ndim
    result = np.array(boundary, dtype=np.result_type(source, boundary, float), copy=True)
    result[center] = 0
    rhs = np.array(source[center], dtype=result.dtype, copy=True)
    for axis in range(source.ndim):
        neighbors = np.roll(result, 1, axis)+np.roll(result, -1, axis)
        rhs -= neighbors[center]/h**2
    result[center] = poisson_zero_boundary(rhs, h)
    return result


def poisson_source(x, y):
    return (x*x+y*y)*np.exp(x*y)


def poisson_exact(x, y):
    return np.exp(x*y)


def poisson_left_boundary(y):
    return poisson_exact(0, y)


def poisson_example(n=32):
    grid = np.linspace(0, 1, n+1)
    X, Y = np.meshgrid(grid, grid, indexing='ij')
    exact = poisson_exact(X, Y)
    solution = poisson_dirichlet(poisson_source(X, Y), exact, 1/n)
    return grid, solution, exact


def heat_periodic(x, time, kappa=.5, period=1):
    x = np.asarray(x)
    if x.ndim != 1 or x.size == 0 or time < 0 or kappa < 0 or period <= 0:
        raise ValueError('Invalid heat equation data.')
    omega = 2*np.pi*np.fft.fftfreq(x.size, d=period/x.size)
    y = np.fft.ifft(np.exp(-kappa*time*omega**2)*np.fft.fft(x))
    return y.real if np.isrealobj(x) else y


def heat_initial(x):
    return ((np.asarray(x) >= .3) & (np.asarray(x) < .7)).astype(float)


def heat_example(n=256, times=(0, .001, .003, .01, .03, .1)):
    grid = np.arange(n)/n
    initial = heat_initial(grid)
    return grid, np.stack([heat_periodic(initial, t) for t in times])
