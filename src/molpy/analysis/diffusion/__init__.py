from molpy.analysis import ComputeContext, Compute
import numpy as np
from molpy.core.logger import get_logger
from molpy.core.config import get_config

logger = get_logger(__name__)
config = get_config()

# Use fastest available fft library
try:
    import pyfftw

    logger.info("Using PyFFTW for FFTs")
    
    pyfftw.config.NUM_THREADS = min(1, config.n_threads)
    logger.info(f"Setting number of threads to {config.n_threads}")

    # Note that currently these functions are defined to match only the parts
    # of the numpy/scipy API that are actually used below. There is no promise
    # that other aspects of the API will be preserved.
    def fft(x, n, axis):
        a = pyfftw.empty_aligned(x.shape, "complex64")
        a[:] = x
        fft_object = pyfftw.builders.fft(a, n=n, axis=axis)
        return fft_object()

    def ifft(x, axis):
        a = pyfftw.empty_aligned(x.shape, "complex64")
        a[:] = x
        fft_object = pyfftw.builders.ifft(a, axis=axis)
        return fft_object()
except ImportError:
    try:
        from scipy.fftpack import fft, ifft

        logger.info("Using SciPy's fftpack for FFTs")
    except ImportError:
        from numpy.fft import fft, ifft

        logger.info("Using NumPy for FFTs")


def _autocorrelation(x):
    r"""Compute the autocorrelation of a sequence"""
    N = x.shape[0]
    F = fft(x, n=2 * N, axis=0)
    PSD = F * F.conjugate()
    res = ifft(PSD, axis=0)
    res = (res[:N]).real
    n = np.arange(1, N + 1)[::-1]  # N to 1
    return res / n[:, np.newaxis]


class DirectMSD(Compute):
    """Direct Mean Square Displacement (MSD) calculation."""

    def compute(self, context: ComputeContext) -> ComputeContext:

        positions = context.frame["atoms"]["xyz"]
        _msd_result = []

        _msd_result.append(
            np.linalg.norm(positions - positions[[0], :, :], axis=-1) ** 2
        )

        context.result[f"{self.name}_msd"] = np.array(_msd_result)

        return context
    
    
class WindowedMSD(Compute):
    """Windowed Mean Square Displacement (MSD) calculation."""

    def __init__(self, name: str, window_size: int = 1):
        super().__init__(name)
        self.window_size = window_size

    def compute(self, context: ComputeContext) -> ComputeContext:

        _msd_result = []

        positions = context.frame["atoms"]["xyz"]

        # First compute the first term r^2(k+m) - r^2(k)
        N = positions.shape[0]
        D = np.square(positions).sum(axis=2)
        D = np.append(D, np.zeros(positions.shape[:2]), axis=0)
        Q = 2 * D.sum(axis=0)
        S1 = np.zeros(positions.shape[:2])
        for m in range(N):
            Q -= D[m - 1, :] + D[N - m, :]
            S1[m, :] = Q / (N - m)

        # The second term can be computed via autocorrelation
        corrs = []
        for i in range(positions.shape[2]):
            corrs.append(_autocorrelation(positions[:, :, i]))
        S2 = np.sum(corrs, axis=0)

        _msd_result.append(S1 - 2 * S2)

        context.result[f"{self.name}_msd"] = np.array(_msd_result)

        return context