import numpy as np


def minmod(x, y):
    """
    制限関数
    :param x,y: 制限される前のdelta二つ。どちらもsize(3,1)のndarray
    :return delta: 制限されたdelta
    """
    sgn = np.sign(x)
    return sgn * np.maximum(np.minimum(np.abs(x), sgn * y), 0.0)


def superbee(x, y):
    sgn = np.sign(x)
    return sgn * np.maximum(np.minimum(2 * np.abs(x), sgn * y), np.minimum(np.abs(x), 2 * sgn * y), 0.0)

def van_leer(x, y):
    sgn = np.sign(x)
    return sgn * np.maximum(np.minimum(2 * np.abs(x), sgn * y), np.minimum(np.abs(x), 2 * sgn * y), 0.0)


def van_albada(x, y):
    sgn = np.sign(x)
    return sgn * np.maximum(np.minimum(2 * np.abs(x), sgn * y), np.minimum(np.abs(x), 2 * sgn * y), 0.0)


def umist(x, y):
    sgn = np.sign(x)
    return sgn * np.maximum(np.minimum(2 * np.abs(x), sgn * y), np.minimum(np.abs(x), 2 * sgn * y), 0.0)
