import numpy as np


def minmod(x, y):
    """
    制限関数
    :param x,y: 制限される前のdelta二つ。どちらもsize(3,1)のndarray
    :return delta: 制限されたdelta
    """
    sgn = np.sign(x)
    return sgn * np.maximum(np.minimum(np.abs(x), sgn * y), 0.0)


def superbee():
    pass
