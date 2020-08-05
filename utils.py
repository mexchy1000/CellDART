#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from keras.utils import to_categorical

def random_mix(Xs, ys, nmix=5, n_samples=10000):
    nclss=len(set(ys))
    Xs_new, ys_new =[], []
    ys_ = to_categorical(ys)
    for i in range(n_samples):
        yy = np.zeros(nclss)
        fraction = np.random.rand(nmix)
        fraction = fraction/np.sum(fraction)
        fraction = np.reshape(fraction, (nmix,1))
        randindex = np.random.randint(len(Xs), size=nmix)
        ymix = ys_[randindex]
        yy = np.sum(ymix*np.reshape(fraction, (nmix,1)), axis=0)
        XX = Xs[randindex] * fraction
        XX_ = np.sum(XX, axis=0)
        ys_new.append(yy)
        Xs_new.append(XX_)
    Xs_new = np.asarray(Xs_new)
    ys_new = np.asarray(ys_new)
    return Xs_new, ys_new