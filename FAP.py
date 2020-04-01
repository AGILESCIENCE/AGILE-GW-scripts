# USAGE: python FAP.py [sigma] [t-T0] [t_bin]
# ex: python FAP.py 4.4 0.87 0.032

import os, sys
import numpy as np

N = 7

far = [[4.0, 4.33E-03, 2.06E-03, 9.16E-04, 6.27E-04],
[4.1, 7.96E-03, 3.73E-03, 1.65E-03, 1.13E-03],
[4.2, 1.09E-02, 5.07E-03, 2.23E-03, 1.53E-03],
[4.3, 1.33E-02, 6.13E-03, 2.69E-03, 1.84E-03],
[4.4, 1.53E-02, 6.98E-03, 3.04E-03, 2.08E-03],
[4.5, 1.69E-02, 7.65E-03, 3.31E-03, 2.27E-03],
[4.6, 1.82E-02, 8.19E-03, 3.52E-03, 2.41E-03],
[4.7, 1.92E-02, 8.61E-03, 3.68E-03, 2.53E-03],
[4.8, 2.01E-02, 8.95E-03, 3.81E-03, 2.62E-03],
[4.9, 2.08E-02, 9.21E-03, 3.91E-03, 2.69E-03],
[5.0, 2.13E-02, 9.42E-03, 3.99E-03, 2.75E-03],
[5.1, 2.18E-02, 9.59E-03, 4.05E-03, 2.79E-03],
[5.2, 2.22E-02, 9.72E-03, 4.09E-03, 2.82E-03],
[5.3, 2.25E-02, 9.82E-03, 4.13E-03, 2.84E-03],
[5.4, 2.27E-02, 9.90E-03, 4.15E-03, 2.86E-03],
[5.5, 2.29E-02, 9.97E-03, 4.17E-03, 2.88E-03],
[5.6, 2.31E-02, 1.00E-02, 4.19E-03, 2.89E-03],
[5.7, 2.32E-02, 1.01E-02, 4.20E-03, 2.89E-03],
[5.8, 2.33E-02, 1.01E-02, 4.21E-03, 2.90E-03],
[5.9, 2.34E-02, 1.01E-02, 4.22E-03, 2.91E-03],
[6.0, 2.35E-02, 1.01E-02, 4.22E-03, 2.91E-03],
[6.1, 2.36E-02, 1.02E-02, 4.23E-03, 2.91E-03],
[6.2, 2.36E-02, 1.02E-02, 4.23E-03, 2.91E-03],
[6.3, 2.37E-02, 1.02E-02, 4.23E-03, 2.92E-03],
[6.4, 2.37E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[6.5, 2.37E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[6.6, 2.38E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[6.7, 2.38E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[6.8, 2.38E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[6.9, 2.38E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.0, 2.38E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.1, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.2, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.3, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.4, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.5, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.6, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.7, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.8, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[7.9, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[8.0, 2.39E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[8.1, 2.40E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[8.2, 2.40E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[8.3, 2.40E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[8.4, 2.40E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[8.5, 2.40E-02, 1.02E-02, 4.24E-03, 2.92E-03],
[8.6, 2.40E-02, 1.03E-02, 4.24E-03, 2.92E-03],
[8.7, 2.40E-02, 1.03E-02, 4.24E-03, 2.92E-03],
[8.8, 2.40E-02, 1.03E-02, 4.24E-03, 2.92E-03],
[8.9, 2.40E-02, 1.03E-02, 4.24E-03, 2.92E-03],
[9.0, 2.40E-02, 1.03E-02, 4.24E-03, 2.92E-03]]

far = np.array(far)

sigma = float(sys.argv[1])

dt = float(sys.argv[2])

Dt = 50.0

t_bin = float(sys.argv[3])

FAR = 0

for i in range(51):
  if abs(float(far[i,0])-sigma) <= 0.5:
      if t_bin == 0.016:
        FAR = float(far[i,1])
      if t_bin == 0.032:
        FAR = float(far[i,2])
      if t_bin == 0.064:
        FAR = float(far[i,3])
      if t_bin == 0.128:
        FAR = float(far[i,4])

FAP = float(N * FAR * dt * 2 * np.log(float(Dt/t_bin)))

import os, sys, math
from scipy.special import erf
import numpy as np

ERF = []

for j in range(10000):
	ERF.append([j*0.001, 1-erf(abs(j*0.001)/math.sqrt(2)), abs(1-erf(abs(j*0.001)/math.sqrt(2))-FAP) ])

ERF = np.array(ERF)

print "FAP(%3.1f sigma) = (%d) * (%.2E Hz) * (%3.2f s) * [2 * ln(%3.2f s / %3.3f s)] = %3.3f = %3.2f sigma" % (sigma, N, FAR, dt, Dt, t_bin, FAP, ERF[int(ERF[:,2].argmin()),0])
