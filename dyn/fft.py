#! python
import numpy as np

fp = open("ani_s1.dat2", "r")
y = []
for i in xrange(601):
    line = fp.readline()
    rec = line.split()
    y.append(float(rec[1]))
fp.close()

y = np.array(y)
yf = np.fft.fft(y)/len(y)
for i in xrange(601):
    if np.abs(yf[i]) > 0.01:
        print i, yf[i]