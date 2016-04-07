import csv
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt

fitdata = genfromtxt('pi0peak_fit_Run13pp510ERT.txt', delimiter=' ')

print fitdata

plt.plot( fitdata[:,1] )

plt.show()
