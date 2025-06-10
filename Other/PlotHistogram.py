import tifffile
import matplotlib.pyplot as plt
import numpy as np

InputFile = tifffile.imread("ElectronReadoutNoise.tif")

InputFileRaveled = np.sqrt(InputFile).ravel()

#plt.plot()
fig = plt.hist(InputFileRaveled, bins = 1000)
plt.yscale('log')
plt.xlim(0,10)
plt.ylim(1,1000000)
plt.xlabel("Readout Noise (electrons)")
plt.ylabel("Pixel Number (Pixels)")
plt.show()