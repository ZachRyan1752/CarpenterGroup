import tifffile
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

GainCalibration = "No" ## Nothing for yes, "No" for no

DataFile1 = tifffile.imread("UltraquietGainCalibration%sCorrection/ElectronReadoutNoise.tif" % GainCalibration)
DataFile2 = tifffile.imread("StandardGainCalibration%sCorrection/ElectronReadoutNoise.tif" % GainCalibration)
DataFile3 = tifffile.imread("FastGainCalibration%sCorrection/ElectronReadoutNoise.tif" % GainCalibration)

DataRaveled1 = DataFile1.ravel()
DataRaveled2 = DataFile2.ravel()
DataRaveled3 = DataFile3.ravel()

#plt.plot()
## General plot settings
bincount = 1000
histformat = "step"
ymax = 1000000
ymin = 1
ytickcount = 6
fig, ax = plt.subplots()

## Accepted colors: https://matplotlib.org/stable/gallery/color/named_colors.html
ax.hist(DataRaveled1, bins = bincount, histtype = histformat, color = "r", label = "Ultraquiet")
ax.hist(DataRaveled2, bins = bincount, histtype = histformat, color = "rebeccapurple", label = "Standard")
ax.hist(DataRaveled3, bins = bincount, histtype = histformat, color = "orange", label = "Fast")


## X-axis settings
ax.set(xlim = (0,10),xlabel = "Readout Noise (electrons)")
ax.locator_params(axis="x",nbins = 10)

## Y-axis settings
ax.set(ylim = (ymin, ymax), ylabel = "Pixel Number (Pixels)", yticks = np.linspace(ymin,ymax,ytickcount), yscale = 'log')

## Both axis settings
formatter = ticker.ScalarFormatter(useMathText = False)
formatter.set_scientific(False)
ax.yaxis.set_major_formatter(formatter)

ax.set_ylabel("Pixel Number (Pixels)", fontsize = 18)
ax.set_xlabel("Readout Noise (electrons)", fontsize = 18)


plt.subplots_adjust(left = 0.2)
plt.legend(loc = "upper right")
## Getting RMS data
hist1, bin_edges1 = np.histogram(DataRaveled1, bins = bincount)
hist2, bin_edges2 = np.histogram(DataRaveled2, bins = bincount)
hist3, bin_edges3 = np.histogram(DataRaveled3, bins = bincount)

bin_centers1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2
bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2
bin_centers3 = (bin_edges3[:-1] + bin_edges3[1:]) / 2


rms1 = np.sqrt(np.sum((bin_centers1**2) * hist1) / np.sum(hist1))
rms2 = np.sqrt(np.sum((bin_centers2**2) * hist2) / np.sum(hist2))
rms3 = np.sqrt(np.sum((bin_centers3**2) * hist3) / np.sum(hist3))


print("RMS 1: %s" % rms1)
print("RMS 2: %s" % rms2)
print("RMS 3: %s" % rms3)


## Displaying the plot
plt.show()