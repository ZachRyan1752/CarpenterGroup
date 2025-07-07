import tifffile
import numpy as np

Mode = "No"
#Region = [350:926,288:864]

Movie = tifffile.imread("R110SM_TrainingScope_30_05_2025_UltraQuiet_100ms(4).tif")
Gain = tifffile.imread("UltraquietGainCalibration%sCorrection/GainPixelSlope.tif" % Mode)

Movie = Movie[0:300,288:864,400:976]
Gain = Gain[288:864,400:976]

ImageOut = Movie * Gain
tifffile.imwrite("GainApplied%sCorrection.tif" % Mode,ImageOut.astype(np.float32))