import GeneralLibrary as GL
import Graphing as Gph
import tifffile

DataPaths = ["UltraquietGainCalibrationNoCorrection/ElectronReadoutNoise.tif","StandardGainCalibrationNoCorrection/ElectronReadoutNoise.tif","FastGainCalibrationNoCorrection/ElectronReadoutNoise.tif"]
#DataPaths.append(GL.AskForFiles()[0])
#DataPaths.append(GL.AskForFiles()[0])
#DataPaths.append(GL.AskForFiles()[0])

DataArray = []
LabelArray = ["Ultraquiet","Standard","Fast"]
for Files in DataPaths:
    Data = tifffile.imread(Files).ravel()
    DataArray.append(Data)
    #LabelArray.append(Files)

Bins = 1000
HistogramFormat = "step"
Gph.Init()
Gph.ACSFormatSingleColumn()
Gph.Histogram(DataArray, LabelArray, Bins, HistogramFormat, AxisLabels = ["Readout Noise (Electrons)","Count (Pixels)"], AxisScales = ["linear","log"], AxisLimits = [(0,5),(0,1e5)])
