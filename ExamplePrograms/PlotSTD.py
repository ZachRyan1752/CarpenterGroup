import GeneralLibrary as GL
import time
import pandas as pd
import Graphing as Gph
import numpy as np





KnownFile = "D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/ZDrift/100nmStepSizeZStack_MMStack_Pos0.ome/"
ExperimentalFile = "D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/ZDrift/StitchedMovie/"

Folder = KnownFile
CSVPaths = GL.ScanForFilesInPathByTag(Folder,".csv")
XDataArray = []
YDataArray1 = []
YDataArray2 = []
LabelArray1 = []
LabelArray2 = []

for files in CSVPaths:
    df = pd.read_csv(Folder+files)

    XDataArray.append(df['Frame']*100)
    YDataArray1.append(df['Sigma X (nm)'])
    YDataArray2.append(df['Sigma Y (nm)'])
    LabelArray1.append("")
    LabelArray1.append("")

YDataArray1 = [np.average(YDataArray1, axis = 0)]
YDataArray2 = [np.average(YDataArray2, axis = 0)]
LabelArray1 = [""]
LabelArray2 = [""]
XDataArray = [df['Frame']*100]

Gph.Init()
Gph.ACSFormatSingleColumn()
Gph.XYPlot(XDataArray, YDataArray1, LabelArray1, AxisLabels = ["Z-Height (nm)", "Standard Deviation (nm)"])
Gph.XYPlot(XDataArray, YDataArray2, LabelArray2, AxisLabels = ["Z-Height (nm)", "Standard Deviation (nm)"])