import pandas as pd
import numpy as np
import GeneralLibrary as GL
import matplotlib.pyplot as plt
import time
from tkinter import filedialog

#CSVFolder = filedialog.askdirectory() + "/"
#DataPaths = GL.ScanForFilesInPathByTag(CSVFolder,".csv")

#CSVPaths = filedialog.askopenfilenames(title = "Select Files", filetypes = [("CSV Files", "*.csv")])
#print(CSVPaths)

#SubDirectories = ["0umStepSize/","0_001umStepSize/","0_01umStepSize/","0_1umStepSize/","1umStepSize/","10umStepSize/"]

SubDirectories = ["0umStepSize/","0_001umStepSize/","0_01umStepSize/","0_1umStepSize/","1umStepSize/","10umStepSize/"]# "0umStepSize/" ,
SubSubDirectory = ["MovieHome/","MovieNegX1/","MovieNegX2/","MovieNegXY/","MovieNegY1/","MovieNegY2/","MoviePosX1/","MoviePosX2/", "MoviePosXY/","MoviePosY1/", "MoviePosY2/"]

CSVPaths = []
for Dirs in SubDirectories:
    for Dir2s in SubSubDirectory:
        DataPaths = GL.ScanForFilesInPathByTag("D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/" + Dirs + Dir2s,".csv")
        for Paths in DataPaths:
            CSVPaths.append("D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/" + Dirs + Dir2s + Paths)

def oneD_FFT(XCoords, YCoords):
    NumberSamplesX = len(XCoords)
    NumberSamplesY = len(YCoords)

    kX = np.arange(NumberSamplesX)
    kY = np.arange(NumberSamplesY)

    Tx = NumberSamplesX / Fs
    Ty = NumberSamplesY / Fs

    FreqX = kX/Tx
    FreqX = FreqX[:len(FreqX)//2]
    
    FreqY = kY/Ty
    FreqY = FreqY[:len(FreqY)//2]

    FFTX = np.fft.fft(XCoords)/NumberSamplesX*2
    FFTY = np.fft.fft(YCoords)/NumberSamplesY*2

    FFTX = FFTX[:NumberSamplesX//2]
    FFTY = FFTY[:NumberSamplesY//2]
    
    XInt = np.sum(FFTX) * 1 / Fs
    YInt = np.sum(FFTY) * 1 / Fs
    
    return FFTX, FFTY, FreqX, FreqY

FFTXArray = []
FFTYArray = []

ErrorXArray = []
ErrorYArray = []

XDeviationArray = []
YDeviationArray = []
Iterator = 0

for Dirs in SubDirectories:
    for Dir2s in SubSubDirectory:
        df = pd.DataFrame()
        df2 = pd.DataFrame()
    
        CSVPaths = GL.ScanForFilesInPathByTag("D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/" + Dirs + Dir2s,".csv")
        
        for DataFiles in CSVPaths:
            CSVData = pd.read_csv("D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/" + Dirs + Dir2s + DataFiles)
            #CSVData = pd.read_csv(CSVFolder + DataFiles)
        
            Fs = 600 # Frequency in Hz
        
            DataX = []
            MeanX = CSVData['Peak X Fit (Pixels)'].mean()
            for DataPoints in CSVData['Peak X Fit (Pixels)']:
                DataOut = DataPoints - MeanX
                DataX.append(DataOut)
        
            DataY = []
            MeanY = CSVData['Peak Y Fit (Pixels)'].mean()
            for DataPoints in CSVData['Peak Y Fit (Pixels)']:
                DataOut = DataPoints - MeanY
                DataY.append(DataOut)
        
            Frames = CSVData['Frame']
        
        
            XCoordinates = np.array(DataX) * 65
            YCoordinates = np.array(DataY) * 65
        
        
            FFTX, FFTY, FreqX, FreqY = oneD_FFT(XCoordinates, YCoordinates)
            
            FFTXArray.append(abs(FFTX))
            FFTYArray.append(abs(FFTY))
        
            XDeviationArray.append(XCoordinates)
            YDeviationArray.append(YCoordinates)
            
            Iterator += 1
    
        AveFFTX = abs(np.average(FFTXArray, axis = 0))
        AveFFTY = abs(np.average(FFTYArray, axis = 0))
    
        AveXDeviation = np.average(XDeviationArray, axis = 0)
        AveYDeviation = np.average(YDeviationArray, axis = 0)


        df['Frequency X'] = FreqX
        df['Frequency Y'] = FreqY
        df["X_FFT"] = AveFFTX
        df["Y_FFT"] = AveFFTY

        df2['Frames'] = CSVData['Frame']
        df2["X_DEVIATION"] = AveXDeviation
        df2["Y_DEVIATION"] = AveYDeviation

        df.to_csv("D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/"+"FFTs/"+Dirs.replace("StepSize/","")+Dir2s.replace("/","")+".csv")
        df2.to_csv("D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/"+"Deviation/"+Dirs.replace("StepSize/","")+Dir2s.replace("/","")+".csv")
#AveFFTX = FFTXArray
#AveFFTY = FFTYArray





import Graphing

#XData = [[],[]]
#for values in AveFFTX:
#    XData[0].append(FreqX)
#    XData[1].append(FreqY)

#YData = [[],[]]
#for Data in AveFFTX:
#    YData[0].append(Data)

#for Data2 in AveFFTY:
#    YData[1].append(Data2)


#Labels = [[],[]]
#for Dirs in SubDirectories:
#    Labels[0].append("X-Axis Response (%s)" % Dirs.replace("umStepSize/",""))
#    Labels[1].append("Y-Axis Response (%s)" % Dirs.replace("umStepSize/",""))
#AxisLabels = ("Frequency (Hz)", "Intensity (nm)")
#Title = "Static Stability"

#Graphing.StackedDoubleXY(XData, YData, Labels, AxisLabels, Title, ShareY = False)

XData = (FreqX, FreqY)
YData = (abs(AveFFTX), abs(AveFFTY))
Labels = ("X-Axis Response", "Y-Axis Response")
AxisLabels = ("Frequency (Hz)", "Intensity (nm)")
Title = "Static Stability"

#Graphing.DoubleXY(XData, YData, Labels, AxisLabels, Title, ShareY = False)


XData = (Frames / Fs, Frames / Fs)
YData = (AveXDeviation, AveYDeviation)
Labels = ("Peak X", "Peak Y")
AxisLabels = ("Time (Sec)", "Deviation from Mean (nm)")
Title = "Static Stability"

#Graphing.DoubleXY(XData, YData, Labels, AxisLabels, Title, AxisLimits = (((0,2),(0,2)),((-50, 50), (-100,100))), ShareY = False, ShareX = False)