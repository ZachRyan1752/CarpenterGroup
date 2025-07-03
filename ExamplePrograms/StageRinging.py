## Goals:
## Sample the movement jitter of the stage at various stepsizes, tested in +x, -x, +y, -x, +xy, -xy directions
## How:
## 1) Pick 4 ROIs at the vertexes of a square, names ROI1, ROI2, ROI3, ROI4
## 2) Move in this pattern 1 --> 2 --> 3 --> 4 --> 3 --> 2 --> 1

## 2) Move to a second ROI, localize beads at every frame --> Examin how peak position changes over time

import MicroscopeController as MC
import ImageProcessing as IP
import GeneralLibrary as GL
import numpy as np
import time
import random
import tifffile

if __name__ == "__main__": ## To prevent recursive program calling
    ## Setting up micromanager control of components
    MicroscopeDeviceAdapterPath = ["C://Program Files//Micro-Manager-2.0//"]
    ExtractedFolderPath = "RawData/"
    ConfigurationPath = "C:/Users/ryan.1752/Desktop/UManagerFiles/Configuration/AbleTable_TrainingScope_HamCam_PIStage_06_05_2025.cfg"
    ImagingSystem1 = MC.ImagingSystem(MicroscopeDeviceAdapterPath, ExtractedFolderPath, ConfigurationPath, 0)

    ## What step sizes (in um) we want to test
    StepSizeArray = [10] ## [0,1e-3,1e-2,1e-1,1,10]
    ## For whatever reason this code breaks with more than 2 step sizes, fix this later

    ## Data collection settings
    FrameRate = 600
    Frames = 1200
    ExposureTime = 1 / FrameRate * 1000
    Duration = Frames / FrameRate

    ## Some formatting to indicate what names to assign to each collection
    NameIndexArray = ["Home","PosX1","PosY1","NegX1","NegY1","PosY2","PosX2","NegY2","NegX2","PosXY","NegXY"]

    ## Creating subfolders if they do not exist
    for StepSize in StepSizeArray: ## This is done separately, PYMMCORE does NOT like freshly created folders
        break
        if StepSize == 0.1: ## Addressing special cases for decimal step sizes
            GL.CreateIfNotExist("RawData/%sumStepSize/" % "0_1")
        elif StepSize == 0.01:
            GL.CreateIfNotExist("RawData/%sumStepSize/" % "0_01")
        elif StepSize == 0.001:
            GL.CreateIfNotExist("RawData/%sumStepSize/" % "0_001")
        else:
            GL.CreateIfNotExist("RawData/888umStepSize%s/" % StepSize)

    ## Generating movie names for each movie at every step size
    for StepSize in StepSizeArray:
        break
        MovieNamesArray = []
        i = 0
        while i < 11:
            if StepSize == 0.1: ## Addressing special cases for decimal step sizes
                MovieNamesArray.append("RawData/%sumStepSize/Movie%s.tif" % ("0_1", NameIndexArray[i]))
            elif StepSize == 0.01:
                MovieNamesArray.append("RawData/%sumStepSize/Movie%s.tif" % ("0_01", NameIndexArray[i]))
            elif StepSize == 0.001:
                MovieNamesArray.append("RawData/%sumStepSize/Movie%s.tif" % ("0_001", NameIndexArray[i]))
            else:
                MovieNamesArray.append("RawData/888umStepSize%s/Movie%s.tif" % (StepSize, NameIndexArray[i]))
            i += 1
        
        print(MovieNamesArray)
        time.sleep(2)
        
        PosArray = ImagingSystem1.GenerateSquareLoopReverseReturnDiagonalSequence(ImagingSystem1.GetStagePosition(),StepSize)
 
        if StepSize == 0: ## If the step size is 0, wait 5 seconds to remove any ringing caused by previous movement
            DelayConstant = 0
        else: 
            DelayConstant = 0


        #for Positions in PosArray:
        #    Movie = ImagingSystem1.CollectMovie(ExposureTime, FrameRate, int(Duration), ROI = [1008, 1008, 288, 288], Position = Positions)
        #    tifffile.imwrite(MovieNamesArray[0],Movie)
        #    MovieNamesArray.pop(0)
        ## Collect the images with the aformentioned settings and file names
        ImagingSystem1.MDAMovieSequence(ExposureTime, int(Duration), FrameRate, MovieNamesArray, ROI = [1008,1008,288,288], Positions = PosArray, Delay = DelayConstant)

    FrameRate = 800
    NumberOfFrames = 2 * FrameRate
    
    ExposureTime = 1 / FrameRate * 1000
    Duration = int(NumberOfFrames/FrameRate)
    if Duration == 0:
        Duration = 1
    MovieNamesArray = ["RawData/Ultrafast.tif"]
    PosArray = [[70,70,70]]
    #RegionOfInterest = [0,0,2304,2304]
    #RegionOfInterest = [0,128,2304,2048]
    #RegionOfInterest = [0,640,2304,1024]
    #RegionOfInterest = [0,896,2304,512]
    RegionOfInterest = [0,1024,2304,256]

    ImagingSystem1.MDAMovieSequence(ExposureTime, int(Duration), FrameRate, MovieNamesArray, ROI = RegionOfInterest, Positions = PosArray)
    MovieNamesArray = ["RawData/Ultrafast.tif"]
    PosArray = [[70,70,70]]
    ImagingSystem1.MDAMovieSequence(ExposureTime, int(Duration), FrameRate, MovieNamesArray, ROI = RegionOfInterest, Positions = PosArray)
    MovieNamesArray = ["RawData/Ultrafast.tif"]
    PosArray = [[70,70,70]]
    ImagingSystem1.MDAMovieSequence(ExposureTime, int(Duration), FrameRate, MovieNamesArray, ROI = RegionOfInterest, Positions = PosArray)
    MovieNamesArray = ["RawData/Ultrafast.tif"]
    PosArray = [[70,70,70]]
    ImagingSystem1.MDAMovieSequence(ExposureTime, int(Duration), FrameRate, MovieNamesArray, ROI = RegionOfInterest, Positions = PosArray)
    MovieNamesArray = ["RawData/Ultrafast.tif"]
    PosArray = [[70,70,70]]
    ImagingSystem1.MDAMovieSequence(ExposureTime, int(Duration), FrameRate, MovieNamesArray, ROI = RegionOfInterest, Positions = PosArray)


