## This library contains code that automatically calculates the gain coefficient, beam profile (For later normalizing data), and other camera calibrations for sCMOS cameras
# Import libraries
import tifffile
import os
import GeneralLibrary as GL
import fits
import numpy as np
from scipy.stats import linregress
import scipy.optimize as opt
import concurrent
import gc
import math
from scipy.fft import fft2, fftshift
import cv2
import ImageProcessing as IP
import lxml.etree as etree
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup
import concurrent
import time
import matplotlib.pyplot as plt

# defs
def LoadMovieAndStitchMovie(MoviePath, MovieStitched):
    Movie = tifffile.imread(MoviePath)
    for Frames in Movie:
        MovieStitched.append(Frames)
        
    return MovieStitched

def LinearFitByColumn(x,XData,YData):
    Slope, Intercept, R2, P, STD_ERR = linregress(XData,YData)
    return x, Slope, Intercept, R2


def FitColumnLinear(Column, LengthOfMovieArray, AverageData, VarianceData):
    ## This function splits a column into pixels for linear fitting, it will be more useful when I implement asynchronous processing
    Iterator2 = 0
    RowSlope = []
    RowIntercept = []
    RowR2 = []
    ### Index changed here as well
    while Iterator2 < np.shape(AverageData)[1]:
        PixelVarianceStack = VarianceData[0:LengthOfMovieArray,Iterator2]
        PixelAverageStack = AverageData[0:LengthOfMovieArray,Iterator2]  ## 0:LengthOfMovieArray,Column,
                
        Column, Slope, Intercept, R2 = LinearFitByColumn(Column, PixelAverageStack, PixelVarianceStack)
                
        RowSlope.append(Slope)
        RowR2.append(R2)
        RowIntercept.append(Intercept)
        Iterator2 += 1

    RowSlope = np.array(RowSlope)
    RowIntercept = np.array(RowIntercept)
    RowR2 = np.array(RowR2)
    
    print("Finished Column: %s" % Column)

    return Column, RowSlope, RowIntercept, RowR2

def FitColumnLinearWrapper(args):
    Column, RowSlope, RowIntercept, RowR2 = FitColumnLinear(args[0], args[1], args[2], args[3])

    return Column, RowSlope, RowIntercept, RowR2

def GetMovieAveragesFromFilePaths(FolderPath,Paths):
    LengthOfPaths = len(Paths)
    AverageMovieArray = []
    MovieIndex = 0
    for Path in Paths:
        MovieData = tifffile.imread(FolderPath+"Calibration/GainCalibration/"+Path)
        MovieAverage = np.average(MovieData)
        MovieShape = np.shape(MovieData)
        AverageMovieArray.append(MovieAverage)
        MovieIndex += 1
        print("Averaging Movie %s out of %s" % (MovieIndex, LengthOfPaths))
    return AverageMovieArray, MovieShape

def StitchMoviesAndReturnAverageAndVarianceArrays(FolderPath,Paths,PathIndexes,AverageValueArray,MovieShape):
    StitchedMovie = []
    IndexIndex = 0
    
    for Indicies in PathIndexes:
        StitchedMovie = LoadMovieAndStitchMovie(FolderPath+"Calibration/GainCalibration/"+Paths[Indicies],StitchedMovie)
        IndexIndex += 1
        print("Stitching Movie %s out of %s" % (IndexIndex, len(PathIndexes)))
        
    StitchedMovie = np.array(StitchedMovie)

    
    StitchedMovieAverage = np.average(StitchedMovie, axis = 0)
    StitchedMovieVariance = np.var(StitchedMovie, axis = 0)
    
    return StitchedMovieAverage, StitchedMovieVariance

def FrameABSubtraction(Frames, FrameAAverage, FrameCounter, MovieShape):
    RatioAB = FrameAAverage / np.average(Frames)
    FrameBAdjusted = Frames * RatioAB
    #FrameBAdjusted = np.expand_dims(FrameBAdjusted, axis = 0)
    print("Starting Frame %s out of %s" % (FrameCounter + 1, MovieShape[0])) ## Technically it is already done at this point

    return FrameBAdjusted

def FrameABSubtractionWrapper(args):
    ## https://stackoverflow.com/questions/6785226/pass-multiple-parameters-to-concurrent-futures-executor-map
    FrameBAdjusted = FrameABSubtraction(args[0],args[1],args[2],args[3])

    return FrameBAdjusted

def StitchMoviesAndReturnAverageAndVarianceArraysAfterBackgroundSubtractionAndABRatio(FolderPath,Paths,PathIndexes,AverageValueArray,MovieShape,BackgroundAverage):
    StitchedMovie = []
    IndexIndex = 0    

    for Indicies in PathIndexes:
        StitchedMovie = LoadMovieAndStitchMovie(FolderPath+"Calibration/GainCalibration/"+Paths[Indicies],StitchedMovie)
        
        IndexIndex += 1
        print("Stitching Movie %s out of %s" % (IndexIndex, len(PathIndexes)))
        
    StitchedMovie = np.array(StitchedMovie)
    StitchedMovie = StitchedMovie - BackgroundAverage
    MovieShape = np.shape(StitchedMovie)
    print(MovieShape,PathIndexes)
    time.sleep(2)
    ## Get the ratio of frames and adjust them to be the same mean
    ## Subtract Frame B from Frame A
    ## Get the variance of the noise
    #StitchedMovieAB = np.empty((0,MovieShape[1],MovieShape[2]))
    StitchedMovieAB = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=61) as executor: 
        ## This is the hard limit for workers on windows: 
        ## See: https://github.com/python/cpython/issues/71090
        ## See: https://github.com/psf/black/issues/564
        FrameCounter = 0
        DataPool = []
        for Frames in StitchedMovie:
            if FrameCounter == 0:
                FrameA = Frames
                FrameAAverage = np.average(FrameA)
            else:
                ProcessData = [Frames,FrameAAverage, FrameCounter, MovieShape]
                DataPool.append(ProcessData)
                #ThreadData = executor.map(FrameABSubtraction, Frames, FrameAAverage, FrameCounter, MovieShape) # .submit
                #ThreadData = {executor.submit(FrameABSubtraction, Frame2, FrameAAverage, FrameCounter, MovieShape) for Frame2 in StitchedMovie}
                #break
            FrameCounter += 1

        ## https://stackoverflow.com/questions/6785226/pass-multiple-parameters-to-concurrent-futures-executor-map
        ThreadReturnedData = executor.map(FrameABSubtractionWrapper, DataPool)
        FinishedProcessCounter = 0
        for FrameBAdjusted in ThreadReturnedData:  
            StitchedMovieAB.append(FrameBAdjusted)
            FinishedProcessCounter += 1
            print("Finished AB Frame %s out of %s" % (FinishedProcessCounter, np.shape(StitchedMovie)[0]))

    StitchedMovieAB = np.array(StitchedMovieAB)
    print("Stitched Array Shape: %s %s %s" % np.shape(StitchedMovieAB))

    StitchedMovieAverage = np.average(StitchedMovieAB, axis = 0)
    StitchedMovieVariance = np.var(StitchedMovieAB, axis = 0)
    
    return StitchedMovieAverage, StitchedMovieVariance


def GainCalibrationFromMovies(FolderPath,IlluminationLevels):
    StartTime = time.time()
    FolderPath = FolderPath.replace("\n  ","")
    FolderPath = FolderPath.replace("\n ","")
    ## Generate a image containing the gain at each individual pixel
    ## Procedure from:
    ## https://mirametrics.com/tech_note_ccdgain.php
    ## Essentially the same as 2-frame calibration, but on a per-pixel basis, this removes one of the two offsets for the fitting
    DataPath = FolderPath + "Calibration/GainCalibration/"
    ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"
    MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif")


    LengthMovieList = len(MoviePaths)
    
    AverageMovieValueArray, MovieShape = GetMovieAveragesFromFilePaths(FolderPath, MoviePaths)

    DarkImageIndex = 0
    MinimumAverageValue = np.min(AverageMovieValueArray)
    IndexArray = []
    for Average in AverageMovieValueArray:
        if Average > MinimumAverageValue * 0.95 and MinimumAverageValue * 1.05 > Average :
            IndexArray.append(DarkImageIndex)
        DarkImageIndex += 1

    print("Getting Dark Movie")
    DarkFrameAverage, DarkFrameVariance = StitchMoviesAndReturnAverageAndVarianceArrays(FolderPath,MoviePaths,IndexArray,AverageMovieValueArray,MovieShape)
    
    AverageFrameOverIlluminationLevels = []
    VarianceFrameOverIlluminationLevels = []

    AverageFrameOverIlluminationLevels.append(np.array(DarkFrameAverage))
    VarianceFrameOverIlluminationLevels.append(np.array(DarkFrameVariance))

    RemovedCount = 0
    for Indexes in IndexArray[::-1]:
        print(AverageMovieValueArray[Indexes],MoviePaths[Indexes],MinimumAverageValue) 
        AverageMovieValueArray.pop(Indexes)
        MoviePaths.pop(Indexes)
        RemovedCount += 1

    print("Getting Movie Averages and Variances")
    print(AverageMovieValueArray)

    

    while len(AverageMovieValueArray) > 0:
        MinimumAverageValue = np.min(AverageMovieValueArray)
        IndexArray = []
        ImageIndex = 0
        for Average in AverageMovieValueArray:
            if Average > MinimumAverageValue * 0.95 and MinimumAverageValue * 1.05 > Average :
                print(Average, MinimumAverageValue, ImageIndex)
                IndexArray.append(ImageIndex)
            ImageIndex += 1

        FrameAverage, FrameVariance = StitchMoviesAndReturnAverageAndVarianceArraysAfterBackgroundSubtractionAndABRatio(FolderPath,MoviePaths,IndexArray,AverageMovieValueArray,MovieShape,DarkFrameAverage)

        RemovedCount = 0
        for Indexes in IndexArray[::-1]:
            print(AverageMovieValueArray[Indexes],MoviePaths[Indexes],MinimumAverageValue) 
            AverageMovieValueArray.pop(Indexes)
            MoviePaths.pop(Indexes)
            RemovedCount += 1
        
        AverageFrameOverIlluminationLevels.append(np.array(FrameAverage-FrameAverage))
        VarianceFrameOverIlluminationLevels.append(np.array(FrameVariance))



    AverageFrameOverIlluminationLevels2 = np.empty(np.shape(MovieShape))
    for Frames in AverageFrameOverIlluminationLevels:
        AverageFrameOverIlluminationLevels2 = np.vstack((AverageFrameOverIlluminationLevels2, Frames))

    VarianceFrameOverIlluminationLevels2 = np.empty(np.shape(MovieShape))
    for Frames in VarianceFrameOverIlluminationLevels:
        VarianceFrameOverIlluminationLevels2 = np.vstack((VarianceFrameOverIlluminationLevels2, Frames))

    AverageFrameOverIlluminationLevels = AverageFrameOverIlluminationLevels2
    VarianceFrameOverIlluminationLevels = VarianceFrameOverIlluminationLevels2
    
    
    #AverageFrameOverIlluminationLevels = np.array(AverageFrameOverIlluminationLevels) ## Convert list to np arrays for easier slicing

    #AverageFrameOverIlluminationLevels = np.array(AverageFrameOverIlluminationLevels) ## Convert list to np arrays for easier slicing
    #VarianceFrameOverIlluminationLevels = np.array(VarianceFrameOverIlluminationLevels) ## Convert list to np arrays for easier slicing

    print("Average and Variance Array Shapes")
    print(np.shape(AverageFrameOverIlluminationLevels),np.shape(VarianceFrameOverIlluminationLevels))
    #time.sleep(5)
    Iterator1 = 0
    ColumnSlope = []
    ColumnIntercept = []
    ColumnR2 = []
        
    BigDataArray = []

    AverageDataToProcessPool = []
    VarianceDataToProcessPool = []
    
    ### index changed here
    while Iterator1 < np.shape(AverageFrameOverIlluminationLevels)[1]: ## Len - 1 converts to index or not
        AverageData = AverageFrameOverIlluminationLevels[0:IlluminationLevels+1,Iterator1] ## Preslicing the arrays, for (theoretically) faster code execution
        VarianceData = VarianceFrameOverIlluminationLevels[0:IlluminationLevels+1,Iterator1]

        AverageDataToProcessPool.append(AverageData)
        VarianceDataToProcessPool.append(VarianceData)
        print("Working on Column Number: %s" % Iterator1)
        Iterator1 += 1
        
    with concurrent.futures.ProcessPoolExecutor(max_workers=61) as executor: 
        ## This is the hard limit for workers on windows: 
        ## See: https://github.com/python/cpython/issues/71090
        ## See: https://github.com/psf/black/issues/564
        ArrayCounter = 0
        DataPool = []
        for Data in AverageDataToProcessPool:
            Column = ArrayCounter
            AverageData = AverageDataToProcessPool[ArrayCounter]
            VarianceData = VarianceDataToProcessPool[ArrayCounter]

            ProcessData = [ArrayCounter, IlluminationLevels, VarianceData, AverageData] ## Variance and average have been swapped
            DataPool.append(ProcessData)
            ArrayCounter += 1

        ## https://stackoverflow.com/questions/6785226/pass-multiple-parameters-to-concurrent-futures-executor-map
        ThreadReturnedData = executor.map(FitColumnLinearWrapper, DataPool)
        FinishedProcessCounter = 0
        for Data in ThreadReturnedData:  
            Column = Data[0] 
            RowSlope = Data[1]
            RowIntercept = Data[2]
            RowR2 = Data[3]
            BigDataArray.append([Column, RowSlope, RowIntercept, RowR2])
            FinishedProcessCounter += 1
            print("Finished Column %s out of %s" % (FinishedProcessCounter, MovieShape[1]))

    ## Sorting the data
    ## This will be used when I implement asynchronous code execution

    Iterator2 = 0
    for elements in BigDataArray:
        ColumnSlope.append(0)
        ColumnIntercept.append(0)
        ColumnR2.append(0)
        print(Iterator2)
        Iterator2 += 1

    Iterator3 = 0
    for elements in BigDataArray:
        Index = elements[0] 
        ColumnSlope[Index] = elements[1]
        ColumnIntercept[Index] = elements[2]
        ColumnR2[Index] = elements[3]
        print(Iterator3)
        Iterator3 += 1

    SlopeImage = np.array(ColumnSlope)
    InterceptImage = np.array(ColumnIntercept)
    R2Image = np.array(ColumnR2)

    VarianceFrame = np.array(VarianceFrameOverIlluminationLevels[0])

    VxG = VarianceFrame * SlopeImage
    VarianceFrameOverIlluminationLevels = VarianceFrameOverIlluminationLevels.astype(np.float32)
    AverageFrameOverIlluminationLevels = AverageFrameOverIlluminationLevels.astype(np.float32)
    print(np.shape(VarianceFrameOverIlluminationLevels),np.shape(AverageFrameOverIlluminationLevels))

    SlopeImage = SlopeImage.astype(np.float32)
    InterceptImage = InterceptImage.astype(np.float32)
    R2Image = R2Image.astype(np.float32)

    tifffile.imwrite(ExtractedDataPath+"Variance.tif",VarianceFrame.astype(np.float32))
    tifffile.imwrite(ExtractedDataPath+"ElectronReadoutNoise.tif",VxG.astype(np.float32))

    tifffile.imwrite(ExtractedDataPath+"GainPixelVariance.tif",VarianceFrameOverIlluminationLevels)
    tifffile.imwrite(ExtractedDataPath+"GainPixelAverage.tif",AverageFrameOverIlluminationLevels)
    tifffile.imwrite(ExtractedDataPath+"GainPixelSlope.tif",SlopeImage)
    tifffile.imwrite(ExtractedDataPath+"GainPixelIntercept.tif",InterceptImage)
    tifffile.imwrite(ExtractedDataPath+"GainPixelR2.tif",R2Image)
        
    print("Average gain: %s" % np.average(SlopeImage))
    print("Average fit R2: %s" % np.average(R2Image))
    StopTime = time.time()
    print("Time To Execute Gain Calibration: %s Seconds" % (StopTime - StartTime))

def SpatialCalibrationFromRuling(FolderPath, LPM, CameraPixelSizeum):
    DataPath = FolderPath + "Calibration/Resolution/"
    ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"
    MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif")

    for MoviePath in MoviePaths:
        Movie = tifffile.imread(DataPath + MoviePath)

    ## https://www.geeksforgeeks.org/how-to-find-the-fourier-transform-of-an-image-using-opencv-python/
    FFT = cv2.dft(np.float32(Movie), flags = cv2.DFT_COMPLEX_OUTPUT)
    FFTCentered = np.fft.fftshift(FFT)

    FFTMagnitude = 20*np.log(cv2.magnitude(FFTCentered[:,:,0],FFTCentered[:,:,1]))

    FFTNormalized = cv2.normalize(FFTMagnitude, None, 0, 255, cv2.NORM_MINMAX, cv2.CV_8UC1)
    
    tifffile.imwrite(ExtractedDataPath+"RonchiRulingFFT.tif",FFTNormalized.astype(np.float32))
    
    shapeFFT = np.shape(FFTNormalized)[0]
    shapeFFTCrop = int(shapeFFT / 2 - 20)
    EditedFFTPeaksPicked, ArrayOfPeaks = IP.FindPeaks(FFTNormalized,8,200,shapeFFTCrop)

    tifffile.imwrite(ExtractedDataPath+"FFTPickedPeaks.tif",EditedFFTPeaksPicked.astype(np.float32))

    CenterPeak = ArrayOfPeaks[0] ## Center should always be the largest peak (no spatial frequency)
    Peak1 = ArrayOfPeaks[1]
    Peak2 = ArrayOfPeaks[2]

    CenterCoords = [CenterPeak[1],CenterPeak[2]]
    Peak1Coords = [Peak1[1],Peak1[2]]
    Peak2Coords = [Peak2[1],Peak2[2]]
    Distance1 = GL.Distance(CenterCoords,Peak1Coords)
    Distance2 = GL.Distance(CenterCoords,Peak2Coords)

    AveDistance = (Distance1 + Distance2) / 2 / np.shape(Movie)[0]

    nmPerLine = 1 / LPM * 1000 * 1000## nmPerLine
    umPerLine = 1 / LPM * 1000

    PixelSize = [AveDistance * nmPerLine, CameraPixelSizeum / (AveDistance * nmPerLine) * 1000]

    Header = "PixelCalibration"
    SubHeaders = ["PixelSize","Magnification"]
    SubHeadersDataName = ["PixelSize","Magnification"]
    GL.MakeXMLWithHeaderbs4(Header, SubHeaders, PixelSize, ExtractedDataPath, "PixelSize.xml")


def GetBeamProfile(FolderPath, BeamIntensityAtSample):
    DataPath = FolderPath + "Calibration/BeamProfile/"
    ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"
    MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif")

    for MoviePath in MoviePaths:
        Movie = tifffile.imread(DataPath + MoviePath)

    MovieAveraged = np.average(Movie, axis = 0)
    tifffile.imwrite(ExtractedDataPath+"AverageBeamProfile.tif",MovieAveraged.astype(np.float32))

    PeakMax = np.max(MovieAveraged)
    PeakCoords = np.where(PeakMax == MovieAveraged)
    PeakCoords = list(zip(PeakCoords))

    
    PixelCalibrationFile = ExtractedDataPath+"PixelSize.xml"
    with open(PixelCalibrationFile) as PCFFile:
        PCF = PCFFile.read()
        PixelSizeXML = BeautifulSoup(PCF, 'xml')

    PixelSizeXMLData = PixelSizeXML.find_all('PixelSize')
    for Data in PixelSizeXMLData:
        PixelSize = np.float32(Data.text)

    ## Finding where the image rises above the background on one size to crop the image for faster fitting
    ## Where the image rised above background - 50 pixels
    OriginalMovieShape = np.shape(MovieAveraged)
    FitBoundsXMin = 0 #np.where(MovieAveraged[PeakCoords[0][0][0]] > 150)[0][0] - 250 
    #FitBoundsXMax = OriginalMovieShape[0] - FitBoundsXMin ## Assumes the peak is perfectly in the center
    FitBoundsYMin = 0 #np.where(MovieAveraged[PeakCoords[1][0][0]] > 150)[0][0] - 250
    #FitBoundsYMax = OriginalMovieShape[1] - FitBoundsXMin ## Assumes the peak is perfectly in the center
    #MovieAveraged = MovieAveraged[FitBoundsXMin:FitBoundsXMax,FitBoundsYMin:FitBoundsYMax]
    #print(FitBoundsXMin,FitBoundsXMax,FitBoundsYMin,FitBoundsYMax)


    ImageShape = np.shape(MovieAveraged)
    xaxis = np.linspace(0,ImageShape[0],ImageShape[0])
    yaxis = np.linspace(0,ImageShape[1],ImageShape[1])
    xaxis, yaxis = np.meshgrid(xaxis, yaxis)

    RaveledData = MovieAveraged.ravel()

    InitialGuess = (PeakMax, PeakCoords[0][0][0] - FitBoundsXMin, PeakCoords[1][0][0] - FitBoundsYMin, 60, 60, 0, 100)
    
    RGLower = [0,0,0,50,50,-0.5,50]
    RGUpper = [10000,2304,2034,70,70,0.5,250]

    ReasonableGuess = (RGLower,RGUpper)

    print(InitialGuess,np.shape(MovieAveraged))

    popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis, yaxis), RaveledData, p0 = InitialGuess, maxfev = 25000, bounds = ReasonableGuess )
    FittedData = fits.twoD_Gaussian((xaxis, yaxis), *popt).reshape(ImageShape[0],ImageShape[1])
    R2 = GL.getR2(MovieAveraged, FittedData)
    print(popt)

    FittingParametersArray = []
    FittingParametersArray.extend(popt)
    FittingParametersArray[1] += FitBoundsXMin
    FittingParametersArray[2] += FitBoundsYMin
    FittingParametersArray.append(R2)
    print(FittingParametersArray)

    FittingParametersArray.extend([OriginalMovieShape[0],OriginalMovieShape[1]])
    print(FittingParametersArray)

    ## Left off here
    w0Pixels = (FittingParametersArray[3] + FittingParametersArray[4]) / 2 * 2
    w0 = (FittingParametersArray[3] + FittingParametersArray[4]) / 2 * 2 * PixelSize / 1e7
    LaserPeakIntensity = BeamIntensityAtSample / (math.pi / 2 * ((w0) ** 2)) / 1000 ## Value in W per cm^2

    FittingParametersArray.extend([LaserPeakIntensity, w0, w0Pixels, PixelSize])
    Header = "BeamProfileFitParameters"
    SubHeaders = ["Amplitude","PeakX","PeakY","SigmaXPixels","SigmaYPixels","ThetaRadians","Offset","R2","MovieShapeX","MovieShapeY","LaserPeakIntensityWPercmPercm","w0cm","w0Pixels","PixelSizenm"]
    GL.MakeXMLWithHeaderbs4(Header, SubHeaders, FittingParametersArray, ExtractedDataPath, "BeamFitParameters.xml")

def BlankAverage(FolderPath):
    BlankMoviesPath = FolderPath + "Calibration/Background/"
    BlankPath = GL.ScanForFilesInPathByTag(BlankMoviesPath,".tif")

    for MoviePath in BlankPath:
        BlankMovie = tifffile.imread(BlankMoviesPath+MoviePath)
    
    BlankAverageFrame = np.average(BlankMovie, axis = 0)

    BlankMovieAveragePath = FolderPath + "Calibration/CalibrationData/BlankAverage.tif"
    tifffile.imwrite(BlankMovieAveragePath,BlankAverageFrame.astype(np.float32))

def ApplyCalibrationsToMovies(FolderPath):
    RawDataPath = FolderPath + "RawData/"
    CalibratedMoviesPath = FolderPath + "CorrectedCroppedData/"
    ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"
    MoviePaths = GL.ScanForFilesInPathByTag(RawDataPath,".tif")


    for MoviePath in MoviePaths:
        MovieExample = tifffile.imread(RawDataPath + MoviePath)
        break


    MovieShape = np.shape(MovieExample)
    FrameCount = MovieShape[0]
    SizeX = MovieShape[1]
    SizeY = MovieShape[2]
    GainCalibrationImage = tifffile.imread(ExtractedDataPath+"GainPixelSlope.tif")
    SizeXGain = np.shape(GainCalibrationImage)[0]
    SizeYGain = np.shape(GainCalibrationImage)[1]

    ## Assume the cropping (If any) is done centered
    CenterX = int(1/2*SizeXGain)
    CenterY = int(1/2*SizeYGain)

    CropXMin = int(CenterX - 1/2 * SizeX)
    CropXMax = int(CenterX + 1/2 * SizeX)
    CropYMin = int(CenterY - 1/2 * SizeY)
    CropYMax = int(CenterY + 1/2 * SizeY)

    #CropXMin = int(SizeXGain * (RatioX - 1) / (2 * RatioX))
    #CropXMax = int(SizeXGain * (RatioX + 1) / (2 * RatioX))
    #CropYMin = int(SizeYGain * (RatioY - 1) / (2 * RatioY))
    #CropYMax = int(SizeYGain * (RatioY + 1) / (2 * RatioY))


    GainCalibrationImageCropped = GainCalibrationImage[CropXMin:CropXMax,CropYMin:CropYMax]

    BeamCalibrationFile = ExtractedDataPath+"BeamFitParameters.xml"
    with open(BeamCalibrationFile) as BFPFile:
        BFP = BFPFile.read()
        BeamFitXML = BeautifulSoup(BFP, 'xml')

    ## Extracting meaningful fit parameters
    ## For whatever reason X and Y are swapped?
    ## Figure this out eventually
    w0XMLData = BeamFitXML.find_all('w0Pixels')
    PeakXData = BeamFitXML.find_all('PeakX')
    PeakYData = BeamFitXML.find_all('PeakY')
    for Data in w0XMLData:
        w0 = np.float32(Data.text)
    for Data in PeakXData:
        PeakY = np.float32(Data.text)
    for Data in PeakYData:
        PeakX = np.float32(Data.text)

    ## Adjusting peaks for pre-cropping of the movies
    print(PeakX, PeakY)
    PeakX = PeakX - CropXMin
    PeakY = PeakY - CropYMin

    CroppedShape = np.shape(GainCalibrationImageCropped)
    
    ## We will use w0 from beam fitting to determine the radius over which the movie will be cropped to
    ## I want the cropping to be around the center of the beams gaussian fit center
    BeamCropXMin = int(PeakX - 3 * w0)
    BeamCropXMax = int(PeakX + 3 * w0)
    BeamCropYMin = int(PeakY - 3 * w0)
    BeamCropYMax = int(PeakY + 3 * w0)

    print(PeakX,PeakY,BeamCropXMin,BeamCropXMax,BeamCropYMin,BeamCropYMax)
    
    MovieCountIterator = 0
    for MoviePath in MoviePaths:
        Movie = tifffile.imread(RawDataPath + MoviePath)
        MovieCalibrated = Movie / GainCalibrationImageCropped ## Applying Gain Calibration

        #MovieCalibrated = Movie * 0.24

        MovieCalibratedCropped = MovieCalibrated[0:FrameCount,BeamCropXMin:BeamCropXMax,BeamCropYMin:BeamCropYMax]

        tifffile.imwrite(CalibratedMoviesPath+"Movie%s.tif" % MovieCountIterator,MovieCalibratedCropped.astype(np.float32))
        MovieCountIterator += 1


    ## Applying gain calibration to the blank file
    BlankFrame = tifffile.imread(FolderPath + "Calibration/CalibrationData/BlankAverage.tif")
    BlankFrame = BlankFrame * 0.24 ## Not the cropped gain calibration image because the code is not yet implemented yet

    tifffile.imwrite(FolderPath + "Calibration/CalibrationData/BlankAverageGainCorrected.tif",BlankFrame.astype(np.float32))