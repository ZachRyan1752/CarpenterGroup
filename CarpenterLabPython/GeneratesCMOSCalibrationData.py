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

# defs
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


def GainCalibrationFromMovies(FolderPath):
    FolderPath = FolderPath.replace("\n  ","")
    FolderPath = FolderPath.replace("\n ","")
    ## Generate a image containing the gain at each individual pixel
    ## Procedure from:
    ## https://mirametrics.com/tech_note_ccdgain.php
    DataPath = FolderPath + "Calibration/GainCalibration/"
    ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"
    MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif")
    MovieDataArray = []
    AverageMovieValueArray = []

    for FilePath in MoviePaths:
        MovieData = tifffile.imread(DataPath + FilePath)
        MovieDataArray.append(MovieData)
        CompleteMovieAverage = np.average(MovieData)
        AverageMovieValueArray.append(CompleteMovieAverage)

    MovieIndex1 = 0
    NewMovieDataArray = []
    for Movie in MovieDataArray:
        MovieIndex2 = 0
        MovieStacked = np.array(np.shape(Movie))
        MovieStacked += Movie
        for Movie2 in MovieDataArray:
            if MovieIndex1 != MovieIndex2: ## If they are different movies
                if AverageMovieValueArray[MovieIndex1] * 0.9 < AverageMovieValueArray[1] and AverageMovieValueArray[1] < AverageMovieValueArray[MovieIndex1] * 1.1:
                    MovieStacked = np.vstack((MovieStacked,Movie2)) ## Append the two movies
                    MovieDataArray.pop(MovieIndex2) ## Remove the compared movie
            MovieIndex2 += 1
        NewMovieDataArray = []
        MovieIndex1 += 1

    del MovieDataArray ## Removing duplicate data

    MovieDataArray = []
    MovieDataArray = NewMovieDataArray

    del NewMovieDataArray ## Removing duplicate data

    AverageMovieValueArray = []
    for Movies in MovieDataArray:
        CompleteMovieAverage = np.average(Movies)
        AverageMovieValueArray.append(CompleteMovieAverage)

    DarkFrameIndex = np.where(np.min(AverageMovieValueArray) == AverageMovieValueArray) ## np.where retrurns an array filled with the positions where the condition is met
    DarkFrameMovie = MovieDataArray[DarkFrameIndex[0][0]]
    DarkFrameAverage = np.average(DarkFrameMovie, axis = 0)
    MovieDataArray.pop(DarkFrameIndex[0][0])
    AverageMovieValueArray.pop(DarkFrameIndex[0][0])

    Looping = True
    while Looping:
        ## Subtract the background (dark image)
        MoviesBackgroundSubtracted = []
        for Movies in MovieDataArray:
            MovieBackgroundSubtracted = []
            for Frames in Movies:
                FrameBackgroundSubtracted = Frames - DarkFrameAverage
                MovieBackgroundSubtracted.append(FrameBackgroundSubtracted)
            MoviesBackgroundSubtracted.append(MovieBackgroundSubtracted)

        ## Get the ratio of frames and adjust them to be the same mean
        ## Subtract Frame B from Frame A
        ## Get the variance of the noise
        AverageFrameOverIlluminationLevels = []
        VarianceFrameOverIlluminationLevels = []
        for Movies in MoviesBackgroundSubtracted:
            FrameA = Movies[0]
            AdjustedFrameArray = []
            for Frames in Movies:
                RatioAB = np.average(FrameA) / np.average(Frames)
                FrameBAdjusted = Frames * RatioAB
                AdjustedFrameArray.append(FrameBAdjusted)
                
            AveragePixelFrame = np.average(AdjustedFrameArray, axis = 0)
            VariancePixelFrame = np.var(AdjustedFrameArray, axis = 0)
            
            AverageFrameOverIlluminationLevels.append(AveragePixelFrame)
            VarianceFrameOverIlluminationLevels.append(VariancePixelFrame)
        
        
        ## Collecting some variables before clearing memory
        LengthOfMovieArray = len(MovieDataArray)-1

        
        ## Explicitely clearing memory
        del MovieDataArray
        del MoviesBackgroundSubtracted


        gc.collect()

        AverageFrameOverIlluminationLevels = np.array(AverageFrameOverIlluminationLevels) ## Convert list to np arrays for easier slicing
        VarianceFrameOverIlluminationLevels = np.array(VarianceFrameOverIlluminationLevels) ## Convert list to np arrays for easier slicing
        
        Iterator1 = 0
        ColumnSlope = []
        ColumnIntercept = []
        ColumnR2 = []
        
        BigDataArray = []

        AverageDataToProcessPool = []
        VarianceDataToProcessPool = []
        ### index changed here
        while Iterator1 < np.shape(AverageFrameOverIlluminationLevels)[1]: ## Len - 1 converts to index or not
            AverageData = AverageFrameOverIlluminationLevels[0:LengthOfMovieArray,Iterator1] ## Preslicing the arrays, for (theoretically) faster code execution
            VarianceData = VarianceFrameOverIlluminationLevels[0:LengthOfMovieArray,Iterator1]

            AverageDataToProcessPool.append(AverageData)
            VarianceDataToProcessPool.append(VarianceData)
            print("Working on Column Number: %s" % Iterator1)
            Iterator1 += 1
        
        Iterator1 = 0   
        for data in AverageDataToProcessPool:
            Column = Iterator1
            AverageData = AverageDataToProcessPool[Iterator1]
            VarianceData = VarianceDataToProcessPool[Iterator1]

            Column, RowSlope, RowIntercept, RowR2 = FitColumnLinear(Iterator1, LengthOfMovieArray, AverageData, VarianceData)
            BigDataArray.append([Column, RowSlope, RowIntercept, RowR2])

            Iterator1 += 1

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

        VarianceFrameOverIlluminationLevels = VarianceFrameOverIlluminationLevels.astype(np.float32)
        AverageFrameOverIlluminationLevels = AverageFrameOverIlluminationLevels.astype(np.float32)
        SlopeImage = SlopeImage.astype(np.float32)
        InterceptImage = InterceptImage.astype(np.float32)
        R2Image = R2Image.astype(np.float32)
        #print(np.shape(VarianceFrameOverIlluminationLevels),np.shape(AverageFrameOverIlluminationLevels))
        #print("Spacer")
        #print(np.shape(SlopeImage),np.shape(InterceptImage),np.shape(R2Image))

        tifffile.imwrite(ExtractedDataPath+"GainPixelVariance.tif",VarianceFrameOverIlluminationLevels)
        tifffile.imwrite(ExtractedDataPath+"GainPixelAverage.tif",AverageFrameOverIlluminationLevels)
        tifffile.imwrite(ExtractedDataPath+"GainPixelSlope.tif",SlopeImage)
        tifffile.imwrite(ExtractedDataPath+"GainPixelIntercept.tif",InterceptImage)
        tifffile.imwrite(ExtractedDataPath+"GainPixelR2.tif",R2Image)
        
        print("Average gain: %s" % np.average(SlopeImage))
        print("Average fit R2: %s" % np.average(R2Image))

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
    RatioX = SizeXGain / SizeX
    RatioY = SizeYGain / SizeY

    CropXMin = int(SizeXGain * (RatioX - 1) / (2 * RatioX))
    CropXMax = int(SizeXGain * (RatioX + 1) / (2 * RatioX))
    CropYMin = int(SizeYGain * (RatioY - 1) / (2 * RatioY))
    CropYMax = int(SizeYGain * (RatioY + 1) / (2 * RatioY))

    

    #print(CropXMin,CropXMax,CropYMin,CropYMax)
    GainCalibrationImageCropped = GainCalibrationImage[CropXMin:CropXMax,CropYMin:CropYMax]

    BeamCalibrationFile = ExtractedDataPath+"BeamFitParameters.xml"
    with open(BeamCalibrationFile) as BFPFile:
        BFP = BFPFile.read()
        BeamFitXML = BeautifulSoup(BFP, 'xml')

    ## Extracting meaningful fit parameters
    w0XMLData = BeamFitXML.find_all('w0Pixels')
    PeakXData = BeamFitXML.find_all('PeakX')
    PeakYData = BeamFitXML.find_all('PeakY')
    for Data in w0XMLData:
        w0 = np.float32(Data.text)
    for Data in PeakXData:
        PeakX = np.float32(Data.text)
    for Data in PeakYData:
        PeakY = np.float32(Data.text)

    ## Adjusting peaks for pre-cropping of the movies
    PeakX = PeakX - CropXMin
    PeakY = PeakY - CropYMin

    CroppedShape = np.shape(GainCalibrationImageCropped)

    ## I want the cropping to be around the center of the beams gaussian fit center
    BeamCropXMin = int(PeakX - w0)
    BeamCropXMax = int(PeakX + w0)
    BeamCropYMin = int(PeakY - w0)
    BeamCropYMax = int(PeakY + w0)

    print(PeakX,PeakY,BeamCropXMin,BeamCropXMax,BeamCropYMin,BeamCropYMax)
    ## We will use w0 from beam fitting to determine the radius over which the movie will be cropped to
    MovieCountIterator = 0
    for MoviePath in MoviePaths:
        Movie = tifffile.imread(RawDataPath + MoviePath)
        #MovieCalibrated = Movie / GainCalibrationImageCropped ## Applying Gain Calibration

        MovieCalibrated = Movie * 0.24

        MovieCalibratedCropped = MovieCalibrated[0:FrameCount,BeamCropXMin:BeamCropXMax,BeamCropYMin:BeamCropYMax]

        tifffile.imwrite(CalibratedMoviesPath+"Movie%s.tif" % MovieCountIterator,MovieCalibratedCropped.astype(np.float32))
        MovieCountIterator += 1


    ## Applying gain calibration to the blank file
    BlankFrame = tifffile.imread(FolderPath + "Calibration/CalibrationData/BlankAverage.tif")
    BlankFrame = BlankFrame * 0.24

    tifffile.imwrite(FolderPath + "Calibration/CalibrationData/BlankAverageGainCorrected.tif",BlankFrame.astype(np.float32))