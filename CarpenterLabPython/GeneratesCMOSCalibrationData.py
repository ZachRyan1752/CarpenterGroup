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
import math
from scipy.fft import fft2, fftshift
import cv2
import ImageProcessing as IP
import lxml.etree as etree
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup
import time
import matplotlib.pyplot as plt
import tkinter
from tkinter import filedialog

# defs
def LoadMovieAndStitchMovie(MoviePath, MovieStitched):
    ## Append a movie to a list of frames
    ## Inputs:
    ## MoviePath: a paths that will be searched and stiched together i.e. String1
    ## MovieStitched: A list containing individual frames that will be appended to the ouput i.e. [Frame1, Frame2, Frame3, ...]

    Movie = tifffile.imread(MoviePath) ## Reading the movie and loading it as a np.array
    for Frames in Movie: ## For every frame in the loaded movie
        MovieStitched.append(Frames) ## Append the frames to the stitched movie
        
    return MovieStitched ## Return the stitched movie

def LinearFitByColumn(x,XData,YData):
    ## Provide a linear fit of data
    Slope, Intercept, R2, P, STD_ERR = linregress(XData,YData)
    return x, Slope, Intercept, R2 ## Return the fit parameters


def FitColumnLinear(Column, LengthOfMovieArray, AverageData, VarianceData):
    ## This function splits a column into pixels for linear fitting
    ColumnIterator = 0 
    ## Initializing some arrays to store data
    RowSlope = [] ## Stores the linear fit slope (gain) of each row of data processed 
    RowIntercept = [] ## The offset of the fits, representing all kinds of noise in our data
    RowR2 = [] ## The goodness of the fit, from the least-squares fitting

    ### Index changed here as well
    while ColumnIterator < np.shape(AverageData)[1]: ## While the current column index is less then the number of total columns
        PixelVarianceStack = VarianceData[0:LengthOfMovieArray,ColumnIterator] ## Slice the variance at all illumination levels at the index of the column, generating a list of variances at illumination levels for a certain pixel
        PixelAverageStack = AverageData[0:LengthOfMovieArray,ColumnIterator] ## Slice the average at all illumination levels at the index of the column, generating a list of averages at illumination levels for a certain pixel
                
        Column, Slope, Intercept, R2 = LinearFitByColumn(Column, PixelAverageStack, PixelVarianceStack) ## Send the data for fitting at the pixel
        

        ## Append the fit parameters to the data lists for later saving
        RowSlope.append(Slope) 
        RowR2.append(R2)
        RowIntercept.append(Intercept)
        ColumnIterator += 1

    ## Convert the lists to np.array, this prevents the bracket characters from being saved along side real data
    RowSlope = np.array(RowSlope)
    RowIntercept = np.array(RowIntercept)
    RowR2 = np.array(RowR2)
    
    ## Output which column is finished, for debugging and watching the code progress
    print("Finished Column: %s" % Column)

    ## Return the fit parameters from each pixel as a list
    return Column, RowSlope, RowIntercept, RowR2

def FitColumnLinearWrapper(args):
    ## This is a wrapper function, it allows the parrallelization process to pass multiple variables
    ## It could also just be put in the main function but that will be done later
    Column, RowSlope, RowIntercept, RowR2 = FitColumnLinear(args[0], args[1], args[2], args[3])

    ## Return the data back to the main code
    return Column, RowSlope, RowIntercept, RowR2


def GetMovieAveragesFromFilePaths(FolderPath,Paths):
    ## Take a list of paths, collect the movies out of them, and find their averages, returning an array

    LengthOfPaths = len(Paths) ## How many movies there are
    AverageMovieArray = [] ## A list containing the averages of ever movie
    MovieIndex = 0 ## The movie currently being examined
    for Path in Paths: ## For every path in the Paths list
        MovieData = tifffile.imread(FolderPath+Path) ## Read the movie
        MovieAverage = np.average(MovieData) ## Find its average
        MovieShape = np.shape(MovieData) ## Find its shape
        AverageMovieArray.append(MovieAverage) ## Append its average to the list of movie averages
        MovieIndex += 1 ## Progress to the next movie
        print("Averaging Movie %s out of %s" % (MovieIndex, LengthOfPaths)) ## Output how many movies are left

    ## Return the list of movie averages, and the movie shape, for later use
    return AverageMovieArray, MovieShape

def StitchMoviesAndReturnAverageAndVarianceArrays(FolderPath,Paths,PathIndexes,AverageValueArray,MovieShape):
    ## Take a list of movies, stitch them together, average them, and find their variances
    StitchedMovie = [] ## A list containing the individual frames of the movie
    PathIndex = 0 ## Iterates for every movie, mainly used for debugging and code monitoring
    
    for Indicies in PathIndexes: ## For every indice (of the list of paths)
        StitchedMovie = LoadMovieAndStitchMovie(FolderPath+Paths[Indicies],StitchedMovie) ## Send the path to stitch all of the movies together
        PathIndex += 1 ## Increment the path index
        print("Stitching Movie %s out of %s" % (PathIndex, len(PathIndexes))) ## Output which movie is being processed
        
    StitchedMovie = np.array(StitchedMovie) ## Convert the stitched movie to a np.array

    
    StitchedMovieAverage = np.average(StitchedMovie, axis = 0) ## Calculate the average, through the time axis
    StitchedMovieVariance = np.var(StitchedMovie, axis = 0) ## Calculate the variance, through the time axis
    
    ## Return the average and variance frames
    return StitchedMovieAverage, StitchedMovieVariance

def FrameABSubtraction(Frames, FrameAAverage, FrameCounter, MovieShape):
    ## Take two frames, A and B, normalize them against eachother, then subtract A from B
    RatioAB = FrameAAverage / np.average(Frames) ## Find their ratios
    FrameBAdjusted = Frames * RatioAB ## Normalize the frames means against one another
    #FrameBAdjusted = np.expand_dims(FrameBAdjusted, axis = 0)
    print("Starting Frame %s out of %s" % (FrameCounter + 1, MovieShape[0])) ## Output some debug information for the user

    ## Remove the subtracted frame
    return FrameBAdjusted

def FrameABSubtractionWrapper(args):
    ## Wrapper for the above function, will later be integrated into the main function, used to pass multiple inputs to the main function
    ## https://stackoverflow.com/questions/6785226/pass-multiple-parameters-to-concurrent-futures-executor-map
    FrameBAdjusted = FrameABSubtraction(args[0],args[1],args[2],args[3])

    ## Return the frame
    return FrameBAdjusted

def StitchMoviesAndReturnAverageAndVarianceArraysAfterBackgroundSubtractionAndABRatio(FolderPath,Paths,PathIndexes,AverageValueArray,MovieShape,BackgroundAverage):
    StitchedMovie = [] ## List containing the frames of each movie
    PathIndex = 0 ## Which index is currently being examined

    for Indicies in PathIndexes: ## For every index listed in the input
        StitchedMovie = LoadMovieAndStitchMovie(FolderPath+Paths[Indicies],StitchedMovie) ## Load and stitch the movies together
        
        PathIndex += 1 ## Increment
        print("Stitching Movie %s out of %s" % (PathIndex, len(PathIndexes))) ## Output progress to user
        
    StitchedMovie = np.array(StitchedMovie) ## Convert to a np.array
    StitchedMovie = StitchedMovie - BackgroundAverage ## Subtract the background average (offset noise in this case)
    MovieShape = np.shape(StitchedMovie) ## Grab the shape of the movie

    print(MovieShape,PathIndexes)

    ## Get the ratio of frames and adjust them to be the same mean
    ## Subtract Frame B from Frame A
    ## Get the variance of the noise
    #StitchedMovieAB = np.empty((0,MovieShape[1],MovieShape[2]))
    StitchedMovieAB = [] ## For every frame in the stitched movie do an AB mean normalization and subtraction to return only the noise
    with concurrent.futures.ProcessPoolExecutor() as executor: 
        ## This is the hard limit for workers on windows: 
        ## See: https://github.com/python/cpython/issues/71090
        ## See: https://github.com/psf/black/issues/564
        FrameCounter = 0 ## Which frame are we on
        DataPool = [] ## The list of data to be sent to the parralizer
        for Frames in StitchedMovie: ## For every frame in the stitched movie
            if FrameCounter == 0: ## For the first frame, designate it FrameA
                FrameA = Frames ## FrameA is the current frame
                FrameAAverage = np.average(FrameA) ## Find the average of FrameA
            
            else: ## For every other frame:
                ProcessData = [Frames,FrameAAverage, FrameCounter, MovieShape] ## Create a list containing the data to process
                DataPool.append(ProcessData) ## Append that list to the bulk data pool to be processed

            FrameCounter += 1 ## Increment to the next frame

        ## https://stackoverflow.com/questions/6785226/pass-multiple-parameters-to-concurrent-futures-executor-map
        ThreadReturnedData = executor.map(FrameABSubtractionWrapper, DataPool) ## This parrallelizes the data processing, for every element in the datapool, a worker is created that processes the data independently
        ## In windows, this is maximized to 61 concurrent processes, and each process is non blocking to the others, bypassing the global interlock
        ## The .map functionality means that the output of this function will be based on the order of input data, not the order of completion of the processing, faciliting later data analysis

        FinishedProcessCounter = 0 ## How many threads have finished
        for FrameBAdjusted in ThreadReturnedData: ## Grab the adjusted frames from the output of the parrallizer
            StitchedMovieAB.append(FrameBAdjusted) ## Append them to the list for stitching
            FinishedProcessCounter += 1 ## Increment
            print("Finished AB Frame %s out of %s" % (FinishedProcessCounter, np.shape(StitchedMovie)[0])) ## Output debug information to the user

    StitchedMovieAB = np.array(StitchedMovieAB) ## Convert it to a np.array
    print("Stitched Array Shape: %s %s %s" % np.shape(StitchedMovieAB)) ## Output the stitched movie size, this should be # of frames in - 1 (The FrameA is NOT included)

    StitchedMovieAverage = np.average(StitchedMovieAB, axis = 0) ## Find the average through time
    StitchedMovieVariance = np.var(StitchedMovieAB, axis = 0) ## Find the variance through time
    
    return StitchedMovieAverage, StitchedMovieVariance ## Return the average and variance frames


def GainCalibrationFromMovies(IlluminationLevels, **kwargs):
    DataType = kwargs.get("Data","Search") ## Extracting parameters from the function inputs
    DataPath = kwargs.get("DataPath","None") ## Folder path if provided

    if DataType == "Find": ## If you want to find the data in a directory by passing it a string
        FolderPath = DataPath.replace("\n  ","") ## This is formatting from extraneous data stored in the XML, will later be removed
        FolderPath = FolderPath.replace("\n ","") ## This is formatting from extraneous data stored in the XML, will later be removed
        
        DataPath = FolderPath + "Calibration/GainCalibration/"
        ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"

    if DataType == "Search": ## If you want to open file explorer and select a folder like that
        ## https://stackoverflow.com/questions/66663179/how-to-use-windows-file-explorer-to-select-and-return-a-directory-using-python
        tkinter.Tk().withdraw() ## No empty tkinter windows
        
        print("What is the data directory?") 
        DataPath = filedialog.askdirectory() + "/" ## Opens file explorer, user selects a directory
        
        print("Where do you want the output to go?")
        ExtractedDataPath = filedialog.askdirectory()  + "/" ## Opens file explorer, user selects a directory

    StartTime = time.time()

    ## Generate a image containing the gain at each individual pixel
    ## Procedure from:
    ## https://mirametrics.com/tech_note_ccdgain.php
    ## Essentially the same as 2-frame calibration, but on a per-pixel basis, this removes one of the two offsets for the fitting

    MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif") ## Scan the above directory for a list of .tif files, agnostic to file name


    LengthMovieList = len(MoviePaths) ## Determine how many movies are there
    
    AverageMovieValueArray, MovieShape = GetMovieAveragesFromFilePaths(DataPath, MoviePaths) ## Scan and collect averages from every movie in the path

    DarkImageIndex = 0 ## The index currently being scanned
    MinimumAverageValue = np.min(AverageMovieValueArray) ## Find the lowest average of any movie
    IndexArray = [] ## An array of indexes, pointing to movies corresponding to the lowest illumination levels
    for Average in AverageMovieValueArray: ## For every average in the movie
        if Average > MinimumAverageValue * 0.95 and MinimumAverageValue * 1.05 > Average : ## If that average is within 5% of the average of the lowest illumination level
            IndexArray.append(DarkImageIndex) ## Append the index of the array of movies considered to be extensions of the "Dark Movie"
        DarkImageIndex += 1 ## Increment

    print("Getting Dark Movie")
    DarkFrameAverage, DarkFrameVariance = StitchMoviesAndReturnAverageAndVarianceArrays(DataPath,MoviePaths,IndexArray,AverageMovieValueArray,MovieShape) ## Extract the average and variance frame from the dark frame
    
    AverageFrameOverIlluminationLevels = [] ## An array containing the averages at every illumination level
    VarianceFrameOverIlluminationLevels = [] ## An array containing the variance at every illumination level

    AverageFrameOverIlluminationLevels.append(np.array(DarkFrameAverage-DarkFrameAverage)) ## Append the dark frame average, minus itself (to reflect there is no signal contribution caused by light)
    VarianceFrameOverIlluminationLevels.append(np.array(DarkFrameVariance)) ## Append the dark frame variance to the array of variance frames

    RemovedCount = 0 ## Index counter for removal
    for Indexes in IndexArray[::-1]: ## Reverse the list of indexes to be removed
        ## The index list is reversed, as if you were to remove lower indexes first, the rest of the indexes would need to be shifted accordingly, as they would shift as well
        ## We start removing from the END of the index list (highest values first)
        print(AverageMovieValueArray[Indexes],MoviePaths[Indexes],MinimumAverageValue)  ## Debug information --> Are we removing what we think we are removing
        ## This prints out the movie average at the index to be removed, the path at that index, and the minimum average, the first and third values should be identical, the movie path indicated should have an average as indicated
        AverageMovieValueArray.pop(Indexes) ## Remove the average from the list of averages
        MoviePaths.pop(Indexes) ## Remove the path from the list of paths, this prevents double checking of already processed data
        RemovedCount += 1 ## Increment removed counter

    print("Getting Movie Averages and Variances")
    print(AverageMovieValueArray) ## Print the list of averages

    while len(AverageMovieValueArray) > 0: ## While there are averages left in the average array
        MinimumAverageValue = np.min(AverageMovieValueArray) ## Find the minimum (next highest illumination level)
        IndexArray = [] ## Array of indexes corresponding to that illumination level
        ImageIndex = 0 ## The index
        for Average in AverageMovieValueArray: ## For every average in the average movie array
            if Average > MinimumAverageValue * 0.95 and MinimumAverageValue * 1.05 > Average : ## If the average is within 5% of the lowest average
                print(Average, MinimumAverageValue, ImageIndex) ## Print some debug information
                IndexArray.append(ImageIndex) ## Append the index corresponding to the array
            ImageIndex += 1 ## Increment the index

        FrameAverage, FrameVariance = StitchMoviesAndReturnAverageAndVarianceArraysAfterBackgroundSubtractionAndABRatio(DataPath,MoviePaths,IndexArray,AverageMovieValueArray,MovieShape,DarkFrameAverage) ## Get the average and variance of the movie

        RemovedCount = 0 ## This isnt in use and can probably be removed at some point
        for Indexes in IndexArray[::-1]: ## Reverse the list, for reasons stated above
            print(AverageMovieValueArray[Indexes],MoviePaths[Indexes],MinimumAverageValue) 
            AverageMovieValueArray.pop(Indexes)
            MoviePaths.pop(Indexes)
            RemovedCount += 1
        
        AverageFrameOverIlluminationLevels.append(np.array(FrameAverage)) ## Append the average frame to the list of average frames
        VarianceFrameOverIlluminationLevels.append(np.array(FrameVariance)) ## Append the variance frame to the list of variance frames


    AverageFrameOverIlluminationLevels2 = np.empty((0,np.shape(AverageFrameOverIlluminationLevels)[1],np.shape(AverageFrameOverIlluminationLevels)[2])) ## Initialize a new numpy array with the dimensions of the list of the average frames

    for Frames in AverageFrameOverIlluminationLevels: ## For every illumination level
        print(np.average(Frames))
        Frames = np.expand_dims(Frames, axis = 0) ## Expand dims, for stacking
        AverageFrameOverIlluminationLevels2 = np.vstack((AverageFrameOverIlluminationLevels2, Frames)) ## Stack the movies
        ## Note: I have tried using a simple np.array(Data) method, but this appears to output an array with incorrect axes

    VarianceFrameOverIlluminationLevels2 = np.empty((0,np.shape(VarianceFrameOverIlluminationLevels)[1],np.shape(VarianceFrameOverIlluminationLevels)[2])) ## Same as above
    for Frames in VarianceFrameOverIlluminationLevels:
        Frames = np.expand_dims(Frames, axis = 0)
        VarianceFrameOverIlluminationLevels2 = np.vstack((VarianceFrameOverIlluminationLevels2, Frames))

    AverageFrameOverIlluminationLevels = AverageFrameOverIlluminationLevels2 ## Resetting the variables to their original designation
    VarianceFrameOverIlluminationLevels = VarianceFrameOverIlluminationLevels2
    
    print("Average and Variance Array Shapes") 
    print(np.shape(AverageFrameOverIlluminationLevels),np.shape(VarianceFrameOverIlluminationLevels)) ## Print the shapes of the average and variance frame stacks, they should be equal to (# of illumination levels, # pixels X, # pixels Y)
    Iterator1 = 0 
    ColumnSlope = [] ## Slope of each column of pixels
    ColumnIntercept = [] ## Noise of each column of pixels
    ColumnR2 = [] ## R2 of fit of each column of pixels
        
    BigDataArray = []

    AverageDataToProcessPool = [] ## Pool of data to chuck into multithreading
    VarianceDataToProcessPool = [] ## Pool of data to chuck into multithreading
    
    ### index changed here
    while Iterator1 < np.shape(AverageFrameOverIlluminationLevels)[1]: ## While we are less then the number of columns
        AverageData = AverageFrameOverIlluminationLevels[0:IlluminationLevels+1,Iterator1] ## Preslicing the arrays, for (theoretically) faster code execution
        VarianceData = VarianceFrameOverIlluminationLevels[0:IlluminationLevels+1,Iterator1] 

        AverageDataToProcessPool.append(AverageData) ## Append data to the data pool to be processed
        VarianceDataToProcessPool.append(VarianceData) ## Append data to the data pool to be processed
        print("Working on Column Number: %s" % Iterator1) ## Outputting debug information for user
        Iterator1 += 1 ## Incrementing
        
    with concurrent.futures.ProcessPoolExecutor() as executor: 
        ## This is the hard limit for workers on windows: 
        ## See: https://github.com/python/cpython/issues/71090
        ## See: https://github.com/psf/black/issues/564
        ColumnIterator = 0 
        DataPool = [] ## The pool of data to be multiprocessed
        for Data in AverageDataToProcessPool: ## For every piece of data within the average column slices
            Column = ColumnIterator
            AverageData = AverageDataToProcessPool[ColumnIterator] 
            VarianceData = VarianceDataToProcessPool[ColumnIterator]

            ProcessData = [ColumnIterator, IlluminationLevels, VarianceData, AverageData] ## Append them to the data pool
            DataPool.append(ProcessData)
            ColumnIterator += 1 ## Iterate

        ## https://stackoverflow.com/questions/6785226/pass-multiple-parameters-to-concurrent-futures-executor-map
        ThreadReturnedData = executor.map(FitColumnLinearWrapper, DataPool) ## Process the data
        FinishedProcessCounter = 0 
        for Data in ThreadReturnedData: ## For ordered returned data from multiprocessing  
            Column = Data[0]  ## Column number
            RowSlope = Data[1] ## Gain of each pixels in the form of a list
            RowIntercept = Data[2] ## Intercept, reflecting noise of each pixels in the form of a list
            RowR2 = Data[3] ## R2, formatted as above
            BigDataArray.append([Column, RowSlope, RowIntercept, RowR2]) ## Chuck all of the data into a massive list
            FinishedProcessCounter += 1 ## Increment
            print("Finished Column %s out of %s" % (FinishedProcessCounter, MovieShape[1])) ## Output information to the user

    ## Sorting the data
    ## This will be used when I implement asynchronous code execution
    ## This actually shouldnt be needed and will be removed later
    ## Basically: it sorts data so that the column in --> column out

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
 
    SlopeImage = np.array(ColumnSlope) ## Converting the slices of column data into a frame, that can be easily applied to future collected movies
    InterceptImage = np.array(ColumnIntercept)
    R2Image = np.array(ColumnR2)

    VxG = np.sqrt(VarianceFrameOverIlluminationLevels[0]) * SlopeImage ## A frame of the variance at no illumination (readout noise) times the gain, converts it into readout noise in electrons

    VarianceFrameOverIlluminationLevels = VarianceFrameOverIlluminationLevels.astype(np.float32) ## Converting the frames to data type of float32, as it can be read by windows and imageJ
    AverageFrameOverIlluminationLevels = AverageFrameOverIlluminationLevels.astype(np.float32)
    print(np.shape(VarianceFrameOverIlluminationLevels),np.shape(AverageFrameOverIlluminationLevels))

    SlopeImage = SlopeImage.astype(np.float32)
    InterceptImage = InterceptImage.astype(np.float32)
    R2Image = R2Image.astype(np.float32)

    ## Saving the collected data to the output path
    tifffile.imwrite(ExtractedDataPath+"ElectronReadoutNoise.tif",VxG.astype(np.float32))

    tifffile.imwrite(ExtractedDataPath+"GainPixelVariance.tif",VarianceFrameOverIlluminationLevels)
    tifffile.imwrite(ExtractedDataPath+"GainPixelAverage.tif",AverageFrameOverIlluminationLevels)
    tifffile.imwrite(ExtractedDataPath+"GainPixelSlope.tif",SlopeImage)
    tifffile.imwrite(ExtractedDataPath+"GainPixelIntercept.tif",InterceptImage)
    tifffile.imwrite(ExtractedDataPath+"GainPixelR2.tif",R2Image)
    
    ## Printing some more debug information
    print("Average gain: %s" % np.average(SlopeImage))
    print("Average fit R2: %s" % np.average(R2Image))
    StopTime = time.time()
    print("Time To Execute Gain Calibration: %s Seconds" % (StopTime - StartTime))

def SpatialCalibrationFromRuling(LPM, CameraPixelSizeum, **kwargs):
    ## Get the spatial calibration information of a imaging system via a Ronchi ruling
    DataType = kwargs.get("Data","Search") ## Extracting parameters from the function inputs
    DataPath = kwargs.get("DataPath","None") ## Folder path if provided

    if DataType == "Find": ## If you want to find the data in a directory by passing it a string
        FolderPath = DataPath.replace("\n  ","") ## This is formatting from extraneous data stored in the XML, will later be removed
        FolderPath = FolderPath.replace("\n ","") ## This is formatting from extraneous data stored in the XML, will later be removed
        
        DataPath = FolderPath + "Calibration/Resolution/"
        ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"

    if DataType == "Search": ## If you want to open file explorer and select a folder like that
        ## https://stackoverflow.com/questions/66663179/how-to-use-windows-file-explorer-to-select-and-return-a-directory-using-python
        
        print("What is the data directory?") 
        DataPath = filedialog.askdirectory() + "/" ## Opens file explorer, user selects a directory
        
        print("Where do you want the output to go?")
        ExtractedDataPath = filedialog.askdirectory()  + "/" ## Opens file explorer, user selects a directory

    MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif") ## Find all of the .tif images in a directory

    for MoviePath in MoviePaths: ## Read the data at each point
        Movie = tifffile.imread(DataPath + MoviePath)

    ## https://www.geeksforgeeks.org/how-to-find-the-fourier-transform-of-an-image-using-opencv-python/
    FFT = cv2.dft(np.float32(Movie), flags = cv2.DFT_COMPLEX_OUTPUT) ## Calculate the 2D FFT
    FFTCentered = np.fft.fftshift(FFT) ## Center it so that 0 frequency is at the center of the image

    FFTMagnitude = 20*np.log(cv2.magnitude(FFTCentered[:,:,0],FFTCentered[:,:,1])) ## Converting to a magnitude spectrum, with no imaginary component

    FFTNormalized = cv2.normalize(FFTMagnitude, None, 0, 255, cv2.NORM_MINMAX, cv2.CV_8UC1) ## Normalizing the FFT
    
    tifffile.imwrite(ExtractedDataPath+"RonchiRulingFFT.tif",FFTNormalized.astype(np.float32)) ## Saving the output of the FFT
    
    shapeFFT = np.shape(FFTNormalized)[0] ## Getting the shape, cropping it to a small region around the center
    shapeFFTCrop = int(shapeFFT / 2 - 20)
    EditedFFTPeaksPicked, ArrayOfPeaks = IP.FindPeaks(FFTNormalized,8,200,shapeFFTCrop) ## Picking peaks

    tifffile.imwrite(ExtractedDataPath+"FFTPickedPeaks.tif",EditedFFTPeaksPicked.astype(np.float32))

    CenterPeak = ArrayOfPeaks[0] ## Center should always be the largest peak (no spatial frequency)
    Peak1 = ArrayOfPeaks[1]
    Peak2 = ArrayOfPeaks[2]

    CenterCoords = [CenterPeak[1],CenterPeak[2]] ## Converting to a more easily processable format
    Peak1Coords = [Peak1[1],Peak1[2]]
    Peak2Coords = [Peak2[1],Peak2[2]]
    Distance1 = GL.Distance(CenterCoords,Peak1Coords) ## Get the frequency space distance between the two peaks
    Distance2 = GL.Distance(CenterCoords,Peak2Coords)

    AveDistance = (Distance1 + Distance2) / 2 / np.shape(Movie)[0] ## Calculating the average distance

    nmPerLine = 1 / LPM * 1000 * 1000 ## Calulating the spacing between lines on the Ronchi ruling
    umPerLine = 1 / LPM * 1000

    PixelSize = [AveDistance * nmPerLine, CameraPixelSizeum / (AveDistance * nmPerLine) * 1000] ## Calculating pixel size and magnification

    ## Outputting the data to a .xml file
    Header = "PixelCalibration"
    SubHeaders = ["PixelSize","Magnification"]
    SubHeadersDataName = ["PixelSize","Magnification"]
    GL.MakeXMLWithHeaderbs4(Header, SubHeaders, PixelSize, ExtractedDataPath, "PixelSize.xml")


def GetBeamProfile(BeamIntensityAtSample, **kwargs):

    DataType = kwargs.get("Data","Search") ## Extracting parameters from the function inputs
    DataPath = kwargs.get("DataPath","None") ## Folder path if provided

    if DataType == "Find": ## If you want to find the data in a directory by passing it a string
        FolderPath = DataPath.replace("\n  ","") ## This is formatting from extraneous data stored in the XML, will later be removed
        FolderPath = FolderPath.replace("\n ","") ## This is formatting from extraneous data stored in the XML, will later be removed
        
        DataPath = FolderPath + "Calibration/BeamProfile/"
        ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"

    if DataType == "Search": ## If you want to open file explorer and select a folder like that
        ## https://stackoverflow.com/questions/66663179/how-to-use-windows-file-explorer-to-select-and-return-a-directory-using-python
        DataPath = GL.AskForDirectory(InputMessage = "What is the data directory?", OutputMessage = "Selected Path: ")
        ExtractedDataPath = GL.AskForDirectory(InputMessage = "Where do you want the output to go?", OutputMessage = "Selected Path: ")

    MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif") ## Scan for all tifs in the directory
    Movie = tifffile.imread(DataPath + MoviePaths[0])

    MovieAveraged = np.average(Movie, axis = 0) ## Averageing the frame stack to get a even beam profile
    tifffile.imwrite(ExtractedDataPath+"AverageBeamProfile.tif",MovieAveraged.astype(np.float32)) ## Saving the average of the beam profile

    PeakMax = np.max(MovieAveraged) ## Finding the most intense spot in the frame
    PeakCoords = np.where(PeakMax == MovieAveraged) ## Finding where that maximum is
    PeakCoords = list(zip(PeakCoords)) ## Converting 

    
    PixelCalibrationFile = ExtractedDataPath+"PixelSize.xml" ## Reading previous configuration data and importing it
    with open(PixelCalibrationFile) as PCFFile:
        PCF = PCFFile.read()
        PixelSizeXML = BeautifulSoup(PCF, 'xml')

    PixelSizeXMLData = PixelSizeXML.find_all('PixelSize')
    for Data in PixelSizeXMLData:
        PixelSize = np.float32(Data.text)

    ## Finding where the image rises above the background on one size to crop the image for faster fitting
    ## Where the image rised above background - 50 pixels
    OriginalMovieShape = np.shape(MovieAveraged)
    FitBoundsXMin = 0
    FitBoundsYMin = 0

    ImageShape = np.shape(MovieAveraged) ## Getting the shape of the averaged frame
    xaxis = np.linspace(0,ImageShape[0],ImageShape[0]) ## Setting up axes for a 2D Gaussian fit
    yaxis = np.linspace(0,ImageShape[1],ImageShape[1])
    xaxis, yaxis = np.meshgrid(xaxis, yaxis)

    RaveledData = MovieAveraged.ravel() ## Flattening the data array 

    InitialGuess = (PeakMax, PeakCoords[0][0][0] - FitBoundsXMin, PeakCoords[1][0][0] - FitBoundsYMin, 60, 60, 0, 100) ## Initial guess for the fit
    ## 1: Maximum area on the frame --> Peak heigh
    ## 2: X position of peak --> Peak X
    ## 3: Y position of peak --> Peak Y
    ## 4: STD X (from a manual fit in imageJ)
    ## 4: STD Y (from a manual fit in imageJ)
    ## 5: Theta --> Assume no rotation
    ## 6: Offet --> 100 from manufacturer specs (minimum signal obtained with zero input light)

    RGLower = [(PeakMax-100)*0.8,0,0,35,35,-0.25,90] ## Setting bounds for the fitting
    RGUpper = [(PeakMax-100)*1.2,2304,2034,70,70,0.25,110]

    ReasonableGuess = (RGLower,RGUpper)

    print(InitialGuess,np.shape(MovieAveraged))

    popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis, yaxis), RaveledData, p0 = InitialGuess, maxfev = 25000, bounds = ReasonableGuess ) ## Doing the fit, with 25k tries
    FittedData = fits.twoD_Gaussian((xaxis, yaxis), *popt).reshape(ImageShape[0],ImageShape[1]) ## Fitting the parameters as a gaussian profile
    R2 = GL.getR2(MovieAveraged, FittedData) ## Calculating the R2
    print(popt)

    Error = MovieAveraged - FittedData
    Error = np.reshape(Error, np.shape(MovieAveraged)) ## Calculating residuals by subtracting the average from the fit

    tifffile.imwrite(ExtractedDataPath+"FitDifference.tif",Error.astype(np.float32)) ## Changing data type to float32
    tifffile.imwrite(ExtractedDataPath+"BeamProfileFit.tif",FittedData.astype(np.float32))

    FittingParametersArray = []
    FittingParametersArray.extend(popt)
    FittingParametersArray[1] += FitBoundsXMin ## Readjusting to real coordinates (fitboundsxmin is zero at this time)
    FittingParametersArray[2] += FitBoundsYMin
    FittingParametersArray.append(R2)
    print(FittingParametersArray)

    FittingParametersArray.extend([OriginalMovieShape[0],OriginalMovieShape[1]]) ## Adding additional information to the parameters array
    print(FittingParametersArray)

    ## Calculating w0, our equivalent to the 1/e^2 distance
    w0Pixels = (FittingParametersArray[3] + FittingParametersArray[4]) / 2 * 2
    w0 = (FittingParametersArray[3] + FittingParametersArray[4]) / 2 * 2 * PixelSize / 1e7
    LaserPeakIntensity = BeamIntensityAtSample / (math.pi / 2 * ((w0) ** 2)) / 1000 ## Value in W per cm^2

    ## Outputting the results to a xml file
    FittingParametersArray.extend([LaserPeakIntensity, w0, w0Pixels, PixelSize])
    Header = "BeamProfileFitParameters"
    SubHeaders = ["Amplitude","PeakX","PeakY","SigmaXPixels","SigmaYPixels","ThetaRadians","Offset","R2","MovieShapeX","MovieShapeY","LaserPeakIntensityWPercmPercm","w0cm","w0Pixels","PixelSizenm"]
    GL.MakeXMLWithHeaderbs4(Header, SubHeaders, FittingParametersArray, ExtractedDataPath, "BeamFitParameters.xml")

def BlankAverage(FolderPath): ## Dont worry about this we wont use it
    BlankMoviesPath = FolderPath + "Calibration/Background/"
    BlankPath = GL.ScanForFilesInPathByTag(BlankMoviesPath,".tif")

    for MoviePath in BlankPath:
        BlankMovie = tifffile.imread(BlankMoviesPath+MoviePath)
    
    BlankAverageFrame = np.average(BlankMovie, axis = 0)

    BlankMovieAveragePath = FolderPath + "Calibration/CalibrationData/BlankAverage.tif"
    tifffile.imwrite(BlankMovieAveragePath,BlankAverageFrame.astype(np.float32))

def ApplyCalibrationsToMovies(**kwargs):
    ## This function applys the calibrations to a selected folder of tiffs

    Calibrations = kwargs.get("Calibrations", ["GainImage", "CropToBeamProfile"])
    GainConstant = kwargs.get("GainConstant", 0)

    DataType = kwargs.get("Data","Search") ## Extracting parameters from the function inputs
    DataPath = kwargs.get("DataPath","None") ## Folder path if provided

    if DataType == "Find": ## If you want to find the data in a directory by passing it a string
        FolderPath = DataPath.replace("\n  ","") ## This is formatting from extraneous data stored in the XML, will later be removed
        FolderPath = FolderPath.replace("\n ","") ## This is formatting from extraneous data stored in the XML, will later be removed
        
        RawDataPath = FolderPath + "RawData/"

        DataPath = FolderPath + "CorrectedCroppedData/"
        ExtractedDataPath = FolderPath + "Calibration/CalibrationData/"
        CalibrationDataPath = FolderPath + "Calibration/CalibrationData/"

    if DataType == "Search": ## If you want to open file explorer and select a folder like that
        RawDataPath = GL.AskForDirectory(InputMessage = "What is the data directory?", OutputMessage = "Selected Path: ")
        ExtractedDataPath = GL.AskForDirectory(InputMessage = "Where do you want the output to go?", OutputMessage = "Selected Path: ")
        CalibrationDataPath = GL.AskForDirectory(InputMessage = "Where is the calibration data?", OutputMessage = "Selected Path: ")

    MoviePaths = GL.ScanForFilesInPathByTag(RawDataPath,".tif")

    for MoviePath in MoviePaths:
        MovieExample = tifffile.imread(RawDataPath + MoviePath)
        break

    MovieShape = np.shape(MovieExample)
    FrameCount = MovieShape[0]
    SizeX = MovieShape[1]
    SizeY = MovieShape[2]
    GainCalibrationImage = tifffile.imread(CalibrationDataPath+"GainPixelSlope.tif")
    SizeXGain = np.shape(GainCalibrationImage)[0]
    SizeYGain = np.shape(GainCalibrationImage)[1]

    ## Assume the cropping (If any) is done centered
    CenterX = int(1/2*SizeXGain)
    CenterY = int(1/2*SizeYGain)

    CropXMin = int(CenterX - 1/2 * SizeX)
    CropXMax = int(CenterX + 1/2 * SizeX)
    CropYMin = int(CenterY - 1/2 * SizeY)
    CropYMax = int(CenterY + 1/2 * SizeY)

    print(CropXMin,CropXMax,CropYMin,CropYMax)
    ## Cropping the gain calibration to the size of the movie
    GainCalibrationImageCropped = GainCalibrationImage[CropXMin:CropXMax,CropYMin:CropYMax]

    BeamCalibrationFile = CalibrationDataPath+"BeamFitParameters.xml"
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
    BeamCropXMin = int(PeakX - 1.5 * w0)
    BeamCropXMax = int(PeakX + 1.5 * w0)
    BeamCropYMin = int(PeakY - 1.5 * w0)
    BeamCropYMax = int(PeakY + 1.5 * w0)

    print(PeakX,PeakY,BeamCropXMin,BeamCropXMax,BeamCropYMin,BeamCropYMax)
    
    MovieCountIterator = 0
    SubHeaders = []
    DataArray = []

    for MoviePath in MoviePaths:
        Movie = tifffile.imread(RawDataPath + MoviePath)

        if "GainImage" in Calibrations:
            Movie = Movie / GainCalibrationImageCropped ## Applying Gain Calibration

        if "GainStatic" in Calibrations: ## Applying a static gain
            Movie = Movie * GainConstant

        if "CropToBeamProfile" in Calibrations: ## Cropping against the profile of the beam
            Movie = Movie[0:FrameCount,BeamCropXMin:BeamCropXMax,BeamCropYMin:BeamCropYMax]
    
        tifffile.imwrite(ExtractedDataPath + MoviePath,Movie.astype(np.float32))
        
        SubHeaders.append("OriginalMovieName")
        SubHeaders.append("NewMovieName")

        DataArray.append(RawDataPath + MoviePath)
        DataArray.append(MoviePath)

        MovieCountIterator += 1


    Header = "MovieTracing"
    GL.MakeXMLWithHeaderbs4(Header, SubHeaders, DataArray,  ExtractedDataPath, "MovieTracing.xml")

def RunGainCalibrationFromMovies(IlluminationLevels, **kwargs): ## Wrapper to enable the use of concurrent.futures based parallelization
    if __name__ != "__main__":
        GainCalibrationFromMovies(IlluminationLevels, **kwargs)