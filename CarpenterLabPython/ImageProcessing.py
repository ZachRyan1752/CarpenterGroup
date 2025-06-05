import numpy as np
import tifffile
import GeneralLibrary as GL
import scipy.optimize as opt
import fits


## Defs
def GetPeaksByFrameMaximumFunction(Frame, BlankFrame, FrameNumber, PeakPickingSettings):
    EditedFrameForPeakPicking = np.zeros(np.shape(Frame))
    FrameToEdit = Frame - BlankFrame
    EditedFrameForPeakPicking += FrameToEdit ## This is important, setting the two equals makes both point to the same data, += makes EditedFrameForPeakPicking reinitiailzae as a separate array, completely unlinked
    
    PeakList = []
    PeaksRemaining = True
    while PeaksRemaining == True:
        FrameMaximum = np.max(EditedFrameForPeakPicking)
        PeakCoordsList = list(zip(np.where(EditedFrameForPeakPicking == FrameMaximum)))
        
        if FrameMaximum > PeakPickingSettings[0]:
            if len(PeakCoordsList[0][0]) == 1:
                PeakCoords = [PeakCoordsList[0][0][0],PeakCoordsList[1][0][0]]
                PeakData = [FrameNumber, FrameMaximum, PeakCoords[0], PeakCoords[1]]
                EditedFrameForPeakPicking[int(PeakData[2]-PeakPickingSettings[1]):int(PeakData[2]+PeakPickingSettings[1]),int(PeakData[3]-PeakPickingSettings[1]):int(PeakData[3]+PeakPickingSettings[1])] = 5 ## Preventing from repicking the same peak

                PeakList.append(PeakData)
            else:
                PeakCoordinateIterator = 0
                for PeakCoords in PeakCoordsList[0][0]:

                    PeakData = [FrameNumber, FrameMaximum, PeakCoordsList[0][0][PeakCoordinateIterator], PeakCoordsList[1][0][PeakCoordinateIterator]]

                    EditedFrameForPeakPicking[int(PeakData[2]-PeakPickingSettings[1]):int(PeakData[2]+PeakPickingSettings[1]),int(PeakData[3]-PeakPickingSettings[1]):int(PeakData[3]+PeakPickingSettings[1])] = 5 ## Preventing from repicking the same peak
                    PeakCoordinateIterator += 1

                    PeakList.append(PeakData)
        
        else:
            PeaksRemaining = False
    return PeakList, EditedFrameForPeakPicking
    ## Output is [Frame #, Peak Maximum, Peak X, Peak Y]



## Classes
class Movie29_05_2025:
    def __init__(self,MovieData,BlankData,DataOutPath,PeakPickingSettings,PeakGroupingSettings,PeakFittingSettings):
        self.MovieData = MovieData
        self.BackgroundBlankAveraged = BlankData
        self.PeakPickingSettings = PeakPickingSettings
        self.PeakGroupingSettings = PeakGroupingSettings
        self.PeakFittingSettings = PeakFittingSettings
        MovieShape = np.shape(self.MovieData)
        self.FrameCount = MovieShape[0]
        self.PixelsX = MovieShape[1]
        self.PixelsY = MovieShape[2]
        self.ExtractedDataPath = DataOutPath

        GL.CreateIfNotExist(DataOutPath)


        ## Debug Options
        self.DebugPeakPicking = True

    def GetPeaks(self):
        ## This is broken into each frame, for parallelization later
        FrameCounterGetPeaks = 0

        self.PickedPeaks = []
        self.MoviePicked = np.empty((0, self.PixelsX, self.PixelsY))
        for FrameData in self.MovieData:
            PeakList, EditedFrame = self.GetPeaksByFrameMaximum(FrameData, self.BackgroundBlankAveraged, FrameCounterGetPeaks)
            
            EditedFrame = np.expand_dims(EditedFrame, axis = 0)
            for PeakData in PeakList:
                self.PickedPeaks.append(PeakData)

            self.MoviePicked = np.vstack((self.MoviePicked,EditedFrame))
            FrameCounterGetPeaks += 1
            print("Picking Peaks at Frame: %s / %s" % (FrameCounterGetPeaks,self.FrameCount))
        
        if self.DebugPeakPicking == True:
            tifffile.imwrite(self.ExtractedDataPath+"/PickedPeaks.tif",self.MoviePicked.astype(np.float32))
            print(len(self.PickedPeaks))


    def GetPeaksByFrameMaximum(self, FrameData, BlankFrameData, FrameNumber):
        ## This calls an external function, so that other libraries (For instance the microscopecontroller library's autofocus function) to use
        PeakList, EditedFrame = GetPeaksByFrameMaximumFunction(FrameData, BlankFrameData, FrameNumber, self.PeakPickingSettings)
        return PeakList, EditedFrame

    def OrganizePeaksThroughTime(self):
        ## Comparing each peak against each other peak, to organize them through time (across frames)
        self.OrganizedPeaks = []

        PickedPeaksLength = len(self.PickedPeaks)
        OrganizePeaksThroughTimeIterator = 0
        for Peak1 in self.PickedPeaks:
            Peak1Array = [Peak1]
            for Peak2 in self.PickedPeaks:
                if np.all(Peak1 != Peak2): ## If they arent the same peak
                    if GL.Distance([Peak1[2],Peak1[3]],[Peak2[2],Peak2[3]]) < self.PeakGroupingSettings[0]: ## If they are within a certain distance of one another
                        Peak1Array.append(Peak2) ## Append it as a peak
                        self.PickedPeaks.remove(Peak2) ## Remove it from the full peak array to not check for it twice

            OrganizePeaksThroughTimeIterator += 1

            print("Organizing Peak: %s out of %s" % (OrganizePeaksThroughTimeIterator,PickedPeaksLength))

            self.PickedPeaks.remove(Peak1) ## Remove the peak that was checked against so it is not checked twice
            self.OrganizedPeaks.append(Peak1Array)
        
        #print(len(self.OrganizedPeaks))

    def FitPeaksOverMovie(self):
        PSF = self.PeakFittingSettings[0]
        BackgroundLevels = np.average(self.MovieData[0,0:PSF,0:PSF])
        BackgroundSum = np.sum(self.MovieData[0,0:PSF,0:PSF])
        self.FittedPeaksOverTime = []
        PeakNumber = 0
        for PeakArrays in self.OrganizedPeaks: ## OrganizedPeaks = [[Peak1Frame1,Peak1Frame2,...],[Peak2Frame1,Peak2Frame2,...],...]
            PeakFitParameters = PeakArrays[0] ## Just taking the first peak fit info (the center should be the same across frames)
            ## PeakFitParameters
            ## Output is [Frame #, Peak Maximum, Peak X, Peak Y]
            
            AreaToFit = self.MovieData[0:self.FrameCount,int(PeakFitParameters[2]-PSF):int(PeakFitParameters[2]+PSF),int(PeakFitParameters[3]-PSF):int(PeakFitParameters[3]+PSF)]
            PeakFittingDataOverTime = self.FitPeaksOverMovieFunction(AreaToFit, PeakFitParameters, BackgroundSum, PSF, BackgroundLevels)
            PeakFittingOverTimeAdjusted = []
            for Data in PeakFittingDataOverTime:
                Data[7] = Data[7] * 65 ## Eventually set this to use data from the calibration files
                Data[8] = Data[8] * 65
                PeakFittingOverTimeAdjusted.append(Data)
            Header = ["Frame","Height","Peak X","Peak Y","Amplitude","Peak X", "Peak Y", "Sigma X (nm)", "Sigma Y (nm)", "Theta (Radians)", "Offset", "R2", "AreaSum5x5", "Photoelectrons per second"]
            GL.WriteCsv2D_Data(PeakFittingOverTimeAdjusted,self.ExtractedDataPath,"Peak%s.csv" % PeakNumber,Header)
            
            PeakNumber += 1
            print("Fitting Peak #%s out of %s" % (PeakNumber,len(self.OrganizedPeaks)))
            self.FittedPeaksOverTime.append(PeakFittingDataOverTime)
    
    def FitPeaksOverMovieFunction(self,AreaToFit, PeakFitParameters, BackgroundSum, PSF, BackgroundLevels):
        ## This could potentiall be fixed with this: https://scipy-cookbook.readthedocs.io/items/FittingData.html
        ## Again this is split as a separate function as it is a very easily parrallizable task
        xaxis = np.linspace(0, 2*PSF, 2*PSF)
        yaxis = np.linspace(0, 2*PSF, 2*PSF)
        xaxis, yaxis = np.meshgrid(xaxis, yaxis)
        #print(BackgroundSum)
        np.sum(AreaToFit[PeakFitParameters[0]])

        PeakFittingDataOverTime = []
        SumOfPSFAtInitialFit = np.sum(AreaToFit[PeakFitParameters[0]]) ## Getting a PSF sized sum of the data
        FrameCounter = 0
        for Frames in AreaToFit: ## For every frame of the movie:
            if np.sum(Frames) > (np.sum(AreaToFit[PeakFitParameters[0]]) * 0.65): ## If the photons collected are different enough from the background try to fit it
                InitialGuess = (PeakFitParameters[1], PSF, PSF, 2, 2, 0, BackgroundLevels)
                FrameDataRaveled = Frames.ravel() ## Raveling is the flattening of an array
                FittingBound = [[PeakFitParameters[1]*0.5,0,0,1,1,-0.5,0],[PeakFitParameters[1]*2,PSF*2,PSF*2,PSF,PSF,0.5,118]]
                try:
                    popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis, yaxis), FrameDataRaveled, p0 = InitialGuess, maxfev = 500, bounds = FittingBound)
                    FittedData = fits.twoD_Gaussian((xaxis, yaxis), *popt)
                    R2 = GL.getR2(FrameDataRaveled,FittedData)

                except:
                    popt = [0,0,0,0,0,0,0]
                    R2 = 0   

                PeakFitParameters[0] = FrameCounter
                popt[1] = popt[1] + PeakFitParameters[2] - PSF ## Adjusting the peak center to real coordinates
                popt[2] = popt[2] + PeakFitParameters[3] - PSF
                if R2 < 0.85:
                    popt = [0,0,0,0,0,0,0]
                    R2 = 0

                AreaSum5x5 = np.sum(Frames[PSF-5:PSF+5,PSF-5:PSF+5])
                PhotoElectronsPerSecond = AreaSum5x5 / 0.1
                FittingData = [*PeakFitParameters,*popt,R2,AreaSum5x5,PhotoElectronsPerSecond]
                PeakFittingDataOverTime.append(FittingData)
                
            FrameCounter += 1
        return PeakFittingDataOverTime






















import numpy as np
import math
import time
import scipy.optimize as opt
import fits
#import fits
import tifffile
import concurrent
import GeneralLibrary as GL

### Functions to add
## A general --> Import all data from a folder
## General --> Create folder if it doesnt already exist function
## Create subclasses for frames, peaks, and peak groups, allow functions on individual objects
## Go through and verify all indexing for peaks, positions, and frames
## Check --> What frame is indexed, what frame is actually checked etc.
## Fix area of peak - background
## Create a graphing library
## Add a queing system for data processing
## --> Save what step in the process and what file it is on.
## --> Allows for resuming of data processing in the event of a failure
## Add a way to check for bistability of peaks (For making good graphs)
## Append the fraction of peaks fitted in each frame
## Why am I even fitting the peaks?
## Create better folder optimization
## Get bulk statistics (Peaks, average persistance percent (frames also), etc.)
## Figure out why many frames are not being fitted
## Set up real time data analysis
## Implement multithreading where it actually matters
## Make more generic functions --> Fit gaussian, rather than the iterative ones I currently have
## Swap to pandas array for indexing of headers

def WriteCsv2D_Data(Data,Path,Filename,Header):
    #if not filename in os.listdir(path):
    with open(Path+Filename,'w') as CsvFile:
        CsvFile.write("")
            
    CsvFile.close()
     
    with open(Path+Filename,'a') as CsvFile:
        for Rows in Header:
            CsvFile.write(str(Rows)+",")
        
        CsvFile.write("\n")   
            
        for Rows in Data:
            for Columns in Rows:
                CsvFile.write(str(Columns)+",")
            CsvFile.write("\n")        
    
    CsvFile.close()
    
def Distance(Coordinates1,Coordinates2):
    distance = math.sqrt((Coordinates2[0]-Coordinates1[0])**2+(Coordinates2[1]-Coordinates1[1])**2)
    return distance

def getR2(rawdata, fitteddata):
    residuals = rawdata - fitteddata
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((rawdata-np.mean(rawdata))**2)
    R2 = 1 - ( ss_res / ss_tot )
    return R2

def FindPeaks(Frame,MinimumDistance,Threshold,EdgeMaskSize): ## Minimum Distance isnt applied correctly here, remove it in the future.
    PeakListSeparate = []
    #PeakList = np.empty((0,3)) ## [Peak Height, Peak X, Peak Y]
    
    MorePeaksDetected = True
    EditedFrame = np.zeros(np.shape(Frame))
    EditedFrame = EditedFrame + Frame
    EditedFrame[0:EdgeMaskSize] = 50
    EditedFrame[EditedFrame.shape[0]-EdgeMaskSize:EditedFrame.shape[0]] = 50
    EditedFrame[0:EditedFrame.shape[0],0:EdgeMaskSize] = 50
    EditedFrame[0:EditedFrame.shape[0],EditedFrame.shape[0]-EdgeMaskSize:EditedFrame.shape[0]] = 50
    
    while MorePeaksDetected:
        if np.max(EditedFrame) > Threshold:
            PeakMax = np.max(EditedFrame)
            PeakCoordsList = np.where(Frame == PeakMax)
            PeakCoordsList = list(zip(PeakCoordsList[0],PeakCoordsList[1]))
        
        else:
            MorePeaksDetected = False
            break
            
        for PeakCoords in PeakCoordsList:
            EditedFrame[PeakCoords[0]-MinimumDistance:PeakCoords[0]+MinimumDistance,PeakCoords[1]-MinimumDistance:PeakCoords[1]+MinimumDistance] = 10
            PeakData = [PeakMax,PeakCoords[0],PeakCoords[1]]
            PeakListSeparate.append(PeakData)
    
    PeakList = np.array(PeakListSeparate)
    return EditedFrame, PeakList
        

def FindPeaksByMaxProjection(Frame, MinimumDistance, Threshold, EdgeMaskSize):
    MaximumFrame = np.max(Frame, axis = 0)
    print(np.shape(MaximumFrame))
    tifffile.imwrite("D:\Microscope Data\ZDRExampleProject\ExtractedData\MovieMax.tif",MaximumFrame)
    AllPeaks = []
    xcoords = 0
    for x in MaximumFrame:
        ycoords = 0
        for y in x:
            if y > Threshold:
                FrameSlice = Frame[0:np.shape(Frame)[0],xcoords,ycoords]
                FrameNumber = list(zip(np.where(FrameSlice == y)))[0][0]
                #print(FrameNumber)
                #time.sleep(5)
                AllPeaks.append([0, y, xcoords, ycoords])
            ycoords = ycoords + 1
        xcoords = xcoords + 1
    #print(len(AllPeaks))
    AllPeaks = np.array(AllPeaks)
    return AllPeaks

def FindPeaksByMovieSum(Frame,MinimumDistance,Threshold,EdgeMaskSize):
    PeakListSeparate = []
    #PeakList = np.empty((0,3)) ## [Peak Height, Peak X, Peak Y]
    Threshold = Threshold * np.shape(Frame)[0]

    MorePeaksDetected = True
    EditedFrame = np.zeros(np.shape(Frame))
    EditedFrame = EditedFrame + Frame
    EditedFrame[0:EdgeMaskSize] = 50
    EditedFrame[EditedFrame.shape[0]-EdgeMaskSize:EditedFrame.shape[0]] = 50
    EditedFrame[0:EditedFrame.shape[0],0:EdgeMaskSize] = 50
    EditedFrame[0:EditedFrame.shape[0],EditedFrame.shape[0]-EdgeMaskSize:EditedFrame.shape[0]] = 50

    MovieSum = np.sum(Frame, axis = 0)
    print(np.shape(MovieSum))

    while MorePeaksDetected:
        if np.max(EditedFrame) > Threshold:
            PeakMax = np.max(EditedFrame)
            PeakCoordsList = np.where(Frame == PeakMax)
            PeakCoordsList = list(zip(PeakCoordsList[0],PeakCoordsList[1]))
        
        else:
            MorePeaksDetected = False
            break
            
        for PeakCoords in PeakCoordsList:
            EditedFrame[PeakCoords[0]-MinimumDistance:PeakCoords[0]+MinimumDistance,PeakCoords[1]-MinimumDistance:PeakCoords[1]+MinimumDistance] = 10
            PeakData = [PeakMax,PeakCoords[0],PeakCoords[1]]
            PeakListSeparate.append(PeakData)
    
    PeakList = np.array(PeakListSeparate)
    return EditedFrame, PeakList

def FilterPeaksByDistance(ArrayOfPeaks,MinimumDistance): ## This would be faster if done in the previous step
    PrimaryIndex = 0
    while PrimaryIndex < ArrayOfPeaks.shape[0]:
        SecondaryIndex = PrimaryIndex ## Algorithm does not recheck already checked peaks
        Peaks = ArrayOfPeaks[PrimaryIndex]
        #print("%s / %s" % (PrimaryIndex,ArrayOfPeaks.shape[0]))
        while SecondaryIndex < ArrayOfPeaks.shape[0]:
            CheckPeaks = ArrayOfPeaks[SecondaryIndex]
            if np.all(Peaks == CheckPeaks) == False:
                distance = Distance((Peaks[2],Peaks[3]),(CheckPeaks[2],CheckPeaks[3]))
                if distance < MinimumDistance and Peaks[0] == CheckPeaks[0]:
                    ArrayOfPeaks = np.delete(ArrayOfPeaks,SecondaryIndex,axis=0)
            SecondaryIndex += 1
            
        PrimaryIndex += 1
            
    return ArrayOfPeaks


def FitPeak(PeakInfo,Frame,FittingProperties,InitialPeakGuesses):
    ViewingWindow = FittingProperties[0]
    
    xaxis = np.linspace(0,ViewingWindow,ViewingWindow)
    yaxis = np.linspace(0,ViewingWindow,ViewingWindow)
    xaxis, yaxis = np.meshgrid(xaxis,yaxis)
    
    #print(PeakInfo)
    RelevantFrame = Frame[int(PeakInfo[2]-ViewingWindow/2):int(PeakInfo[2]+ViewingWindow/2),int(PeakInfo[3]-ViewingWindow/2):int(PeakInfo[3]+ViewingWindow/2)]
    
    RaveledData = RelevantFrame.ravel()
    
    InitialGuess = (PeakInfo[1], *InitialPeakGuesses)
    RGlower = []
    RGupper = []
    
    for parameters in InitialGuess:
        RGlower.append(parameters*0.3)
        RGupper.append(parameters*3)
    
    RGlower[0] = InitialGuess[0]*0.35
    RGlower[3] = 0.25
    RGlower[4] = 0.25
    RGlower[5] = 0
    RGlower[6] = 0
    
    RGupper[0] = InitialGuess[0]*3
    RGupper[3] = 5
    RGupper[4] = 5
    RGupper[5] = 0.5
    RGupper[6] = 1000 
    
    ReasonableGuessses = (RGlower,RGupper)
    
    BasePeak = RelevantFrame
    #print(BasePeak.shape)
    if RaveledData.size < 100:
        popt = 0
        R2 = 0
        FittedData = 0
        PeakFit = [0,0,0,0,0,0,0]
        
    else:
        try:
            popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis,yaxis), RaveledData, p0= InitialGuess, bounds = ReasonableGuessses, maxfev = 10000)
            FittedData = fits.twoD_Gaussian((xaxis,yaxis),*popt).reshape(ViewingWindow,ViewingWindow)
            R2 = getR2(RelevantFrame,FittedData)
            PeakFit = [R2,*popt]
        
        except Exception as Error:
            print(Error)
            popt = 0
            R2 = 0
            FittedData = 0
            PeakFit = [0,0,0,0,0,0,0]
            #print("here")
        
    FittedPeak = FittedData
    
    
    return FittedPeak, PeakFit, PeakInfo, BasePeak

def FastGaussianFit(Frame,PeakInfo,ViewingWindow):
    ## Quick Fit Function for Gaussian Peaks:
    pass


        
    

class Movie:
    def __init__(self,Identifier,Path,Filename,PeakPickingSettings,FittingProperties,FitProperties,PeakGroupingSettings,ConversionFactors,**kwargs):
        self.ForceLoadNonSavedDataTrue = False
        
        self.ForceLoadNonSavedDataTrue = kwargs.get("ForceLoad",False)
        self.ForceLoadNonSavedData = kwargs.get("Movie",0)

        self.ID = Identifier
        self.Path = Path
        self.Filename = Filename
        self.PeakPickingSettings = PeakPickingSettings
        self.FittingProperties = FittingProperties
        self.FitProperties = FitProperties
        self.PeakGroupingSettings = PeakGroupingSettings
        self.ConversionFactors = ConversionFactors
        if self.ForceLoadNonSavedDataTrue == False:
            self.MovieData = tifffile.imread(Filename)
        else:
            self.MovieData = self.ForceLoadNonSavedData
    
        self.MovieLength = np.shape(self.MovieData)[0]
        self.AllPeaks = np.empty([0,0])
        self.FilteredPeaks = np.empty([0,0])
        self.FittedPeaks = np.empty([0,0])
        ## Change these for debug options
        self.DebugFrameRate = 120
        self.DebugPickPeaks = False
        self.DebugSortPeaks = False
        self.DebugFitPeaks = False
        self.DebugOrganizePeaksThroughTime = False 
        ## Dont Touch Theses
        self.FrameNumber = 0
        GL.CreateIfNotExist(self.Path)

    def GetPeaksByTotalArea(self):
        EditedFrame, self.FilteredPeaks = FindPeaksByMovieSum(self.MovieData,self.PeakPickingSettings[0],self.PeakPickingSettings[1],self.PeakPickingSettings[2])
        tifffile.imwrite("%s/PeaksSelected.tif" % (self.Path),EditedFrame)

    def GetPeaks(self):
        self.AllPeaks = np.empty([0,4]) ## [Frame, Peak Height, Peak X, Peak Y]
        self.FrameNumber = 0 ## Starts here to align with what ImageJ says.
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            FrameNumber = 0
            for FramesToPick in self.MovieData:
                ThreadData = {executor.submit(self.GetPeaksFunct,FramesToPick,FrameNumber): FramesToPick}
                FrameNumber += 1
            #ThreadData = {executor.submit(self.GetPeaksFunct,FramesToPick): FramesToPick for FramesToPick in self.MovieData}        
            for FinishedProcesses in concurrent.futures.as_completed(ThreadData):        
                Peet = FinishedProcesses.result()

        #WriteCsv2D_Data(self.AllPeaks,"Data/%s/Peaks/" % self.ID,"PeakDataAllFrames.csv" ,["Frame","Height","Peak X","Peak Y"])
        #print(self.AllPeaks)
    def SortPeakArrayByHeight(self):
        self.AllPeaks = self.AllPeaks[self.AllPeaks[:,2].argsort()[::-1]] ## Sort the numpy array in reverse order by hieght --> Tallest peaks prioritized, regardless of frame number

    def GetPeaksByMaxFunction(self):
        self.AllPeaks = FindPeaksByMaxProjection(self.MovieData,self.PeakPickingSettings[0],self.PeakPickingSettings[1],self.PeakPickingSettings[2] )

    def GetPeaksFunct(self,Frames,FrameNumber):
        
        
        EditedFrame, PeakData = FindPeaks(Frames,self.PeakPickingSettings[0],self.PeakPickingSettings[1],self.PeakPickingSettings[2])
        
        if self.DebugPickPeaks == True:
            if PeakData.size != 0: 
                WriteCsv2D_Data(PeakData,"%s/" % (self.Path),"PeakDataFrame%s.csv" % (FrameNumber + 1),["Height","Peak X","Peak Y"])
            tifffile.imwrite("%s/PeaksSelected%s.tif" % (self.Path,FrameNumber + 1),EditedFrame)
            print("Picking Peaks at Frame:%s " % (FrameNumber + 1))
            time.sleep(1/self.DebugFrameRate)
            
        if PeakData.shape != (0,):
            FramePeakData = np.zeros((PeakData.shape[0],1)) + FrameNumber + 1    
            DataToAllPeaks = np.hstack((FramePeakData,PeakData))
                
            self.AllPeaks = np.vstack((self.AllPeaks,DataToAllPeaks))
            
        
        return "peet"
    
    def FilterPeaks(self):
        self.FilteredPeaks = FilterPeaksByDistance(self.AllPeaks,self.PeakPickingSettings[2])

        if self.DebugSortPeaks == True:
            WriteCsv2D_Data(self.FilteredPeaks,"%s%s/" % (self.Path,self.ID),"FilteredPeakData.csv",["Frame","Height","Peak X","Peak Y"])
    
        #print(self.FilteredPeaks)

    def TruncatePeakList(self):
        #print(np.shape(self.FittedPeakArray))
        self.FittedPeakArray = self.FittedPeakArray[0:5]
        #print(np.shape(self.FittedPeakArray))

    def FitPeaks(self):
        PeakNumber = 0
        
        AllFittedPeaks = np.empty([0,12])
        for Peaks in self.FilteredPeaks:
            FrameNumber = Peaks[0] - 1
            FrameForFitting = self.MovieData[int(FrameNumber)]
            FittedPeak, PeakFit, Peaks, BasePeak = FitPeak(Peaks,FrameForFitting,self.FittingProperties,self.FitProperties)
            if np.all(PeakFit == [0,0,0,0,0,0,0]):
                PeakFitNegatives = np.zeros([8]) - 1000
                PeakFitNegatives = np.hstack((PeakFitNegatives,Peaks))                
                AllFittedPeaks = np.vstack((AllFittedPeaks,PeakFitNegatives))
                
            else:
                PeakFit = np.hstack((PeakFit,Peaks)) 
                AllFittedPeaks = np.vstack((AllFittedPeaks,PeakFit))
            
            if self.DebugFitPeaks == True:
                tifffile.imwrite("%sPeakFittedFrame%sPeak%s.tif" % (self.Path,int(Peaks[0]),PeakNumber),FittedPeak)
                tifffile.imwrite("%sPeakRawDataFrame%sPeak%s.tif" % (self.Path,int(Peaks[0]),PeakNumber),BasePeak)
                PeakNumber += 1            
        
        self.FittedPeaks = AllFittedPeaks
        #WriteCsv2D_Data(self.FittedPeaks,"Data/%s/Peaks/" % self.ID,"FittedPeaks.csv",["R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
    
    def GetFittedPeaks(self):
        return self.FittedPeaks

    def SortData(self):
        pass
    
    def FilterPeakFits(self):
        self.GoodFits = np.empty([0,12])
        self.BadFits = np.empty([0,12])
        for Peaks in self.FittedPeaks:
            if Peaks[0] < self.FittingProperties[1]:
                self.BadFits = np.vstack((self.BadFits,Peaks))
            else:
                self.GoodFits = np.vstack((self.GoodFits,Peaks))


        #WriteCsv2D_Data(self.GoodFits,"Data/%s/Peaks/" % self.ID,"GoodFits.csv",["R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
        #WriteCsv2D_Data(self.BadFits,"Data/%s/Peaks/" % self.ID,"BadFits.csv",["R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
    
    def RefitBadFits(self): ## Apply processes to refit badly fit peaks. Figure out why then reprocess them here.
        IntialGuessSmallPeaks = (10,10,0.8,0.8,0.001,100)
        RefittedBadPeaks = np.empty([0,12])
        AllRefittedPeaks = np.empty([0,12])
        
        for BadFits in self.BadFits:
            FittedPeak, PeakFit, PeakInfo, BasePeak = FitPeak(BadFits[8:12],self.MovieData[int(BadFits[8]-1)],self.FittingProperties,IntialGuessSmallPeaks)
        
            if np.all(PeakFit == [0,0,0,0,0,0,0]):
                PeakFitNegatives = np.zeros([8]) - 1000
                PeakFitNegatives = np.hstack((PeakFitNegatives,PeakInfo))                
                AllRefittedPeaks = np.vstack((AllRefittedPeaks,PeakFitNegatives))
                
            else:
                PeakFit = np.hstack((PeakFit,PeakInfo)) 
                AllRefittedPeaks = np.vstack((AllRefittedPeaks,PeakFit))
            
            if self.DebugFitPeaks == True:
                tifffile.imwrite("%sPeakFittedFrame%sPeak%s.tif" % (self.Path,int(PeakInfo[0]),PeakNumber),FittedPeak)
                tifffile.imwrite("%sPeakRawDataFrame%sPeak%s.tif" % (self.Path,int(PeakInfo[0]),PeakNumber),BasePeak)
                PeakNumber += 1            
        
        self.RefittedBadPeaks = AllRefittedPeaks
        #WriteCsv2D_Data(self.RefittedBadPeaks,"Data/%s/Peaks/" % self.ID,"RefittedBadPeaks.csv",["R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
    
    def OrganizePeaksThroughTime(self):
        PeakNumber = 1
        self.FittedPeakArray = []
        for Peaks in self.GoodFits:
            GroupedFits = np.empty([0,12])
            IndexGroupingA = 0
            for OtherPeaks in self.GoodFits:
                if np.all(OtherPeaks == Peaks) == False and Peaks[8] != OtherPeaks[8]:
                    if Distance((Peaks[10],Peaks[11]),(OtherPeaks[10],OtherPeaks[11])) < self.PeakGroupingSettings[0]:
                        GroupedFits = np.vstack((GroupedFits,OtherPeaks))
                        self.GoodFits = np.delete(self.GoodFits,IndexGroupingA,axis=0)
                        IndexGroupingA -= 1
                
                IndexGroupingA += 1
            PeakNumberArray = np.empty([np.shape(GroupedFits)[0],1]) + PeakNumber
            GroupedFits = np.hstack((PeakNumberArray,GroupedFits))
            PeakNumber += 1
            if GroupedFits.shape != (0,13):
                self.FittedPeakArray.append(GroupedFits)
        
        if self.DebugOrganizePeaksThroughTime == True:
            IndexGroupingB = 1
            for GroupedPeaks in self.FittedPeakArray:
                WriteCsv2D_Data(GroupedPeaks,"%s" % (self.Path),"GroupedPeak%s.csv" % IndexGroupingB,["Peak Number","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
                IndexGroupingB += 1
    
    def FinalizePeaks(self):
        TotalFrames = np.shape(self.MovieData)[0] - 1
        self.FinalPeakInfo = []
        BackgroundNoise = np.average(self.MovieData[0:10]) ## This a bad way of doing this for data that starts with laser active
        for PeakGroupings in self.FittedPeakArray: ## Grouped Peaks
            CompletePeakInfo = np.empty([0,13])
            ExtraInfoArray = np.empty([0,6])
            #0: "Peak Number"
            #1: "R2"
            #2: "Fit Height"
            #3: "Fit Peak X"
            #4: "Fit Peak Y"
            #5: "Sigma X"
            #6: "Sigma Y"
            #7: "Fit Theta"
            #8: "Fit Offset"
            #9: "Frame"
            #10: "Height"
            #11: "Peak X"
            #12: "Peak Y"
            #13: "Photoelectrons Total"
            #14: "Photoelectrons per Second"
            #15: "Area of Peak"
            #16: "Photoelectrons - Background"
            #17: "Photoelectrons per Second - Background"
            #18: "Area of Peak - Background"
            PercentFitted = 0
            FrameNumber = 0
            for Peaks in PeakGroupings: ## Each Peak Item
                while Peaks[9] != FrameNumber and FrameNumber != self.MovieLength: ## If there is no data for a particular frame
                    PeakHeight = self.MovieData[FrameNumber,int(Peaks[11]),int(Peaks[12])]
                    BlankInfo = np.array([Peaks[0],0,0,0,0,0,0,0,0,FrameNumber,PeakHeight,0,0])
                    CompletePeakInfo = np.vstack((CompletePeakInfo,BlankInfo))
                    PE = PeakHeight # Photoelectrons
                    PES = PE / self.ConversionFactors[1] # Photoelectrons / Second
                    PeakArea = np.sum(self.MovieData[FrameNumber,int(Peaks[11]-self.PeakGroupingSettings[1]):int(Peaks[11]+self.PeakGroupingSettings[1]),int(Peaks[12]-self.PeakGroupingSettings[1]):int(Peaks[12]+self.PeakGroupingSettings[1])]) # Peak Area
                    
                    PEBC = PE - BackgroundNoise ## Background corrected Photoelectonrs
                    PESBC = PEBC / self.ConversionFactors[1] ## Background corrected PE/S
                    PeakAreaBC = PeakArea - (self.PeakGroupingSettings[1]*2)**2*BackgroundNoise ## Background corrected Peak Area
                    
                    ExtraInfoArray = np.vstack((ExtraInfoArray,np.array([PE,PES,PeakArea,PEBC,PESBC,PeakAreaBC])))
                    FrameNumber += 1
                
                if Peaks[9] == FrameNumber: ## Appends the info if their is data at a certain frame
                    CompletePeakInfo = np.vstack((CompletePeakInfo,Peaks))
                    
                    PE = Peaks[2] # Photoelectrons
                    PES = PE * self.ConversionFactors[1] # Photoelectrons / Second
                    PeakArea = np.sum(self.MovieData[FrameNumber,int(Peaks[11]-self.PeakGroupingSettings[1]):int(Peaks[11]+self.PeakGroupingSettings[1]),int(Peaks[12]-self.PeakGroupingSettings[1]):int(Peaks[12]+self.PeakGroupingSettings[1])]) # Peak Area
                    
                    PEBC = PE - BackgroundNoise ## Background corrected Photoelectonrs
                    PESBC = PEBC * self.ConversionFactors[1] ## Background corrected PE/S
                    PeakAreaBC = PeakArea - (self.PeakGroupingSettings[1]*2)**2*BackgroundNoise ## Background corrected Peak Area
                    ExtraInfoArray = np.vstack((ExtraInfoArray,np.array([PE,PES,PeakArea,PEBC,PESBC,PeakAreaBC])))
                    FrameNumber += 1
                    
            while FrameNumber < self.MovieLength:
                PeakHeight = self.MovieData[FrameNumber,int(Peaks[11]),int(Peaks[12])]
                BlankInfo = np.array([Peaks[0],0,0,0,0,0,0,0,0,FrameNumber,PeakHeight,0,0])
                CompletePeakInfo = np.vstack((CompletePeakInfo,BlankInfo))
                PE = PeakHeight # Photoelectrons
                PES = PE * self.ConversionFactors[1] # Photoelectrons / Second
                PeakArea = np.sum(self.MovieData[FrameNumber,int(Peaks[11]-self.PeakGroupingSettings[1]):int(Peaks[11]+self.PeakGroupingSettings[1]),int(Peaks[12]-self.PeakGroupingSettings[1]):int(Peaks[12]+self.PeakGroupingSettings[1])]) # Peak Area
                    
                PEBC = PE - BackgroundNoise ## Background corrected Photoelectonrs
                PESBC = PEBC * self.ConversionFactors[1] ## Background corrected PE/S
                PeakAreaBC = PeakArea - (self.PeakGroupingSettings[1]*2)**2*BackgroundNoise ## Background corrected Peak Area
                        
                ExtraInfoArray = np.vstack((ExtraInfoArray,np.array([PE,PES,PeakArea,PEBC,PESBC,PeakAreaBC])))
                FrameNumber += 1
                
            AllPeakData = np.hstack((CompletePeakInfo,ExtraInfoArray))
            self.FinalPeakInfo.append(AllPeakData)
        


        IndexFinalizeA = 1
        BulkPhotonsPerSecond = []
        BulkTotalCollectedPhotons = []
        TotalIndex = 1
        for FinalizedData in self.FinalPeakInfo:
            PercentFitted = 0
            TotalPhotonsCollected = 0
            for items in FinalizedData:
                if items[1] != 0:
                    PercentFitted += 1
                
                TotalPhotonsCollected += items[16]
                
            TotalIndex += 1
            DataToAppend = [TotalIndex,TotalPhotonsCollected]
            BulkTotalCollectedPhotons.append(DataToAppend)

            PercentFitted = PercentFitted / self.MovieLength * 100
            PercentFitted = np.round(PercentFitted,3)
            WriteCsv2D_Data(FinalizedData,"%s/" %self.Path,"GroupedPeakFinal%sFitted%sPercent.csv" % (IndexFinalizeA,PercentFitted),["Peak Number","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y", "Photoelectrons Total", "Photoelectrons per Second", "Area of Peak", "Photoelectrons - Background", "Photoelectrons per Second - Background", "Area of Peak - Background"])
            IndexFinalizeA += 1
            
            WriteCsv2D_Data(BulkTotalCollectedPhotons,"%s/" %self.Path,"TotalCollectedPhotons.csv",["Index","Total Photons Collected - Background"])


    def OrganizePeaksThroughTimeNew(self):
        PeakNumber = 1
        self.FilteredPeaksForFittingArray = []
        for Peaks in self.FilteredPeaks:
            GroupedFits = np.empty([0,4]) # ["Frame","Height","Peak X","Peak Y"]
            IndexGroupingA = 0
            for OtherPeaks in self.FilteredPeaks:
                if np.all(OtherPeaks == Peaks) == False and Peaks[0] != OtherPeaks[0]:
                    if Distance((Peaks[2],Peaks[3]),(OtherPeaks[2],OtherPeaks[3])) < self.PeakGroupingSettings[0]:
                        GroupedFits = np.vstack((GroupedFits,OtherPeaks))
                        self.FilteredPeaks = np.delete(self.FilteredPeaks,IndexGroupingA,axis=0)
                        IndexGroupingA -= 1
                
                IndexGroupingA += 1
            PeakNumberArray = np.empty([np.shape(GroupedFits)[0],1]) + PeakNumber
            GroupedFits = np.hstack((PeakNumberArray,GroupedFits))
            PeakNumber += 1
            if GroupedFits.shape != (0,13):
                self.FilteredPeaksForFittingArray.append(GroupedFits)
        
        if self.DebugOrganizePeaksThroughTime == True:
            IndexGroupingB = 1
            for GroupedPeaks in self.FilteredPeaks:
                WriteCsv2D_Data(GroupedPeaks,self.Path,"GroupedPeak%s.csv" % IndexGroupingB,["Peak Number","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
                IndexGroupingB += 1

        self.FilteredPeaksForFittingArray = self.FilteredPeaks
        self.FittedPeakArray = self.FilteredPeaksForFittingArray

        print(self.FittedPeakArray)

    def FitPeaksOverMovie(self):
        PeakNumber = 0

        AllFittedPeaks = np.empty([0,12])
        #print(self.FittedPeakArray)
        #time.sleep(5)

        for Peaks in self.FittedPeakArray:
            FrameNumber = 0
            for Frames in self.MovieData:
                
                FittedPeak, PeakFit, Peaks, BasePeak = FitPeak(Peaks,Frames,self.FittingProperties,self.FitProperties)
                
                if np.all(PeakFit == [0,0,0,0,0,0,0]):
                    PeakFitNegatives = np.zeros([8]) - 1000
                    PeakFitNegatives = np.hstack((PeakFitNegatives,Peaks))
                    AllFittedPeaks = np.vstack((AllFittedPeaks,PeakFitNegatives))

                else:
                    Peaks[0] = FrameNumber
                    PeakFit = np.hstack((PeakFit,Peaks)) 
                    #print(PeakFit)
                    #time.sleep(5)
                    AllFittedPeaks = np.vstack((AllFittedPeaks,PeakFit))
                
                if self.DebugFitPeaks == True:
                    tifffile.imwrite("%sPeakFittedFrame%sPeak%s.tif" % (self.Path,int(Peaks[0]),PeakNumber),FittedPeak)
                    tifffile.imwrite("%sPeakRawDataFrame%sPeak%s.tif" % (self.Path,int(Peaks[0]),PeakNumber),BasePeak)
                     
                
                FrameNumber += 1
            #print(PeakNumber)
            PeakNumber += 1   
        self.FittedPeaks = AllFittedPeaks

    def GetFinalizedPeakArray(self):
        return self.FinalPeakInfo
    
    def ReturnMovieData(self):
        return self.MovieData
    
    def ReturnPeaks(self):
        return self.AllPeaks
    
    def ReturnID(self):
        print(self.ID)
        return self.ID
    
    def GetPeakArrayShapes(self):
        return self.AllPeaks.shape, self.FilteredPeaks.shape, self.FittedPeaks.shape, self.GoodFits.shape, self.BadFits.shape, self.RefittedBadPeaks.shape