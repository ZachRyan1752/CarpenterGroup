import numpy as np
import tifffile
import GeneralLibrary as GL
import scipy.optimize as opt
import fits
import scipy.ndimage as ndimage
import concurrent

## Defs
def GetPeaksByFrameMaximumFunction(Frame, BlankFrame, FrameNumber, PeakPickingSettings):
    ## A method of getting peaks by:
    ## 1) Find the maximum in the frame
    ## 2) Finding where that maximum is
    ## 3) Saving the peak maximum, frame #, and coords
    ## 4) Blank the data around the peak (set to 0)
    ## 5) Repeat until some threshold is met
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

def GetPeaksByGradient(FrameToPick, FrameNumber, FrameCount, PeakPickingSettings, PeakFittingSettings):
    ## A method of picking peaks via gradients
    ## 1) Generate the gradient of the frame
    ## 2) Find the maximum in the X gradient frame
    ## 3) in a ROI surrounding the X-Maximum:
    ## 3a) Find the X minimum
    ## 3b) Find the Y maxmimum
    ## 3c) Find the Y minimum
    ## 4) Find where these minima and maxima are, then average their positions
    ## Note: By definition a peak is located where the gradient changes signs
    ## But since it is very unlikely to actually see a zero value in our data, the average of minima and maxima positions will suffice
    print("Starting Frame %s out of %s" % (FrameNumber, FrameCount))
    Frame = np.zeros(np.shape(FrameToPick)) + FrameToPick
    #Frame += FrameToPick
    
    ## Denoising the peak
    DenoiseSTD = PeakPickingSettings[0] #1 # 0.45
    Threshold = PeakPickingSettings[1] #7
    ViewWindow = PeakFittingSettings[0] #4
    
    FrameShape = np.shape(Frame)

    ## Denoising the input data by applying a gaussian filter to it
    FrameDenoised = np.zeros(np.shape(Frame)) + ndimage.gaussian_filter(Frame, DenoiseSTD)
    
    ## Calculating the gradient
    GradientXY = np.gradient(FrameDenoised)
    GradientX = GradientXY[0]
    GradientY = GradientXY[1]
    
    ## Setting the edges to a low value so peak picking doesnt occur here
    GradientX[0:ViewWindow] = 0
    GradientX[FrameShape[0]-ViewWindow:] = 0
    GradientX[0:FrameShape[0],0:ViewWindow] = 0
    GradientX[0:FrameShape[0],FrameShape[1]-ViewWindow:] = 0

    GradientY[0:ViewWindow] = 0
    GradientY[FrameShape[0]-ViewWindow:] = 0
    GradientY[0:FrameShape[0],0:ViewWindow] = 0
    GradientY[0:FrameShape[0],FrameShape[1]-ViewWindow:] = 0

    PeakList = []
    PeaksLeft = True
    while PeaksLeft:
        MaximaX = np.max(GradientX)
        
        if MaximaX < Threshold:
            
            PeaksLeft = False
            break
        
        MaxXCoordsList = np.where(MaximaX == GradientX)
        MaxXCoords = [MaxXCoordsList[0][0],MaxXCoordsList[1][0]]
        
        if MaxXCoords[0] < FrameShape[0] - ViewWindow and MaxXCoords[0] > ViewWindow:
            if MaxXCoords[1] < FrameShape[1] - ViewWindow and MaxXCoords[1] > ViewWindow:
                GradXCropped = GradientX[MaxXCoords[0]-ViewWindow:MaxXCoords[0]+ViewWindow,MaxXCoords[1]-ViewWindow:MaxXCoords[1]+ViewWindow]
                GradYCropped = GradientY[MaxXCoords[0]-ViewWindow:MaxXCoords[0]+ViewWindow,MaxXCoords[1]-ViewWindow:MaxXCoords[1]+ViewWindow]

                XOffset = MaxXCoords[0]-ViewWindow
                YOffset = MaxXCoords[1]-ViewWindow
        
                MinimaX = np.min(GradXCropped)
                MinXCoordsList = np.where(MinimaX == GradXCropped)
                MinXCoords = [MinXCoordsList[0][0]+XOffset,MinXCoordsList[1][0]+YOffset]
        
                MaximaY = np.max(GradYCropped)
                MaxYCoordsList = np.where(MaximaY == GradYCropped)
                MaxYCoords = [MaxYCoordsList[0][0]+XOffset,MaxYCoordsList[1][0]+YOffset]
                
                MinimaY = np.min(GradYCropped)
                MinYCoordsList = np.where(MinimaY == GradYCropped)
                MinYCoords = [MinYCoordsList[0][0]+XOffset,MinYCoordsList[1][0]+YOffset]

                AveXCoords = MaxXCoords[0] + MinXCoords[0] + MaxYCoords[0] + MinYCoords[0]
                AveYCoords = MaxXCoords[1] + MinXCoords[1] + MaxYCoords[1] + MinYCoords[1]
                AveCoords = [int(AveXCoords / 4), int(AveYCoords / 4)]
                FrameMax = FrameToPick[AveCoords[0],AveCoords[1]]
                
                GradientX[AveCoords[0]-ViewWindow:AveCoords[0]+ViewWindow,AveCoords[1]-ViewWindow:AveCoords[1]+ViewWindow] = 0
                Frame[AveCoords[0]-ViewWindow:AveCoords[0]+ViewWindow,AveCoords[1]-ViewWindow:AveCoords[1]+ViewWindow] = 0
                
                PeakList.append([FrameNumber, FrameMax, AveCoords[0], AveCoords[1]])
                
            else:
                GradientX[MaxXCoords[0],MaxXCoords[1]] = 10
                Frame[MaxXCoords[0],MaxXCoords[1]] = 10
        else:
            GradientX[MaxXCoords[0],MaxXCoords[1]] = 10
            Frame[MaxXCoords[0],MaxXCoords[1]] = 10
    
    return Frame, PeakList, GradientX

def GetPeaksByGradientWrapper(args):
    ## A wrapper function for the get peaks by gradient method, used for compatibility with concurrent.futures
    ## Will be removed in the future
    Frame, PeakList, GradientX = GetPeaksByGradient(args[0], args[1], args[2], args[3], args[4])
    return Frame, PeakList, GradientX

def OrganizePeaksThroughTime(PickedPeaks, PeakGroupingSettings):
    ## Comparing each peak against each other peak, to organize them through time (across frames)
    OrganizedPeaks = []
    PickedPeaksLength = len(PickedPeaks)
    OrganizePeaksThroughTimeIterator = 0
    for Peak1 in PickedPeaks:
        Peak1Array = [Peak1]
        for Peak2 in PickedPeaks:
            if np.all(Peak1 != Peak2): ## If they arent the same peak
                if GL.Distance([Peak1[2],Peak1[3]],[Peak2[2],Peak2[3]]) < PeakGroupingSettings[0]: ## If they are within a certain distance of one another
                    Peak1Array.append(Peak2) ## Append it as a peak
                    PickedPeaks.remove(Peak2) ## Remove it from the full peak array to not check for it twice

        OrganizePeaksThroughTimeIterator += 1

        print("Organizing Peak: %s out of %s" % (OrganizePeaksThroughTimeIterator,PickedPeaksLength))

        PickedPeaks.remove(Peak1) ## Remove the peak that was checked against so it is not checked twice
        OrganizedPeaks.append(Peak1Array)
        
    return OrganizedPeaks




def FitPeaksOverMovieFunction2(AreaToFit, PeakFitParameters, PSF, BackgroundLevels, PeakNumber, NumbOfPeaks, MinR2):
    print("Starting Fit #%s out of %s" % (PeakNumber + 1,NumbOfPeaks))
    ## This could potentially be fixed with this: https://scipy-cookbook.readthedocs.io/items/FittingData.html
    ## Again this is split as a separate function as it is a very easily parrallizable task

    xaxis = np.linspace(0, 2*PSF, int(2*PSF))
    yaxis = np.linspace(0, 2*PSF, int(2*PSF))
    xaxis, yaxis = np.meshgrid(xaxis, yaxis)  

    PeakXAdjusted = PSF
    PeakYAdjusted = PSF

    AreaSum = np.sum(AreaToFit[PeakFitParameters[0]])
    PeakFittingDataOverTime = []
    SumOfPSFAtInitialFit = np.sum(AreaToFit[PeakFitParameters[0]]) ## Getting a PSF sized sum of the data
    FrameCounter = 0
    for Frames in AreaToFit: ## For every frame of the movie:
        if True:#np.sum(Frames) > (100*PSF*PSF) * 1.2: ## If the photons collected are different enough from the background try to fit it
            
            

            InitialGuess = (np.max(Frames), PSF, PSF, 1.7, 1.7, 0.005, 100)
            FrameDataRaveled = Frames.ravel() ## Raveling is the flattening of an array
            FittingBound = [[np.max(Frames)*0.5,0,0,0.25,0.25,0,0],[np.max(Frames)*2,PSF*2,PSF*2,5,5,0.5,2**12]]
            PeakFitParameters[1] = np.max(Frames)

            try:
                ## Fitting the data, wrapped in a try/except loop as sometimes it picks peaks wrong or there is just no data to fit
                popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis, yaxis), FrameDataRaveled, maxfev = 500, p0 = InitialGuess, bounds = FittingBound)
                FittedData = fits.twoD_Gaussian((xaxis, yaxis), *popt)
                R2 = GL.getR2(FrameDataRaveled,FittedData)

            except Exception as printout:
                print(printout)
                print(FittingBound)
                print(PeakFitParameters)
                print(AreaSum * 0.85, np.sum(Frames))
                popt = [1,1,1,0,0,0,0]
                R2 = 0

            PeakFitParameters[0] = FrameCounter
            popt[1] = popt[1] + PeakFitParameters[2] - PSF ## Adjusting the peak center to real coordinates
            popt[2] = popt[2] + PeakFitParameters[3] - PSF

            if R2 < MinR2:
                pass
                #popt = [0,0,0,0,0,0,0]

            AreaSum5x5 = np.sum(Frames[PSF-5:PSF+5,PSF-5:PSF+5])
            PhotoElectronsPerSecond = AreaSum5x5 / 0.1
            FittingData = [*PeakFitParameters,*popt,R2,AreaSum5x5,PhotoElectronsPerSecond]
            PeakFittingDataOverTime.append(FittingData)

        else:
            popt = [2,2,2,0,0,0,0]
            R2 = 0
            AreaSum5x5 = np.sum(Frames[PSF-5:PSF+5,PSF-5:PSF+5])
            PhotoElectronsPerSecond = AreaSum5x5 / 0.1
            FittingData = [*PeakFitParameters,*popt,R2,AreaSum5x5,PhotoElectronsPerSecond]
            PeakFittingDataOverTime.append(FittingData)           
            ## [Frame #, Peak Height, Peak X, Peak Y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset, R2, Area sum 5x5, Photoelectrons Per Second]

        FrameCounter += 1

    return PeakFittingDataOverTime

def FitPeaksOverMovieFunctionWrapper2(args):
    ## The function must be called as a method of the class explicetely, since it is a static method
    PeakFittingDataOverTime = FitPeaksOverMovieFunction2(args[0], args[1], args[2], args[3], args[4], args[5], args[6])

    return PeakFittingDataOverTime




def FitPeaksOverMovieFunction(AreaToFit, PeakFitParameters, BackgroundSum, PSF, BackgroundLevels, xaxis, yaxis, PeakCounter, FrameCount):
    ## This could potentially be fixed with this: https://scipy-cookbook.readthedocs.io/items/FittingData.html
    ## Again this is split as a separate function as it is a very easily parrallizable task
    
    xaxis = np.linspace(0, 2*PSF, 2*PSF)
    yaxis = np.linspace(0, 2*PSF, 2*PSF)
    xaxis, yaxis = np.meshgrid(xaxis, yaxis)  

    AreaSum = np.sum(AreaToFit[PeakFitParameters[0]])

    PeakFittingDataOverTime = []
    SumOfPSFAtInitialFit = np.sum(AreaToFit[PeakFitParameters[0]]) ## Getting a PSF sized sum of the data
    FrameCounter = 0
    for Frames in AreaToFit: ## For every frame of the movie:
        if np.sum(Frames) > (AreaSum * 0.65): ## If the photons collected are different enough from the background try to fit it
            InitialGuess = (PeakFitParameters[1], PSF, PSF, 2, 2, 0, BackgroundLevels)
            FrameDataRaveled = Frames.ravel() ## Raveling is the flattening of an array
            FittingBound = [[PeakFitParameters[1]*0.5,0,0,0.25,0.25,-0.5,0],[PeakFitParameters[1]*2,PSF*2,PSF*2,PSF,PSF,0.5,118]]
            
            try:
                popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis, yaxis), FrameDataRaveled, maxfev = 10, p0 = InitialGuess, bounds = FittingBound)
                FittedData = fits.twoD_Gaussian((xaxis, yaxis), *popt)
                R2 = GL.getR2(FrameDataRaveled,FittedData)
            except:
                try:
                    popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis, yaxis), FrameDataRaveled, maxfev = 100, p0 = InitialGuess, bounds = FittingBound)
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

        else:
            popt = [0,0,0,0,0,0,0]
            R2 = 0
            AreaSum5x5 = np.sum(Frames[PSF-5:PSF+5,PSF-5:PSF+5])
            PhotoElectronsPerSecond = AreaSum5x5 / 0.1
            FittingData = [*PeakFitParameters,*popt,R2,AreaSum5x5,PhotoElectronsPerSecond]
            PeakFittingDataOverTime.append(FittingData)           
            
            ## [Frame #, Peak Height, Peak X, Peak Y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset, R2, Area sum 5x5, Photoelectrons Per Second]

        FrameCounter += 1


    return PeakFittingDataOverTime


## Classes
class Movie29_05_2025:
    def __init__(self,MovieData,BlankData,DataOutPath,PeakPickingSettings,PeakGroupingSettings,PeakFittingSettings, **kwargs):
        ## Initializing some internal parameters
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
        self.DebugPeakPicking = kwargs.get("DebugPicking",False)

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
            tifffile.imwrite(self.ExtractedDataPath+"PickedPeaks.tif",self.MoviePicked.astype(np.float32))
            print(len(self.PickedPeaks))


    def GetPeaksByGradient(self):
        ## Getting peaks via the gradient method for each frame
        CompletePeakList = []
        EditedPickedMovie = []
        EditedGradientX = []

        with concurrent.futures.ProcessPoolExecutor() as executor: 
            ## This is the hard limit for workers on windows: 
            ## See: https://github.com/python/cpython/issues/71090
            ## See: https://github.com/psf/black/issues/564
            ArrayCounter = 0
            DataPool = []
            for Frames in self.MovieData:
                ProcessData = [Frames, ArrayCounter, self.FrameCount, self.PeakPickingSettings, self.PeakFittingSettings]
                DataPool.append(ProcessData)
                ArrayCounter += 1

            ## https://stackoverflow.com/questions/6785226/pass-multiple-parameters-to-concurrent-futures-executor-map
            ThreadReturnedData = executor.map(GetPeaksByGradientWrapper, DataPool)
            
            FinishedProcessCounter = 0
            for Data in ThreadReturnedData:
                Frame = Data[0]
                PeakList = Data[1]
                GradientX = Data[2]
                
                EditedPickedMovie.append(Frame)
                EditedGradientX.append(GradientX)
                
                for Peaks in PeakList:
                    CompletePeakList.append(Peaks)

                FinishedProcessCounter += 1
                print("Finished Frame %s out of %s" % (FinishedProcessCounter, self.FrameCount))

        EditedPickedMovie = np.array(EditedPickedMovie)
        EditedGradientX = np.array(EditedGradientX)

        if self.DebugPeakPicking == True:
            tifffile.imwrite(self.ExtractedDataPath + "PickedPeaks.tif", EditedPickedMovie.astype(np.float32))   
            tifffile.imwrite(self.ExtractedDataPath + "PickedGradientX.tif", EditedGradientX.astype(np.float32))
            print(len(CompletePeakList))
        
        self.PickedPeaks = CompletePeakList

    def TruncateOrganizedPeaksList(self,Number):
        ## Truncate the peak list after they were organized through time
        ## This is useful for testing data, or the auto focus function, where only a few peaks need to be analyzed to see a clearer picture
        self.OrganizedPeaks = self.OrganizedPeaks[0:Number]


    
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

    def FitPeaksOverMovie(self):
        ## Fitting ever peak picked and organized at every frame in the movie
        PSF = self.PeakFittingSettings[0]  
        
        BackgroundLevels = 100
        BackgroundSum = np.sum(self.MovieData[0,0:PSF,0:PSF])
        self.FittedPeaksOverTime = []
       
        DataPool = []
        PeakNumber = 0
        NumbOfPeaks = len(self.OrganizedPeaks)
        for PeakArrays in self.OrganizedPeaks: ## OrganizedPeaks = [[Peak1Frame1,Peak1Frame2,...],[Peak2Frame1,Peak2Frame2,...],...]
            print("Fitting Peak #%s out of %s" % (PeakNumber + 1,len(self.OrganizedPeaks)))
            PeakFitParameters = PeakArrays[0]
            ## Output is [Frame #, Peak Maximum, Peak X, Peak Y]

            AreaToFit = self.MovieData[0:self.FrameCount,int(PeakFitParameters[2]-PSF):int(PeakFitParameters[2]+PSF),int(PeakFitParameters[3]-PSF):int(PeakFitParameters[3]+PSF)]

            #Data = [AreaToFit, PeakFitParameters, BackgroundSum, PSF, BackgroundLevels, xaxis, yaxis, PeakNumber, self.FrameCount]
            Data = [AreaToFit, PeakFitParameters, PSF, BackgroundLevels, PeakNumber, NumbOfPeaks, self.PeakFittingSettings[1]]

            DataPool.append(Data)
            PeakNumber += 1

        with concurrent.futures.ProcessPoolExecutor() as executor:  ## ProcessPoolExecutor
            ## This is the hard limit for workers on windows: 
            ## See: https://github.com/python/cpython/issues/71090
            ## See: https://github.com/psf/black/issues/564
            ## https://stackoverflow.com/questions/6785226/pass-multiple-parameters-to-concurrent-futures-executor-map
            ThreadReturnedData = executor.map(FitPeaksOverMovieFunctionWrapper2, DataPool)
            FinishedProcessCounter = 0
            for PeakFittingDataOverTime in ThreadReturnedData:  
                PeakFittingOverTimeAdjusted = []
                FitPercent = 0
                for Data in PeakFittingDataOverTime:
                    Data[7] = Data[7] * 65 ## Eventually set this to use data from the calibration files
                    Data[8] = Data[8] * 65
                    PeakFittingOverTimeAdjusted.append(Data)
                    if Data[11] > self.PeakFittingSettings[1]:
                        FitPercent += 1

                FitPercent = np.round(FitPercent / self.FrameCount * 100, 2)
                
                if FitPercent > self.PeakFittingSettings[2]:
                    ## Saving the output of the fitting to a .CSV file
                    Header = ["Frame","Height","Peak X","Peak Y","Amplitude","Peak X Fit (Pixels)", "Peak Y Fit (Pixels)", "Sigma X (nm)", "Sigma Y (nm)", "Theta (Radians)", "Offset", "R2", "AreaSum5x5", "Photoelectrons per second"]
                    GL.WriteCsv2D_Data(PeakFittingOverTimeAdjusted,self.ExtractedDataPath,"Peak%s_%sPercentFitted.csv" % (FinishedProcessCounter,FitPercent),Header)
                
                    print("Finished Peak %s out of %s" % (FinishedProcessCounter + 1, len(self.OrganizedPeaks)))

                else:
                    print("Fitting failed, not enough fits")

                FinishedProcessCounter += 1

                self.FittedPeaksOverTime.append(PeakFittingDataOverTime)

    def GetFittedPeaks(self):
        return self.FittedPeaksOverTime    
    
    def RunGetPeaksByGradient(self):
        ## Wrapper for functions, this enables the use of concurrent futures for asynchronous code execution without recursively calling the entire library thousands of times
        if __name__ != "__main__":
            self.GetPeaksByGradient()
    
    def RunOrganizePeaksThroughTime(self):
        if __name__ != "__main__":
            self.OrganizePeaksThroughTime()

    def RunFitPeaksOverMovie(self):
        if __name__ != "__main__":
            self.FitPeaksOverMovie()