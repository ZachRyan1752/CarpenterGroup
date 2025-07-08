## Import Library
import pymmcore_plus as pymmcore
import tifffile
import os
from datetime import datetime
from dicttoxml import dicttoxml
import numpy as np
import useq
from useq import MDAEvent, MDASequence
#import pycromanager as PM
import time
import pandas as pd
import statsmodels.api as sm
import scipy
import random
import ImageProcessing as IP

## Relevant literature
## https://pymmcore-plus.github.io/pymmcore-plus/api/cmmcoreplus/



## Defs
def CreateIfNotExist(Path):
    if os.path.exists(Path):
        pass
    else:
        os.mkdir(Path)

def GetFormattedYMDHSAP():
    DateAndTime = datetime.now()
    DateAndTimeFormatted = datetime.strftime(DateAndTime,"20%y_%m_%d_%I_%M_%p")
    return DateAndTimeFormatted

def GetFormattedYMD():
    DateAndTime = datetime.now()
    DateAndTimeFormatted = datetime.strftime(DateAndTime,"20%y_%m_%d")
    return DateAndTimeFormatted

def GetFormattedHSAP():
    DateAndTime = datetime.now()
    DateAndTimeFormatted = datetime.strftime(DateAndTime,"%I_%M_%p")
    return DateAndTimeFormatted

## Classes
class ImagingSystem():
    def __init__(self, MicroscopeDeviceAdapterPath, ExtractedFolderPath, ConfigurationPath, Iterator, **kwargs):

        self.MovieNumber = Iterator
        
        CenterStageOnStart = kwargs.get("CenterStageOnStart", False)
        MicroscopeSetupDebug = kwargs.get("MicroscopeSetupDebug", False)

        self.MicroscopeDeviceAdapterPath = MicroscopeDeviceAdapterPath
        self.ExtractedFolderPath = ExtractedFolderPath # + "/RawData/"
        self.ConfigurationPath = ConfigurationPath

        ## Microscope Setup (DO NOT TOUCH)
        self.mmc = pymmcore.CMMCorePlus()  # Instance micromanager core
        self.mmc.setDeviceAdapterSearchPaths(self.MicroscopeDeviceAdapterPath)
        self.mmc.loadSystemConfiguration(self.ConfigurationPath) # Loading the configuration file  

        self.Position = self.GetStagePosition() ## Getting the position of the stage
        self.Stage = self.mmc.getXYStageDevice()
        self.mmc.setCircularBufferMemoryFootprint(16384) ## Size of the circular buffer (Stores images from continuous acquisitions)

        ## Indexes for where in the stage scanning algorithm it is
        self.MoveInPatternIndicies = [0,0,0]

        ## Getting the devices
        CameraDevice = self.mmc.getCameraDevice()
        self.CameraDevice = CameraDevice

        self.ShutterMode = "Toggle"

        if CenterStageOnStart == True:
            self.mmc.setXYPosition(70,70)
            self.mmc.setZPosition(70) ## Starts at the middle so it can focus up or down
            self.Position = [70,70,70]
        
        if MicroscopeSetupDebug == True:
            print(self.mmc.getVersionInfo()) # gives: 'MMCore version 11.5.0'
            print(self.mmc.getAPIVersionInfo()) # gives: 'Device API version 73, Module API version 10'
            print(self.mmc.getLoadedDevices()) # Returns loaded devices
            print(self.mmc.getCameraChannelNames())

        ## Various "Fixes"
        ## https://forum.image.sc/t/micro-manager-become-extremely-slow/104862
        #self.mmc.setChannelGroup("")

        @self.mmc.mda.events.frameReady.connect ## When a MDA frameReady event occurs
        def on_frame(image: np.ndarray, event: useq.MDAEvent):
            if self.FPS != 1:
                time.sleep(self.MovementDelay)
                #print(self.MDApositions)
                self.StartTimens = time.time_ns()
                self.Shutter(ShutterMode = "Set", ShutterState = True)
                Movie = self.CollectMovie(self.ExposureTime, self.FPS, self.Duration, **self.kwargspassthrough) ## Directly pass the movement coordinates to the collect movie function to reduce delay
                self.Shutter(ShutterMode = "Set", ShutterState = False)
                tifffile.imwrite(self.MovieNamesTemp[0],Movie)
                self.MovieNamesTemp.pop(0)

            else:
                time.sleep(self.MovementDelay)
                self.Shutter(ShutterMode = "Set", ShutterState = True)
                Movie = self.CollectImage(self.ExposureTime, **self.kwargspassthrough) ## Directly pass the movement coordinates to the collect movie function to reduce delay
                self.Shutter(ShutterMode = "Set", ShutterState = False)
                tifffile.imwrite(self.MovieNamesTemp[0],Movie)
                self.MovieNamesTemp.pop(0)

    def GetStagePosition(self): ## Saves the new position internally, then outputs a list containing the X, Y, and Z positions of the stage
        XPos = self.mmc.getXPosition()
        YPos = self.mmc.getYPosition()
        ZPos = self.mmc.getZPosition()
        self.Position = [XPos,YPos,ZPos]
        return self.Position

    def MovieIteratorIncrement(self):
        self.MovieNumber += 1

    def OffsetStage(self,SetPosition): ## 
        try:
            Dx = SetPosition
            self.mmc.setRelativeXYZPosition(Dx[0],Dx[1],Dx[2])
            self.GetStagePosition()

        except:
            print(SetPosition,self.GetStagePosition())


    def MoveStageWithErrorCorrection(self,SetPosition):
        OldPosition = self.GetStagePosition()
        
        Delta = [0,0,0]
        Delta[0] = SetPosition[0] - OldPosition[0]
        Delta[1] = SetPosition[1] - OldPosition[1]
        Delta[2] = SetPosition[2] - OldPosition[2]

        self.mmc.setXYPosition(SetPosition[0]+Delta[0],SetPosition[1]+Delta[1])
        self.mmc.setZPosition(SetPosition[2]+Delta[2])
        OldPosition = self.GetStagePosition()


    def MoveStage(self,SetPosition):
        self.mmc.setXYPosition(SetPosition[0],SetPosition[1])
        self.mmc.setZPosition(SetPosition[2])
        OldPosition = self.GetStagePosition()

    def Unload(self):
        self.mmc.unloadAllDevices()

    def GetCameraDevice(self):
        return self.CameraDevice

    def SetInternalCollectionSettings(self,**kwargs):
        self.ExposureTime = kwargs.get("Exposure",100)
        self.FrameSize = kwargs.get("FrameSize", [2304,2304])
        self.CropType = kwargs.get("CropType", "Centered")

    def CaptureImage(self,SaveState,**kwargs):
        ExposureTime = kwargs.get("Exposure",100)
        FrameSize = kwargs.get("FrameSize", [2304,2304])
        CropType = kwargs.get("CropType", "Centered")

        Settings = kwargs.get("SettingsFrom", "External")

        if Settings == "Internal":
            ExposureTime = self.ExposureTime
            FrameSize = self.FrameSize
            CropType = self.CropType

        ## Save State:
        ## 1) "Save"
        ## Saves the image to the listed directory
        ## 2) "Return"
        ## Returns the image and relevant metadata
        ## 3) SaveAndReturn
        ## Does both

        CameraDevice = self.mmc.getCameraDevice()
        self.CameraDevice = CameraDevice
        self.mmc.setCameraDevice(CameraDevice)
        self.mmc.setExposure(CameraDevice,ExposureTime) ## IN MILLISECONDS!!
        #self.mmc.setProperty(CameraDevice,"HamamatsuHam_DCAM-ScanMode",1) ## Fix this eventually
 

        self.mmc.snapImage()
        ImageData = self.mmc.getTaggedImage()

        Time = GetFormattedHSAP()
        
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.ExtractedFolderPath + "_" + Time + "_" + PositionString + "1.tiff"
    
        Image = ImageData[0]
        MetaData = ImageData[1]


        if CropType == "Origin":
            Image = Image[FrameSize[0],FrameSize[1]]
    
        if CropType == "Centered":
            Shape = np.shape(Image)
            CenterX = int(1/2*Shape[0])
            CenterY = int(1/2*Shape[1])

            CropXMin = int(CenterX - 1/2 * FrameSize[0])
            CropXMax = int(CenterX + 1/2 * FrameSize[0])
            CropYMin = int(CenterY - 1/2 * FrameSize[1])
            CropYMax = int(CenterY + 1/2 * FrameSize[1])

            Image = Image[CropXMin:CropXMax,CropYMin:CropYMax]

        
        if SaveState == "Save":
            OME_Xml = dicttoxml(MetaData)
            tifffile.imwrite(ImageName, Image, description = OME_Xml, metadata = None)

        if SaveState == "Return":
            return Image, MetaData
    
        if SaveState == "SaveAndReturn":
            OME_Xml = dicttoxml(MetaData)
            tifffile.imwrite(ImageName, Image, description = OME_Xml, metadata = None)
            return Image, MetaData

    #def RunMDA(self, MDASequence):
    #    self.FrameBuffer = []
    #    @self.mmc.mda.events.frameReady.connect
    #    def on_frame(image: np.ndarray, event: useq.MDAEvent):
    #        self.FrameBuffer.append(image)

    #    self.mmc.run_mda(MDASequence)

    def Shutter(self,**kwargs):
        ShutterMode = kwargs.get("ShutterMode",self.ShutterMode)
        State = kwargs.get("ShutterState", True)

        if ShutterMode == "Auto":
            self.mmc.setAutoShutter(State)

        if ShutterMode == "Set":
            self.mmc.setShutterOpen(State)

    def CollectImage(self, ExposureTime, **kwargs):
        self.mmc.setExposure(ExposureTime) ## IN MILLISECONDS!!
        
        start = time.time()
        self.mmc.snapImage()
        stop = time.time()
        delta = stop - start
        print(delta, 1/delta)

        ImageData = self.mmc.getTaggedImage()

        Image = ImageData[0]
        MetaData = ImageData[1]

        return Image

    def CollectMovie(self, ExposureTime, FPS, Duration, **kwargs):
        
        self.mmc.setExposure(ExposureTime)
        #ROI = kwargs.get("ROI", [0,0,2304,2304])


        #ROI = kwargs.get("ROI", "Default")
        #if ROI != "Default": ## This takes almost 180 ms to execute, do it before we call to start the move to prevent delays
        #    self.mmc.setROI(ROI[0],ROI[1],ROI[2],ROI[3])


        #self.mmc.setExposure(ExposureTime) ## Set exposure time now to reduce delay
        FramesToCollect = Duration * FPS

        #if StepTo != ["X","X","X"]:
        #    self.MoveStage(StepTo)
            

        
        self.EndTimens = time.time_ns()
        print("%s us Movement delay" % ((self.EndTimens - self.StartTimens) / 1000))
        #print(self.EndTimens - self.StartTimens)
        #print(Duration*FPS, ExposureTime)
        self.mmc.startSequenceAcquisition(Duration*FPS, ExposureTime, True)
        start = time.time()
        while self.mmc.isSequenceRunning(): ## https://sourceforge.net/p/micro-manager/mailman/micro-manager-general/thread/03A4DB6AE86963478C9EB15787545DA23E39DC1A%40BE-EXCH-01.andortech.net/
            pass
        
        stop = time.time()
        delta = stop - start
        print(delta, FramesToCollect/delta)
        #self.mmc.stopSequenceAcquisition()
        
        MovieSlices = list()
        while self.mmc.getRemainingImageCount() > 0: ## This needs to be done while acquisition happens or the "circular buffer" will over flow
            if len(MovieSlices) == FramesToCollect:
                break
            else:
                ImageData = self.mmc.popNextImage()
                ImageData = ImageData.astype(np.int32)
                MovieSlices.append(ImageData.astype(np.int32))
        
        AssembledMovie = np.stack(MovieSlices,axis=0)
        #AssembledMovie = AssembledMovie.reshape((FramesToCollect,ROI[2],ROI[3]))
        
        print("Movie Shape: %s" % ([np.shape(AssembledMovie)]))
        
        self.mmc.clearCircularBuffer() ## Remove any old frames that arent from the current movie
        #time.sleep(2)

        return AssembledMovie
        


    def MDAMovieSequence(self, ExposureTime, Duration, FPS, MovieNames, **kwargs):
        CL = self.GetStagePosition()
        CurrentLocation = [(CL[0],CL[1],CL[2])] ## Wrapped tuple, so it can be unpacked correctly

        self.ExposureTime = ExposureTime
        self.Duration = Duration
        self.FPS = FPS

        ROI = kwargs.get("ROI", "Default")
        if ROI != "Default": ## This takes almost 180 ms to execute, do it before we call to start the move to prevent delays
            self.mmc.setROI(ROI[0],ROI[1],ROI[2],ROI[3])
        
        self.MDApositions = kwargs.get("Positions", CurrentLocation)
        self.MovementDelay = kwargs.get("Delay", 0.025)

        #self.MovieNamesTemp = []
        self.MovieNamesTemp = [*MovieNames]

        self.mmc.setExposure(ExposureTime)

        self.kwargspassthrough = kwargs
        DefaultSequence = MDASequence(
            config = {"config": "HammamatsuCam", "exposure": 0.0}, ## Exposure time is set to zero, as we arent actually intending to collect any data at this time
            time_plan = {"interval": 0.00, "loops": len(self.MDApositions)}, ## Collect at every position
            #grid_plan = {'rows': 4, "columns": 4}
            #stage_positions = [*positions]
            #time_plan = {"interval": 1 / FPS, "loops": Duration * FPS},
            )

        Sequence = kwargs.get("Sequence", DefaultSequence)

        #self.mmc.setExposure(0)
        self.mmc.setExposure(ExposureTime) ## Set exposure time now to reduce delay



        
        self.mmc.mda.run(Sequence) ## We dont want this running in another thread based on my method


        print("Sequence Done")

    def GenerateSquareLoopReverseReturnDiagonalSequence(self, Position, StepSize):
        Positions = []
        
        Pos = Position

        Positions.append((Pos[0],Pos[1],Pos[2]))
        Positions.append((Pos[0]+StepSize,Pos[1],Pos[2]))
        Positions.append((Pos[0]+StepSize,Pos[1]+StepSize,Pos[2]))
        Positions.append((Pos[0],Pos[1]+StepSize,Pos[2]))
        Positions.append((Pos[0],Pos[1],Pos[2]))
        Positions.append((Pos[0],Pos[1]+StepSize,Pos[2]))
        Positions.append((Pos[0]+StepSize,Pos[1]+StepSize,Pos[2]))
        Positions.append((Pos[0]+StepSize,Pos[1],Pos[2]))
        Positions.append((Pos[0],Pos[1],Pos[2]))
        Positions.append((Pos[0]+StepSize,Pos[1]+StepSize,Pos[2]))
        Positions.append((Pos[0],Pos[1],Pos[2]))

        return Positions








        
    def CaptureMovie(self,SaveState,duration,Iterator,**kwargs):
        Time = GetFormattedHSAP()
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))

        ExposureTime = kwargs.get("Exposure",100)
        FrameSize = kwargs.get("FrameSize", [2304,2304])
        CropType = kwargs.get("CropType", "Centered")
        MovieName = kwargs.get("MovieName", "_" + Time + "_" + PositionString  + str(Iterator) + ".tiff")
        Settings = kwargs.get("SettingsFrom", "External")

        ROI = kwargs.get("ROI", "Default")

        if ROI != "Default":

            self.mmc.setROI(ROI[0],ROI[1],ROI[2],ROI[3])

        MovieName = self.ExtractedFolderPath + MovieName

        if Settings == "Internal":
            ExposureTime = self.ExposureTime
            FrameSize = self.FrameSize
            CropType = self.CropType

        Iterator = str(Iterator)
        ## Save State:
        ## 1) "Save"
        ## Saves the image to the listed directory
        ## 2) "Return"
        ## Returns the image and relevant metadata
        ## 3) SaveAndReturn
        ## Does both


        #PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        #ImageName = self.ExtractedFolderPath + "_" + Time + "_" + PositionString  + Iterator + ".tiff"

        Interval = ExposureTime
        FramesToCollect = int(1/ExposureTime * 1000 * duration)
        print(FramesToCollect)
        CameraDevice = self.mmc.getCameraDevice()
        self.mmc.setCameraDevice(CameraDevice)

        self.mmc.setExposure(CameraDevice,ExposureTime) ## IN MILLISECONDS!!
        #self.mmc.setCircularBufferMemoryFootprint(16384) ## Size of the circular buffer (Stores images from continuous acquisitions)


        start_time = time.time()
        self.mmc.startSequenceAcquisition(FramesToCollect, Interval, True)
        while self.mmc.isSequenceRunning(): ## https://sourceforge.net/p/micro-manager/mailman/micro-manager-general/thread/03A4DB6AE86963478C9EB15787545DA23E39DC1A%40BE-EXCH-01.andortech.net/
            pass
        #self.mmc.startContinuousSequenceAcquisition(1)
        #while time.time() - start_time < duration + Interval: ## Ensures that the movie is exactly the right length, but this can also be done concurrently
        #    pass
        delta = time.time() - start_time
        self.mmc.stopSequenceAcquisition()
        
        MovieSlices = list()

        while self.mmc.getRemainingImageCount() > 0: ## This needs to be done while acquisition happens or the "circular buffer" will over flow
            if len(MovieSlices) == FramesToCollect:
                break
            else:
                ImageData = self.mmc.popNextImage()
                ImageData = ImageData.astype(np.int32)
                MovieSlices.append(ImageData.astype(np.int32))
        
        AssembledMovie = np.stack(MovieSlices,axis=0)
        if np.shape(AssembledMovie)[0] != FramesToCollect:
            print(np.shape(AssembledMovie))
        else:
            AssembledMovie = AssembledMovie.reshape((FramesToCollect,ROI[2],ROI[3]))
        print("Movie Shape: %s" % ([np.shape(AssembledMovie)]))
        
        self.mmc.clearCircularBuffer() ## Remove any old frames that arent from the current movie
        
        self.mmc.snapImage() ## THIS RETURNS A FLOAT64 DATA TYPE ARRAY THAT IMAGE J CANNOT READ WHEN IT IS STACKED INTO A MOVIE!!
        ImageData = self.mmc.getTaggedImage()
        MetaData = ImageData[1] ## Collecting an image for the metadata




        #PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        #ImageName = self.ExtractedFolderPath + "_" + Time + "_" + PositionString  + Iterator + ".tiff"

        if CropType == "Origin":
            Image = Image[FrameSize[0],FrameSize[1]]
    
        if CropType == "Centered":
            Shape = np.shape(AssembledMovie)
            CenterX = int(1/2*Shape[1])
            CenterY = int(1/2*Shape[2])

            CropXMin = int(CenterX - 1/2 * FrameSize[0])
            CropXMax = int(CenterX + 1/2 * FrameSize[0])
            CropYMin = int(CenterY - 1/2 * FrameSize[1])
            CropYMax = int(CenterY + 1/2 * FrameSize[1])

            AssembledMovie = AssembledMovie[:,CropXMin:CropXMax,CropYMin:CropYMax]



        if SaveState == "Save":
            OME_Xml = dicttoxml(MetaData)
            #print(OME_Xml)
            tifffile.imwrite(MovieName, AssembledMovie.astype(np.float32))

        if SaveState == "Return":
            return AssembledMovie, MetaData
    
        if SaveState == "SaveAndReturn":
            OME_Xml = dicttoxml(MetaData)
            #print(OME_Xml)
            tifffile.imwrite(MovieName, AssembledMovie.astype(np.float32))
            return AssembledMovie, MetaData
        
        AverageFrameTime = delta / np.shape(AssembledMovie)[0]
        print("Total Time: %s, Expected Time: %s, Average Frame Time: %s, Average Frame Rate: %s" % (delta, duration, AverageFrameTime, 1 / AverageFrameTime))


    def FocusStep(self,Offset,Iterator,stepsize,steps,name,**kwargs): ## Find focus via maximum peak brightness
        Method = kwargs.get("Method","BestGaussianFit")
        
        ## Method based on: https://www.csl.cornell.edu/~cbatten/pdfs/batten-image-processing-sem-slides-scanning2001.pdf
        ## https://www.statology.org/autocorrelation-python/
        ## http://hahnlab.com/downloads/protocols/2006%20Methods%20Enzym-Shen%20414%20Chapter%2032-opt.pdf
        ## https://pmc.ncbi.nlm.nih.gov/articles/PMC5119614/        

        self.FocusStepDebug = kwargs.get("Debug","False")
        
        SettingsPeakPicking = (5, 500, 20) ## Minimum Distance (Pixels), ## Minimum Peak Height, ## Minimum Distance From the Edge
        SettingsFitting = (20, 0.80) ## Viewing Window for Peaks (Choose so that function --> 0), ## Threshold R2
        SettingsFit = (10,10,1.7,1.7,0.001,100) ## X, Y, Sigma X, Sigma Y, Theta, offset
        SettingsPeakGrouping = (2,3) ## Maximum Distance, ## Area for peak comparisons
        SettingsConversionFactors = (0,10) ## Photoelectrons per incoming photon, ## Framerate of data collection

        ImageStack = np.empty([0,2304,2304])
        Iterator = str(Iterator)

        PreviousPosition = self.GetStagePosition()
        self.OffsetStage(Offset) ## Moving the stage for focus the image
        FocusArray = np.empty([0,1]) #pd.DataFrame(columns = ['Z-Position','Variance'])
        VarianceArray = np.empty([0,1])
        StepsToCheck = steps ## Scanning the entire distance of the stage
        self.mmc.setZPosition(PreviousPosition[2]-(1/2)*StepsToCheck*stepsize) ## Check lower and higher than current focus setting
        Unfocused = True
        x = 0

        Time = GetFormattedHSAP()
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.ExtractedFolderPath + "/AutoFocus" + Iterator + "/" + name + ".csv"
        
        while Unfocused: 
            
            Position = self.GetStagePosition() ## Get the position

            ## Averaging 5 images
            Image1, MD = self.CaptureImage("Return", Settings = "Internal") ## Get an image
            Image = Image1

            ImageStack = np.vstack((ImageStack,Image[None])) ## [None] adds a third array to enable this kind of stacking
            Variance = np.var(Image)
            FocusMetric = np.array([Position[2]])
            Variance = np.array([Variance])
            FocusArray = np.vstack((FocusArray, FocusMetric))
            VarianceArray = np.vstack((VarianceArray, Variance))

            if x < StepsToCheck:
                MovePosition = [0,0,stepsize]
                self.OffsetStage(MovePosition)
                x += 1
                #print(x)
            
            else:
                if Method == "BestGaussianFit":
                    PeakPickingSettings = (100, 13) ## Minimum peak height, minimum distance between peaks
                    PeakGroupingSettings = (4,0) ## How close can peaks be across frames before they are considered the same peak
                    PeakFittingSettings = (13,0) ## PSF distance (will be fit over this 2*distance)

                    ImageStack = ImageStack.astype(np.int32)
                    #CreateIfNotExist(self.ExtractedFolderPath + "Movie%s/" % self.MovieNumber)
                    #tifffile.imwrite(self.ExtractedFolderPath + "Movie%s/Image.tif" % self.MovieNumber, ImageStack)

                    TotalPeakList = []
                    for Frames in ImageStack:
                        Frame, PeakList, GradientX = IP.GetPeaksByGradient(Frames, 0)
                        for Peaks in PeakList:
                            TotalPeakList.append(Peaks)

                    OrganizedPeaks = IP.OrganizePeaksThroughTime(TotalPeakList, PeakGroupingSettings)
                
                    OrganizedPeaks = OrganizedPeaks[0:5]

                    PSF = 5
                    xaxis = np.linspace(0,PSF,PSF)
                    yaxis = np.linspace(0,PSF,PSF)
                    xaxis, yaxis = np.meshgrid(xaxis, yaxis)
                    FinalizedData = []
                    for Peaks in OrganizedPeaks:
                        PeakData = []
                        MovieToFit = ImageStack[:][int(Peaks[2]-PSF):int(Peaks[2]+PSF),int(Peaks[3]-PSF):int(Peaks[3]+PSF)]
                        for Frames in MovieToFit:                    
                            try:
                                popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis, yaxis), Frames.ravel())
                                PeakData.append([*popt])
                            except:
                                pass
                        FinalizedData.append(PeakData)

                    FinalizedData = np.array(FinalizedData)
                    print(FinalizedData)
                
             
                
                    ZValuesArray = FocusArray
                    IteratorFocus = 0
                    while IteratorFocus < np.shape(FinalizedData)[0] / np.shape(FocusArray)[0] - 1:
                        ZValuesArray = np.vstack((ZValuesArray,FocusArray))
    
                        IteratorFocus += 1
    
                    FinalizedData = np.hstack((ZValuesArray,FinalizedData))
    
    
                    if self.FocusStepDebug == True:
                        IP.WriteCsv2D_Data(FinalizedData,self.ExtractedFolderPath + "/AutoFocus" + Iterator + "/",name + ".csv",["Z-Position","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
                    #IP.WriteCsv2D_Data(FinDataArrayStacked,self.DataFolderPath + "/AutoFocus" + Iterator + "/",name + ".csv",["Z-Position","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
                    
    
                    AverageFocus = np.empty((0))
                    #["Z-Height","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"]
                    DataArray = []
                    #print("Big Data: %s" % (np.shape(FinalizedData)[0]))
                    for DataPoints in FinalizedData: ## Discard poorly fitted data
                        Discard = False
                        DataToCheck = [DataPoints[0],1,DataPoints[3],DataPoints[4],(DataPoints[3]+DataPoints[4])/2] ## Z-Position, R2, Sigma X, Sigma Y, Average of Sigma X and Sigma Y
                        
                        if DataToCheck[1] < 0.85:
                            #print(DataToCheck)
                            Discard = True
                        
                        if DataToCheck[1] > 1:
                            Discard = True
    
                        if Discard == False:
                            DataArray.append(DataToCheck)
                            #print(DataToCheck)
                        #print(Discard,DataToCheck)
    
                    if self.FocusStepDebug == True:
                        IP.WriteCsv2D_Data(DataArray,self.ExtractedFolderPath + "/AutoFocus" + Iterator + "/",name + "CheckedFits.csv",["Z-Position","R2","Sigma X","Sigma Y","Average Sigma"])
    
                    PreviousData = [0,0,0,0,0]
                    LowestAveSigmaArray = []
                    LowestAveSigma = [0,100000]
                    
                    #print("Length of DataArray: %s" % len(DataArray))
                    for Data in DataArray:
                        if Data[0] < PreviousData[0]:
                            #print(LowestAveSigma)
                            LowestAveSigmaArray.append(LowestAveSigma)
                            LowestAveSigma = [0,100000]
                        
                        #print(Data[4]<LowestAveSigma[1])
                        if Data[4] < LowestAveSigma[1]:
                            LowestAveSigma = [Data[0],Data[4]]
                        PreviousData = Data
    
                    LowestAveSigmaArray.append(LowestAveSigma)
                    
                    AverageZPosition = 0
                    FociCount = 0
                    for Foci in LowestAveSigmaArray:
                        AverageZPosition = AverageZPosition + Foci[0]
                        FociCount += 1
                    
                    AverageZPosition = AverageZPosition / FociCount
                    print("Average Z Position: %s" % AverageZPosition)
    
                    FocusArrayVariance = np.hstack((FocusArray,VarianceArray))
                    
                    
                    MaxVariance = np.max(VarianceArray)
                    MaxFocusPos = np.where(MaxVariance == VarianceArray)
                    MaxFocus = FocusArray[MaxFocusPos]
                
                if Method == "HighestVariance":
                    FocusPos = np.where(np.max(VarianceArray) == VarianceArray)
                    AverageZPosition = FocusArray[FocusPos]


                #print(MaxFocus[0])

                #np.savetxt(ImageName,FocusArrayVariance,delimiter = ",")
                self.mmc.setZPosition(AverageZPosition)
                FOCUSSPOT = AverageZPosition
                #self.mmc.setZPosition(MaxFocus[0]) ## Reenable if this breaks
                #FOCUSSPOT = MaxFocus[0]
                Unfocused = False
        
        if self.FocusStepDebug == True:
            Time = GetFormattedHSAP()
            IntermediateImage, MD = self.CaptureImage("Return") ## Before Picture
            PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
            ImageName = self.ExtractedFolderPath + "/AutoFocus" + Iterator + "/" + name + ".tiff"
            tifffile.imwrite(ImageName, IntermediateImage)

        Offset = [-Offset[0],-Offset[1],0] ## Dont change Z once focused
        self.OffsetStage(Offset) ## Returning to the originial posisition

        
        return FOCUSSPOT, ImageStack


    def FocusInStagesByOffset(self,Offset,Iterator,**kwargs):
        self.FocusDebug = kwargs.get("Debug","False")
        MethodName = kwargs.get("Method","BestGaussianFit")

        Iterator = str(Iterator)
        if self.FocusDebug == True:
            Time = GetFormattedHSAP()
            BeforeImage, MD = self.CaptureImage("Return",Settings = "Internal") ## Before Picture
            PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
            ImageName = self.ExtractedFolderPath + "/AutoFocus" + Iterator + "/BeforeFocus.tiff"
            FolderName = self.ExtractedFolderPath + "/AutoFocus" + Iterator + '/'
            CreateIfNotExist(FolderName)

            tifffile.imwrite(ImageName, BeforeImage, compression = None)
        
        stepsize = 0.3 ## In micron, keep this relatively narrow, and hope your sample is flat otherwise inversion of contrast effects will ruin this function.
        steps = 19 ## In count
        FocusRough, ImageStackRough = self.FocusStep(Offset,Iterator,stepsize,steps,"Rough",Debug = self.FocusDebug, Method = MethodName)
        if self.FocusDebug == True:
            ImageName = self.ExtractedFolderPath + "/AutoFocus" + Iterator + "/RoughImageStack.tiff"
            tifffile.imwrite(ImageName, ImageStackRough, compression = None, imagej = True)        
            print(FocusRough)

        stepsize = 0.05 ## In micron
        steps = 24 ## In count
        FocusCoarse, ImageStackCoarse = self.FocusStep(Offset,Iterator,stepsize,steps,"Coarse",Debug = self.FocusDebug, Method = MethodName)
        if self.FocusDebug == True:
            ImageName = self.ExtractedFolderPath + "/AutoFocus" + Iterator + "/CoarseImageStack.tiff"
            tifffile.imwrite(ImageName, ImageStackCoarse, compression = None, imagej = True)
            print(FocusCoarse)
        
        self.FocusPosition = FocusCoarse
        return FocusCoarse

    def MoveByPositions(self, PositionArray):

        self.MoveStage(PositionArray[0])
        PositionArray.pop(0)

        return PositionArray

    def MoveInPattern(self, Function, Vars, RandomOffsetRange, **kwargs):

        AutofocusZ = kwargs.get("Autofocus", False)
        AutofocusMethod = kwargs.get("Method","BestGaussianFit")

        XVar = Vars[0]
        YVar = Vars[1]
        ZVar = Vars[2]




        Inx = self.MoveInPatternIndicies

        ## Vars should take the form of
        ## [[X_start,X_Stop,X_steps],[Y_start,Y_Stop,Y_steps],[Z_start,Z_Stop,Z_steps]]
        
        if Function == "set":
            RandomOffset = random.uniform(-1/2*RandomOffsetRange,1/2*RandomOffsetRange) ## This random offset should only be generated once
            self.XRange = np.linspace(XVar[0]+RandomOffset,XVar[1]+RandomOffset,XVar[2])
            self.YRange = np.linspace(YVar[0]+RandomOffset,YVar[1]+RandomOffset,YVar[2])
            self.ZRange = np.linspace(ZVar[0],ZVar[1],ZVar[2])

            print(self.GetStagePosition)

            if AutofocusZ == True:
                Offset = [0.5,0.5,0]
                FocusSpot, ImageStack = self.FocusInStagesByOffset(Offset,self.MovieNumber,Method = AutofocusMethod, **kwargs)
                return Incomplete, FocusSpot

            self.MoveInPatternLengths = [len(self.XRange), len(self.YRange), len(self.ZRange)]
            self.Iterator = 1

        if Function == "move":
            ## Check if overflowing

            if self.MoveInPatternIndicies[2] == self.MoveInPatternLengths[2]:
                self.MoveInPatternIndicies[2] = 0
                self.MoveInPatternIndicies[1] += 1

            if self.MoveInPatternIndicies[1] == self.MoveInPatternLengths[1]:
                self.MoveInPatternIndicies[1] = 0
                self.MoveInPatternIndicies[0] += 1

            if self.MoveInPatternIndicies[0] == self.MoveInPatternLengths[0]:
                Complete = "Complete"
                return Complete
            
            if np.all(self.MoveInPatternIndicies == self.MoveInPatternLengths): ## Checks if the movements are completed first
                Complete = "Complete"
                return Complete
            
            if self.MoveInPatternLengths == [1,1,1] and self.MoveInPatternIndicies[0] == 1:
                Complete = "Complete"
                return Complete
            
            if self.MoveInPatternLengths == [1,1,1] and self.MoveInPatternIndicies[1] == 1:
                Complete = "Complete"
                return Complete
            
            if self.MoveInPatternLengths == [1,1,1] and self.MoveInPatternIndicies[2] == 1:
                Complete = "Complete"
                return Complete

            if AutofocusZ == True:
                Offset = [0.5,0.5,0]
                FocusSpot, ImageStack = self.FocusInStagesByOffset(Offset,self.MovieNumber,Method = AutofocusMethod,**kwargs)
                self.MoveStageWithErrorCorrection([self.XRange[Inx[0]],self.YRange[Inx[1]],FocusSpot])
                self.Iterator += 1
                return Incomplete, FocusSpot

            else:
                self.MoveStageWithErrorCorrection([self.XRange[Inx[0]],self.YRange[Inx[1]],self.ZRange[Inx[2]]])

            self.MoveInPatternIndicies[2] += 1 ## Step the Z coordinate
            
        Incomplete = "Incomplete"
        return Incomplete















    def CaptureMDATimeSequence(self,ExposureTime, Duration, FPS, **kwargs):
        ROI = kwargs.get("ROI", "Default")

        if ROI != "Default":

            self.mmc.setROI(ROI[0],ROI[1],ROI[2],ROI[3])

        self.mmc.setExposure(ExposureTime)
        
        ## https://pymmcore-plus.github.io/pymmcore-plus/guides/mda_engine/
        self.FrameBuffer = []

        @self.mmc.mda.events.frameReady.connect
        def on_frame(image: np.ndarray, event: useq.MDAEvent):
            self.FrameBuffer.append(image)


        Sequence = MDASequence(
            config = {"config": "HammamatsuCam", "exposure": ExposureTime},
            time_plan = {"interval": 1 / FPS, "loops": Duration * FPS},
        )

        FramesToCollect = Duration * FPS

        start_time = time.time()
        
        self.mmc.run_mda(Sequence)
        
        
        while len(self.FrameBuffer) != FramesToCollect:
            pass

        stop_time = time.time()
        delta = stop_time - start_time

        print("Time: %s, Average FPS: %s, Desired FPS: %s" % (delta, delta / FramesToCollect, FPS))

        AssembledMovie = np.array(self.FrameBuffer)
        #print(np.shape(AssembledMovie))
        #while self.mmc.isSequenceRunning() == True:
        #    pass    
            
        #AssembledMovie = np.stack(self.FrameBuffer,axis=0)
        #AssembledMovie = self.FrameBuffer.reshape((FramesToCollect,ROI[2],ROI[3]))

        #self.FrameBuffer = []

        return AssembledMovie.astype(np.float32)





import numpy as np
import tifffile
import GeneralLibrary as GL
import scipy.optimize as opt
import fits
import scipy.ndimage as ndimage
import concurrent

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

def GetPeaksByGradient(FrameToPick, FrameNumber):
    Frame = np.zeros(np.shape(FrameToPick))
    Frame += FrameToPick
    
    ## Denoising the peak
    DenoiseSTD = 1 # 0.45
    Threshold = 7
    ViewWindow = 4
    
    FrameShape = np.shape(Frame)
    
    ## Setting the edges to a low value so peak picking doesnt occur here
    Frame[0:ViewWindow] = 20
    Frame[FrameShape[0]-ViewWindow:] = 20
    Frame[0:FrameShape[0],0:ViewWindow] = 20
    Frame[0:FrameShape[0],FrameShape[1]-ViewWindow:] = 20
    
    ## Denoising the input data by applying a gaussian filter to it
    FrameDenoised = np.zeros(np.shape(Frame)) + ndimage.gaussian_filter(Frame, DenoiseSTD)
    
    ## Calculating the gradient
    GradientXY = np.gradient(FrameDenoised)
    GradientX = GradientXY[0]
    GradientY = GradientXY[1]

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
        
                GradientX[AveCoords[0]-ViewWindow:AveCoords[0]+ViewWindow,AveCoords[1]-ViewWindow:AveCoords[1]+ViewWindow] = 0
                FrameMax = Frame[AveCoords[0],AveCoords[1]]
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
    Frame, PeakList, GradientX = GetPeaksByGradient(args[0], args[1])
    return Frame, PeakList, GradientX

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
            tifffile.imwrite(self.ExtractedDataPath+"PickedPeaks.tif" %s ,self.MoviePicked.astype(np.float32))
            print(len(self.PickedPeaks))

    def GetPeaksByGradient(self):


        CompletePeakList = []
        EditedPickedMovie = []
        EditedGradientX = []

        with concurrent.futures.ProcessPoolExecutor(max_workers=61) as executor: 
            ## This is the hard limit for workers on windows: 
            ## See: https://github.com/python/cpython/issues/71090
            ## See: https://github.com/psf/black/issues/564
            ArrayCounter = 0
            DataPool = []
            for Frames in self.MovieData:
                ProcessData = [Frames, ArrayCounter]
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



        #EditedPickedMovie = np.empty((0,self.PixelsX,self.PixelsY))
        #EditedGradientX = np.empty((0,self.PixelsX,self.PixelsY))
        #FrameCounter = 0
        #for Frames in self.MovieData:
        #    EditedPeak, PeakListFromFrame, GradientX = GetPeaksByGradient(Frames, FrameCounter)
        #    FrameCounter += 1
        #    for Peaks in PeakListFromFrame:
        #        CompletePeakList.append(Peaks)
        
        #    EditedPeak = np.expand_dims(EditedPeak, axis = 0)
        #    GradientX = np.expand_dims(GradientX, axis = 0)

        #    EditedPickedMovie = np.vstack((EditedPickedMovie, EditedPeak))
        #    EditedGradientX = np.vstack((EditedGradientX,GradientX))

        EditedPickedMovie = np.array(EditedPickedMovie)
        EditedGradientX = np.array(EditedGradientX)

        if self.DebugPeakPicking == True:
            tifffile.imwrite(self.ExtractedDataPath + "PickedPeaks.tif", EditedPickedMovie.astype(np.float32))   
            tifffile.imwrite(self.ExtractedDataPath + "PickedGradientX.tif", EditedGradientX.astype(np.float32))
            print(len(CompletePeakList))
        
        self.PickedPeaks = CompletePeakList


    def TruncatePeakList(self,Count):
        self.PickedPeaks = self.PickedPeaks[0:Count]
    
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
        
        xaxis = np.linspace(0, 2*PSF, 2*PSF)
        yaxis = np.linspace(0, 2*PSF, 2*PSF)
        xaxis, yaxis = np.meshgrid(xaxis, yaxis)  
        
        BackgroundLevels = np.average(self.MovieData[0,0:PSF,0:PSF])
        BackgroundSum = np.sum(self.MovieData[0,0:PSF,0:PSF])
        self.FittedPeaksOverTime = []
        PeakNumber = 0
        for PeakArrays in self.OrganizedPeaks: ## OrganizedPeaks = [[Peak1Frame1,Peak1Frame2,...],[Peak2Frame1,Peak2Frame2,...],...]
            print("Fitting Peak #%s out of %s" % (PeakNumber,len(self.OrganizedPeaks)-1))
            PeakFitParameters = PeakArrays[0] ## Just taking the first peak fit info (the center should be the same across frames)
            ## PeakFitParameters
            ## Output is [Frame #, Peak Maximum, Peak X, Peak Y]
            
            AreaToFit = self.MovieData[0:self.FrameCount,int(PeakFitParameters[2]-PSF):int(PeakFitParameters[2]+PSF),int(PeakFitParameters[3]-PSF):int(PeakFitParameters[3]+PSF)]
            PeakFittingDataOverTime = self.FitPeaksOverMovieFunction(AreaToFit, PeakFitParameters, BackgroundSum, PSF, BackgroundLevels, xaxis, yaxis)
            PeakFittingOverTimeAdjusted = []
            FitPercent = 0
            for Data in PeakFittingDataOverTime:
                Data[7] = Data[7] * 65 ## Eventually set this to use data from the calibration files
                Data[8] = Data[8] * 65
                PeakFittingOverTimeAdjusted.append(Data)
                if Data[11] > 0:
                    FitPercent += 1
                    
            FitPercent = np.round(FitPercent / self.FrameCount * 100, 2)
                
            Header = ["Frame","Height","Peak X","Peak Y","Amplitude","Peak X", "Peak Y", "Sigma X (nm)", "Sigma Y (nm)", "Theta (Radians)", "Offset", "R2", "AreaSum5x5", "Photoelectrons per second"]
            GL.WriteCsv2D_Data(PeakFittingOverTimeAdjusted,self.ExtractedDataPath,"Peak%s_%sPercentFitted.csv" % (PeakNumber,FitPercent),Header)
            
            PeakNumber += 1

            self.FittedPeaksOverTime.append(PeakFittingDataOverTime)

    def GetFittedPeaks(self):
        return self.FittedPeaksOverTime
    
    def FitPeaksOverMovieFunction(self,AreaToFit, PeakFitParameters, BackgroundSum, PSF, BackgroundLevels, xaxis, yaxis):
        ## This could potentially be fixed with this: https://scipy-cookbook.readthedocs.io/items/FittingData.html
        ## Again this is split as a separate function as it is a very easily parrallizable task

        #print(BackgroundSum)
        AreaSum = np.sum(AreaToFit[PeakFitParameters[0]])

        PeakFittingDataOverTime = []
        SumOfPSFAtInitialFit = np.sum(AreaToFit[PeakFitParameters[0]]) ## Getting a PSF sized sum of the data
        FrameCounter = 0
        for Frames in AreaToFit: ## For every frame of the movie:
            if np.sum(Frames) > (AreaSum * 0.65): ## If the photons collected are different enough from the background try to fit it
                InitialGuess = (PeakFitParameters[1], PSF, PSF, 2, 2, 0, BackgroundLevels)
                FrameDataRaveled = Frames.ravel() ## Raveling is the flattening of an array
                FittingBound = [[PeakFitParameters[1]*0.5,0,0,1,1,-0.5,0],[PeakFitParameters[1]*2,PSF*2,PSF*2,PSF,PSF,0.5,118]]
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
                        try:
                            popt, pcov = opt.curve_fit(fits.twoD_Gaussian, (xaxis, yaxis), FrameDataRaveled, maxfev = 1000, p0 = InitialGuess, bounds = FittingBound)
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



































































































class Microscope():
    def __init__(self,Paths,ProjectName,Debug):
        self.ProjectName = ProjectName
        self.ConfigurationPath = Paths[0]
        self.MicroscopeDataPath = Paths[1] # Where the folder containing your data folder will be
        self.MicroscopeDeviceAdapterPath = [Paths[2]]

        self.MicroscopeSetupDebug = Debug

        Date = GetFormattedYMD()
        self.FileNameCore = "%s_%s" % (self.ProjectName,Date)

        ## Microscope Setup (DO NOT TOUCH)
        self.mmc = pymmcore.CMMCorePlus()  # Instance micromanager core
        self.mmc.setDeviceAdapterSearchPaths(self.MicroscopeDeviceAdapterPath)


        if self.MicroscopeSetupDebug == True:
            print(self.mmc.getDeviceAdapterSearchPaths()) ## Where is it looking for device adapters, should be a list of chars

        self.mmc.loadSystemConfiguration(self.ConfigurationPath) # Loading the configuration file  

        if self.MicroscopeSetupDebug == True:
            print(self.mmc.getVersionInfo()) # gives: 'MMCore version 11.5.0'
            print(self.mmc.getAPIVersionInfo()) # gives: 'Device API version 73, Module API version 10'
            print(self.mmc.getLoadedDevices()) # Returns loaded devices

        self.DataFolderPath = self.MicroscopeDataPath + ProjectName + "/"
        self.DataFolderPathWFileNameCore = self.MicroscopeDataPath + ProjectName + "/" + self.FileNameCore

        CreateIfNotExist(self.DataFolderPath)
        self.mmc.setXYPosition(0,0)
        self.mmc.setZPosition(75) ## Starts at the middle so it can focus up or down
        self.Position = [0,0,100]
        #self.mmcCore = PM.Core()
        #self.mmStudio = PM.Studio()
        self.MoveInPatternIndicies = [0,0,0]

        CameraDevice = self.mmc.getCameraDevice()
        self.CameraDevice = CameraDevice

    def CaptureImage(self,SaveState,*args):
        ExposureTime = args.get("Exposure",100)
        ## Save State:
        ## 1) "Save"
        ## Saves the image to the listed directory
        ## 2) "Return"
        ## Returns the image and relevant metadata
        ## 3) SaveAndReturn
        ## Does both

        CameraDevice = self.mmc.getCameraDevice()
        self.CameraDevice = CameraDevice
        self.mmc.setCameraDevice(CameraDevice)
        self.mmc.setExposure(CameraDevice,ExposureTime) ## IN MILLISECONDS!!
        #self.mmc.setProperty(CameraDevice,"HamamatsuHam_DCAM-ScanMode",1) ## Fix this eventually
 

        self.mmc.snapImage()
        ImageData = self.mmc.getTaggedImage()

        Time = GetFormattedHSAP()
        
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPathWFileNameCore + "_" + Time + "_" + PositionString + "1.tiff"
    
        Image = ImageData[0]
        MetaData = ImageData[1]

        if SaveState == "Save":
            OME_Xml = dicttoxml(MetaData)
            tifffile.imwrite(ImageName, Image, description = OME_Xml, metadata = None)

        if SaveState == "Return":
            return Image, MetaData
    
        if SaveState == "SaveAndReturn":
            OME_Xml = dicttoxml(MetaData)
            tifffile.imwrite(ImageName, Image, description = OME_Xml, metadata = None)
            return Image, MetaData

    def CaptureMovie(self,SaveState,framerate,exposuretime,duration,Iterator):
        Iterator = str(Iterator)
        ## Save State:
        ## 1) "Save"
        ## Saves the image to the listed directory
        ## 2) "Return"
        ## Returns the image and relevant metadata
        ## 3) SaveAndReturn
        ## Does both
        AssembledMovie = np.empty([0,2304,2304])
        Time = GetFormattedHSAP()
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        MovieName = self.DataFolderPathWFileNameCore + "_" + Time + "_" + PositionString  + Iterator + ".tiff"
    
        self.mmc.snapImage() ## THIS RETURNS A FLOAT64 DATA TYPE ARRAY THAT IMAGE J CANNOT READ WHEN IT IS STACKED INTO A MOVIE!!
        ImageData = self.mmc.getTaggedImage()
        MetaData = ImageData[1] ## Collecting an image for the metadata

        MovieSlices = list()

        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPathWFileNameCore + "_" + Time + "_" + PositionString  + Iterator + ".tiff"

        Interval = 1/framerate
        FramesToCollect = int(framerate * duration)
        CameraDevice = self.mmc.getCameraDevice()
        self.mmc.setCameraDevice(CameraDevice)
        self.mmc.setExposure(CameraDevice,Interval*1000) ## IN MILLISECONDS!!
        self.mmc.setCircularBufferMemoryFootprint(16384) ## Size of the circular buffer (Stores images from continuous acquisitions)

        self.mmc.startContinuousSequenceAcquisition(0)
        start_time = time.time()
        while time.time() - start_time < duration + Interval: ## Ensures that the movie is exactly the right length, but this can also be done concurrently
            pass
        self.mmc.stopSequenceAcquisition()
        #print(self.mmc.getRemainingImageCount())
    
        while self.mmc.getRemainingImageCount() > 0: ## This needs to be done while acquisition happens or the "circular buffer" will over flow
            if len(MovieSlices) == FramesToCollect:
                break
            else:
                #print(self.mmc.getRemainingImageCount())
                #t1 = time.time()
                ImageData = self.mmc.popNextImage()
                #t2 = time.time()
                ImageData = ImageData.astype(np.int16)
                #t3 = time.time()
                #AssembledMovie = np.vstack((AssembledMovie,ImageData[None]))
                #t4 = time.time()
                #print("T1: %s, T2: %s, T3: %s" % (t2-t1,t3-t2,t4-t3))
                #time.sleep(5)
                MovieSlices.append(ImageData.astype(np.int16))
        
        AssembledMovie = np.stack(MovieSlices,axis=0)
        if np.shape(AssembledMovie)[0] != FramesToCollect:
            print(np.shape(AssembledMovie))
        else:
            AssembledMovie = AssembledMovie.reshape((FramesToCollect,2304,2304))
        print("Movie Shape: %s" % ([np.shape(AssembledMovie)]))
        #for slices in MovieSlices:
        #    AssembledMovie = np.vstack((AssembledMovie,slices[None])) ## This is more efficient in bulk, do not do this step 1 by 1 for each frame
        #print(self.mmc.getRemainingImageCount())

        #print("done",self.mmc.getRemainingImageCount())
        #print(type(AssembledMovie[0,0,0]))
        
        self.mmc.clearCircularBuffer() ## Remove any old frames that arent from the current movie
        
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPathWFileNameCore + "_" + Time + "_" + PositionString  + Iterator + ".tiff"

        #print(len(MovieSlices))
        #for Frames in MovieSlices:
        #    AssembledMovie = np.vstack((AssembledMovie,Frames[None]))
       
        #AssembledMovie = np.array(MovieSlices)
        #AssembledMovie = AssembledMovie.astype(np.int16)
        #print(np.shape(AssembledMovie))
        #print(type(AssembledMovie[0,0,0]))
        #print(np.shape(AssembledMovie))

        if SaveState == "Save":
            OME_Xml = dicttoxml(MetaData)
            tifffile.imwrite(ImageName, AssembledMovie, imagej = True)

        if SaveState == "Return":
            return AssembledMovie, MetaData
    
        if SaveState == "SaveAndReturn":
            OME_Xml = dicttoxml(MetaData)
            tifffile.imwrite(ImageName, AssembledMovie, imagej = True)
            return AssembledMovie, MetaData
    
    def TestReturnString():
        String = "String"
        return String

    def GetStagePosition(self):
        XPos = self.mmc.getXPosition()
        YPos = self.mmc.getYPosition()
        ZPos = self.mmc.getZPosition()
        self.Position = [XPos,YPos,ZPos]
        return self.Position

    def MoveStage(self,SetPosition):
        try:
            Dx = SetPosition
            self.mmc.setRelativeXYZPosition(Dx[0],Dx[1],Dx[2])
            self.GetStagePosition()

        except:
            print(SetPosition,self.GetStagePosition())

    def SetStage(self,SetPosition):
        self.mmc.setXYPosition(SetPosition[0],SetPosition[1])
        self.mmc.setZPosition(SetPosition[2])
        OldPosition = self.GetStagePosition()

    def CorrectPosition(self,SetPosition):
        OldPosition = self.GetStagePosition()
        
        Delta = [0,0,0]
        Delta[0] = SetPosition[0] - OldPosition[0]
        Delta[1] = SetPosition[1] - OldPosition[1]
        Delta[2] = SetPosition[2] - OldPosition[2]

        self.mmc.setXYPosition(SetPosition[0]+Delta[0],SetPosition[1]+Delta[1])
        self.mmc.setZPosition(SetPosition[2]+Delta[2])
        OldPosition = self.GetStagePosition()

    def FocusStep(self,Offset,Iterator,stepsize,steps,name,**kwargs): ## Find focus via maximum peak brightness
        self.FocusStepDebug = kwargs.get("Debug","False")
        
        SettingsPeakPicking = (5, 500, 20) ## Minimum Distance (Pixels), ## Minimum Peak Height, ## Minimum Distance From the Edge
        SettingsFitting = (20, 0.80) ## Viewing Window for Peaks (Choose so that function --> 0), ## Threshold R2
        SettingsFit = (10,10,1.7,1.7,0.001,100) ## X, Y, Sigma X, Sigma Y, Theta, offset
        SettingsPeakGrouping = (2,3) ## Maximum Distance, ## Area for peak comparisons
        SettingsConversionFactors = (0,10) ## Photoelectrons per incoming photon, ## Framerate of data collection

        ImageStack = np.empty([0,2304,2304])
        Iterator = str(Iterator)

        PreviousPosition = self.GetStagePosition()
        self.MoveStage(Offset) ## Moving the stage for focus the image
        FocusArray = np.empty([0,1]) #pd.DataFrame(columns = ['Z-Position','Variance'])
        VarianceArray = np.empty([0,1])
        StepsToCheck = steps ## Scanning the entire distance of the stage
        self.mmc.setZPosition(PreviousPosition[2]-(1/2)*StepsToCheck*stepsize) ## Check lower and higher than current focus setting
        Unfocused = True
        x = 0

        Time = GetFormattedHSAP()
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/" + name + ".csv"
        
        while Unfocused: 
            
            Position = self.GetStagePosition() ## Get the position

            ## Averaging 5 images
            Image1, MD = self.CaptureImage("Return") ## Get an image
            Image = Image1

            #print(type(Image))
            ImageStack = np.vstack((ImageStack,Image[None])) ## [None] adds a third array to enable this kind of stacking
            Variance = np.max(Image)
            #Variance = np.sum(scipy.fft.fft2(Image[864:1440,864:1440])[144:432,144:432])
            #print(np.shape(Variance))
            #scipy.fft.ifft2(scipy.signal.fftconvolve(Image[864:1440,864:1440],Image[864:1440,864:1440],"full")*Image[864:1440,864:1440])
            #Variance = np.max(scipy.signal.correlate2d(Image[864:1440,864:1440],Image[864:1440,864:1440],"full","fill"))
            #print("FFT")
            #Variance = sm.tsa.acf(Image, nlags = 10)
            #Variance = np.var(Image[864:1440,864:1440]) ## Get the normalized variance
            FocusMetric = np.array([Position[2]])
            Variance = np.array([Variance])
            FocusArray = np.vstack((FocusArray, FocusMetric))
            VarianceArray = np.vstack((VarianceArray, Variance))

            if x < StepsToCheck:
                MovePosition = [0,0,stepsize]
                self.MoveStage(MovePosition)
                x += 1
                #print(x)
            
            else:
                ImageStack = ImageStack.astype(np.int16)
                PeaksMovie = IP.Movie("FocusMovie", name + "FocusMovie",self.DataFolderPath + "/AutoFocus" + Iterator + "/",SettingsPeakPicking,SettingsFitting,SettingsFit,SettingsPeakGrouping,SettingsConversionFactors,ForceLoad = True, Movie = ImageStack )
                PeaksMovie.GetPeaks()
                PeaksMovie.SortPeakArrayByHeight()
                PeaksMovie.FilterPeaks()
                
                PeaksMovie.OrganizePeaksThroughTimeNew()

                PeaksMovie.TruncatePeakList()
                PeaksMovie.FitPeaksOverMovie()
                
                PeaksMovie.FilterPeakFits()
                PeaksMovie.RefitBadFits()
                
                
                #PeaksMovie.FinalizePeaks()

                #FinalizedData = PeaksMovie.GetFinalizedPeakArray()
                FinalizedData = PeaksMovie.GetFittedPeaks()

                Shape1, Shape2, Shape3, Shape4, Shape5, Shape6 = PeaksMovie.GetPeakArrayShapes()
                #print(Shape1, Shape2, Shape3, Shape4, Shape5, Shape6)
                #print(FinalizedData)
                #print(np.shape(FocusArray),np.shape(FinalizedData))
                ##FinDataArrayStacked = np.empty((0,20))
                #print(np.shape(FinalizedData)[0],np.shape(FinalizedData)[0]/50-1)
                ##BigData = np.zeros((13)) ## This grabs garbage data for whatever reason if you allow it to be np.empty
                ##IteratorFunction = 0
                #print(np.shape(FinalizedData)[0]/np.shape(ImageStack)[0]-1)
                #print(np.shape(ImageStack)[0]-1)
                #print("Array Shape")
                ##while IteratorFunction < np.shape(FinalizedData)[0]/np.shape(ImageStack)[0]:
                    ##IteratorTwo = 0
                    ##SmallData = np.empty((12))
                    #print(IteratorFunction)
                    ##while IteratorTwo < np.shape(ImageStack)[0]-1:
                        #print(50*IteratorFunction+IteratorTwo)
                        #print(FinalizedData[(np.shape(ImageStack)[0]-1)*IteratorFunction+IteratorTwo])
                        #print("Data")
                        #time.sleep(1/5)
                        
                        ##SmallData = np.vstack((SmallData,FinalizedData[(np.shape(ImageStack)[0])*IteratorFunction+IteratorTwo]))
                        ##IteratorTwo += 1
                        ##print((np.shape(ImageStack)[0])*IteratorFunction+IteratorTwo,np.shape(FinalizedData))

                        
                        
                    #print(IteratorFunction)
                    #print("Iterator Function")
                    #print(SmallData)
                    
                    ##SmallData = np.hstack((FocusArray,SmallData))
                    ##SmallData = np.vstack((SmallData,np.zeros(13)))
                    ##BigData = np.vstack((BigData,SmallData))
                    ##IteratorFunction += 1 
                #print((np.shape(ImageStack)[0])*IteratorFunction+IteratorTwo)
                #time.sleep(6000)
                ##FinDataArrayStacked = BigData
                #for Data in FinalizedData:
                    #Data = np.hstack((FocusArray,Data))    
                    #Data = np.vstack((Data,np.zeros(20)))
                    #FinDataArrayStacked = np.vstack((FinDataArrayStacked,Data))
                #FinalizedData = np.hstack((FocusArray,FinalizedData))
                #print(np.shape(FinDataArrayStacked))

                #ZValuesArray = np.vstack((FocusArray,FocusArray,FocusArray,FocusArray,FocusArray))
                ZValuesArray = FocusArray
                IteratorFocus = 0
                while IteratorFocus < np.shape(FinalizedData)[0] / np.shape(FocusArray)[0] - 1:
                    ZValuesArray = np.vstack((ZValuesArray,FocusArray))

                    IteratorFocus += 1
                FinalizedData = np.hstack((ZValuesArray,FinalizedData))


                if self.FocusStepDebug == True:
                    IP.WriteCsv2D_Data(FinalizedData,self.DataFolderPath + "/AutoFocus" + Iterator + "/",name + ".csv",["Z-Position","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
                #IP.WriteCsv2D_Data(FinDataArrayStacked,self.DataFolderPath + "/AutoFocus" + Iterator + "/",name + ".csv",["Z-Position","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"])
                

                AverageFocus = np.empty((0))
#["Z-Height","R2","Fit Height","Fit Peak X","Fit Peak Y","Sigma X","Sigma Y","Fit Theta","Fit Offset","Frame","Height","Peak X","Peak Y"]
                DataArray = []
                #print("Big Data: %s" % (np.shape(FinalizedData)[0]))
                for DataPoints in FinalizedData: ## Discard poorly fitted data
                    Discard = False
                    DataToCheck = [DataPoints[0],DataPoints[1],DataPoints[5],DataPoints[6],(DataPoints[5]+DataPoints[6])/2] ## Z-Position, R2, Sigma X, Sigma Y, Average of Sigma X and Sigma Y
                    
                    if DataToCheck[1] < 0.85:
                        #print(DataToCheck)
                        Discard = True
                    
                    if DataToCheck[1] > 1:
                        Discard = True

                    if Discard == False:
                        DataArray.append(DataToCheck)
                        #print(DataToCheck)
                    #print(Discard,DataToCheck)

                if self.FocusStepDebug == True:
                    IP.WriteCsv2D_Data(DataArray,self.DataFolderPath + "/AutoFocus" + Iterator + "/",name + "CheckedFits.csv",["Z-Position","R2","Sigma X","Sigma Y","Average Sigma"])

                PreviousData = [0,0,0,0,0]
                LowestAveSigmaArray = []
                LowestAveSigma = [0,100000]
                
                #print("Length of DataArray: %s" % len(DataArray))
                for Data in DataArray:
                    if Data[0] < PreviousData[0]:
                        #print(LowestAveSigma)
                        LowestAveSigmaArray.append(LowestAveSigma)
                        LowestAveSigma = [0,100000]
                    
                    #print(Data[4]<LowestAveSigma[1])
                    if Data[4] < LowestAveSigma[1]:
                        LowestAveSigma = [Data[0],Data[4]]
                        #print(Data[0],Data[1],Data[4])
                        #print("Yes")
                        #print(LowestAveSigma,Data[0],Data[4])
                    PreviousData = Data

                LowestAveSigmaArray.append(LowestAveSigma)
                
                #print(LowestAveSigmaArray)
                #print("Length of DataArray: %s" % len(DataArray))
                AverageZPosition = 0
                FociCount = 0
                for Foci in LowestAveSigmaArray:
                    AverageZPosition = AverageZPosition + Foci[0]
                    FociCount += 1
                
                AverageZPosition = AverageZPosition / FociCount
                print("Average Z Position: %s" % AverageZPosition)

                FocusArrayVariance = np.hstack((FocusArray,VarianceArray))
                
                
                MaxVariance = np.max(VarianceArray)
                MaxFocusPos = np.where(MaxVariance == VarianceArray)
                MaxFocus = FocusArray[MaxFocusPos]
                #print(MaxFocus[0])

                #np.savetxt(ImageName,FocusArrayVariance,delimiter = ",")
                self.mmc.setZPosition(AverageZPosition)
                FOCUSSPOT = AverageZPosition
                #self.mmc.setZPosition(MaxFocus[0]) ## Reenable if this breaks
                #FOCUSSPOT = MaxFocus[0]
                Unfocused = False
        
        if self.FocusStepDebug == True:
            Time = GetFormattedHSAP()
            IntermediateImage, MD = self.CaptureImage("Return") ## Before Picture
            PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
            ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/" + name + ".tiff"
            tifffile.imwrite(ImageName, IntermediateImage)

        Offset = [-Offset[0],-Offset[1],0] ## Dont change Z once focused
        self.MoveStage(Offset) ## Returning to the originial posisition

        
        return FOCUSSPOT, ImageStack

    def FocusInStagesByOffset(self,Offset,Iterator,**kwargs):
        self.FocusDebug = kwargs.get("Debug","False")

        Iterator = str(Iterator)
        if self.FocusDebug == True:
            Time = GetFormattedHSAP()
            BeforeImage, MD = self.CaptureImage("Return") ## Before Picture
            PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
            ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/BeforeFocus.tiff"
            FolderName = self.DataFolderPath + "/AutoFocus" + Iterator + '/'
            CreateIfNotExist(FolderName)

            tifffile.imwrite(ImageName, BeforeImage, compression = None)
        
        stepsize = 0.3 ## In micron, keep this relatively narrow, and hope your sample is flat otherwise inversion of contrast effects will ruin this function.
        steps = 19 ## In count
        FocusRough, ImageStackRough = self.FocusStep(Offset,Iterator,stepsize,steps,"Rough",Debug = self.FocusDebug)
        if self.FocusDebug == True:
            ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/RoughImageStack.tiff"
            tifffile.imwrite(ImageName, ImageStackRough, compression = None, imagej = True)        
            print(FocusRough)

        stepsize = 0.05 ## In micron
        steps = 24 ## In count
        FocusCoarse, ImageStackCoarse = self.FocusStep(Offset,Iterator,stepsize,steps,"Coarse",Debug = self.FocusDebug)
        if self.FocusDebug == True:
            ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/CoarseImageStack.tiff"
            tifffile.imwrite(ImageName, ImageStackCoarse, compression = None, imagej = True)
            print(FocusCoarse)
        
        return FocusCoarse
        #stepsize = 0.05 ## In micron
        #steps = 24 ## In count
        #FocusFine, ImageStackFine = self.FocusStep(Offset,Iterator,stepsize,steps,"Fine")
        #ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/FineImageStack.tiff"
        #tifffile.imwrite(ImageName, ImageStackFine, compression = None, imagej = True)
        #print(FocusFine)

    def MoveInPattern(self,Function,Vars,RandomOffsetRange):
        XVar = Vars[0]
        YVar = Vars[1]
        ZVar = Vars[2]
        Inx = self.MoveInPatternIndicies

        ## Vars should take the form of
        ## [[X_start,X_Stop,X_steps],[Y_start,Y_Stop,Y_steps],[Z_start,Z_Stop,Z_steps]]
        
        if Function == "set":
            RandomOffset = random.uniform(-1/2*RandomOffsetRange,1/2*RandomOffsetRange) ## This random offset should only be generated once
            self.XRange = np.linspace(XVar[0]+RandomOffset,XVar[1]+RandomOffset,XVar[2])
            self.YRange = np.linspace(YVar[0]+RandomOffset,YVar[1]+RandomOffset,YVar[2])
            self.ZRange = np.linspace(ZVar[0],ZVar[1],ZVar[2])
            self.MoveInPatternLengths = [len(self.XRange), len(self.YRange), len(self.ZRange)]

        if Function == "move":
            ## Check if overflowing

            if self.MoveInPatternIndicies[2] == self.MoveInPatternLengths[2]:
                self.MoveInPatternIndicies[2] = 0
                self.MoveInPatternIndicies[1] += 1
            if self.MoveInPatternIndicies[1] == self.MoveInPatternLengths[1]:
                self.MoveInPatternIndicies[1] = 0
                self.MoveInPatternIndicies[0] += 1
            if self.MoveInPatternIndicies[0] == self.MoveInPatternLengths[0]:
                Complete = "Complete"
                return Complete
            
            if np.all(self.MoveInPatternIndicies == self.MoveInPatternLengths): ## Checks if the movements are completed first
                Complete = "Complete"
                return Complete
            
            if self.MoveInPatternLengths == [1,1,1] and self.MoveInPatternIndicies[0] == 1:
                Complete = "Complete"
                return Complete
            
            if self.MoveInPatternLengths == [1,1,1] and self.MoveInPatternIndicies[1] == 1:
                Complete = "Complete"
                return Complete
            
            if self.MoveInPatternLengths == [1,1,1] and self.MoveInPatternIndicies[2] == 1:
                Complete = "Complete"
                return Complete

            self.SetStage([self.XRange[Inx[0]],self.YRange[Inx[1]],self.ZRange[Inx[2]]])
            #self.CorrectPosition([self.XRange[Inx[0]],self.YRange[Inx[1]],self.ZRange[Inx[2]]])

            self.MoveInPatternIndicies[2] += 1 ## Step the Z coordinate
            
        Incomplete = "Incomplete"
        return Incomplete

    def Unload(self):
        self.mmc.unloadAllDevices()

    def GetCameraDevice(self):
        return self.CameraDevice
    
    def FocusMicroscopeFromOffsetPosition(self,Offset,Iterator): 
        ## Method based on: https://www.csl.cornell.edu/~cbatten/pdfs/batten-image-processing-sem-slides-scanning2001.pdf
        ## https://www.statology.org/autocorrelation-python/
        ## http://hahnlab.com/downloads/protocols/2006%20Methods%20Enzym-Shen%20414%20Chapter%2032-opt.pdf
        ## https://pmc.ncbi.nlm.nih.gov/articles/PMC5119614/

        
        ## Split This function up into a Focus (Coarse,Fine,Finest) + Focus Function (Stepsize, Steps)
        Iterator = str(Iterator)
        Time = GetFormattedHSAP()
        BeforeImage, MD = self.CaptureImage("Return") ## Before Picture
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/Before.tiff"
        FolderName = self.DataFolderPath + "/AutoFocus" + Iterator + '/'
        CreateIfNotExist(FolderName)

        tifffile.imwrite(ImageName, BeforeImage)
        PreviousPosition = self.GetStagePosition()
        self.MoveStage(Offset) ## Moving the stage for focus the image
        FocusArray = np.empty([0,1]) #pd.DataFrame(columns = ['Z-Position','Variance'])
        VarianceArray = np.empty([0,1])
        stepsize = 4 ## Stepping 4 um at a time ## Scanning the entire array of values
        StepsToCheck = 49 ## Scanning the entire distance of the stage
        self.mmc.setZPosition(PreviousPosition[2]-(1/2)*StepsToCheck*stepsize) ## Check lower and higher than current focus setting
        Unfocused = True
        x = 0

        Time = GetFormattedHSAP()
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/CoarseFocus.csv"
        
        while Unfocused: ## Coarse Focus
            
            Position = self.GetStagePosition() ## Get the position
            Image, MD = self.CaptureImage("Return") ## Get an image
            Variance = np.var(Image) ## Get the variance
            FocusMetric = np.array([Position[2]])
            Variance = np.array([Variance])
            FocusArray = np.vstack((FocusArray, FocusMetric))
            VarianceArray = np.vstack((VarianceArray, Variance))
            #pd.concat([FocusArray,FocusMetric], ignore_index = True) ## Append Data to the array

            if x < StepsToCheck:
                MovePosition = [0,0,stepsize]
                self.MoveStage(MovePosition)
                x += 1
                #print(x)
            
            else:
                #print(np.shape(FocusArray))
                #FocusArray = np.delete(FocusArray,0)
                FocusArrayVariance = np.hstack((FocusArray,VarianceArray))
                np.savetxt(ImageName, FocusArrayVariance, delimiter = ',')
                #highestvalue = 0
                #for values in FocusArray:
                #    if np.all(values[1] > highestvalue):
                #        highestvalue = values
                #    else:
                #        pass

                #MaxVariance = highestvalue

                MaxVariance = np.max(VarianceArray)
                MaxFocusPos = np.where(MaxVariance == VarianceArray)

                #print(MaxVariance)
                #print(FocusArray[MaxFocusPos])
                #print(MaxVariance)
                #print(FocusArray[1])
                #for X in FocusArray:
                #    if np.all(X == MaxVariance):
                #        MaxFocus = X[0]

                #print(FocusArray[0],FocusArray[1])
                #MaxFocusPosition = np.where(MaxVariance == FocusArray[1])
                #MaxFocus = FocusArray[0,MaxFocusPosition]
                MaxFocus = FocusArray[MaxFocusPos]
                print(MaxFocus[0])
                self.mmc.setZPosition(MaxFocus[0])
                Unfocused = False

        Time = GetFormattedHSAP()
        IntermediateImage, MD = self.CaptureImage("Return") ## Before Picture
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/CoarseFocus.tiff"
        tifffile.imwrite(ImageName, IntermediateImage)
        
        ## Fine Focusing
        IntermediatePosition = self.GetStagePosition()
        FineFocusArray = np.empty([0,1]) #pd.DataFrame(columns = ['Z-Position','Variance'])
        FineVarianceArray = np.empty([0,1])
        stepsize = 0.4 ## Stepping 4 um at a time ## Scanning the entire array of values
        StepsToCheck = 9 ## Scanning the entire distance of the stage
        self.mmc.setZPosition(IntermediatePosition[2]-(1/2)*StepsToCheck*stepsize) ## Check lower and higher than current focus setting
        Unfocused = True
        x = 0

        Time = GetFormattedHSAP()
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/FineFocus.csv"
        
        while Unfocused: ## Fine Focus
            
            Position = self.GetStagePosition() ## Get the position
            Image, MD = self.CaptureImage("Return") ## Get an image
            Variance = np.var(Image) ## Get the variance
            FineFocusMetric = np.array([Position[2]])
            FineVariance = np.array([Variance])
            FineFocusArray = np.vstack((FocusArray, FocusMetric))
            FineVarianceArray = np.vstack((VarianceArray, FineVariance))
            #pd.concat([FocusArray,FocusMetric], ignore_index = True) ## Append Data to the array

            if x < StepsToCheck:
                MovePosition = [0,0,stepsize]
                self.MoveStage(MovePosition)
                x += 1
                #print(x)
            
            else:
                #print(np.shape(FineFocusArray))
                #FocusArray = np.delete(FocusArray,0)
                np.savetxt(ImageName, FineFocusArray, delimiter = ',')
                MaxVariance = np.max(FineVarianceArray)
                MaxFocusPos = np.where(MaxVariance == FineVarianceArray)
                #print(MaxVariance)
                #print(FocusArray[1])
                #for X in FineFocusArray:
                #    if X[1] == MaxVariance:
                #        MaxFocus = X[0]

                #print(MaxFocusPos)
                MaxFocus = FineFocusArray[MaxFocusPos]
                #print(FocusArray[0],FocusArray[1])
                #MaxFocusPosition = np.where(MaxVariance == FocusArray[1])
                #MaxFocus = FocusArray[0,MaxFocusPosition]
                print(MaxFocus[0])
                self.mmc.setZPosition(MaxFocus[0])
                Unfocused = False

        ## Fine Focusing End


        Time = GetFormattedHSAP()
        FineImage, MD = self.CaptureImage("Return") ## Before Picture
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/FineFocus.tiff"
        tifffile.imwrite(ImageName, FineImage)


        ## Finer Focusing
        IntermediatePosition = self.GetStagePosition()
        FinerFocusArray = np.empty([0,1]) #pd.DataFrame(columns = ['Z-Position','Variance'])
        FinerVarianceArray = np.empty([0,1])
        stepsize = 0.040 ## Stepping 5 nm at a time ## Scanning the entire array of values
        StepsToCheck = 9 ## Scanning the entire distance of the stage
        self.mmc.setZPosition(IntermediatePosition[2]-(1/2)*StepsToCheck*stepsize) ## Check lower and higher than current focus setting
        Unfocused = True
        x = 0

        Time = GetFormattedHSAP()
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/FinestFocus.csv"
        
        while Unfocused: ## Fine Focus
            
            Position = self.GetStagePosition() ## Get the position
            Image, MD = self.CaptureImage("Return") ## Get an image
            Variance = np.var(Image) ## Get the variance
            FocusMetric = np.array([Position[2]]) #pd.DataFrame({'Z-Position': [Position[2]],'Variance': [Variance]})
            VarianceArray = np.array([Variance])
            FinerFocusArray = np.vstack((FinerFocusArray, FocusMetric))
            FinerVarianceArray = np.vstack((FinerVarianceArray, VarianceArray))
            #pd.concat([FocusArray,FocusMetric], ignore_index = True) ## Append Data to the array

            if x < StepsToCheck:
                MovePosition = [0,0,stepsize]
                self.MoveStage(MovePosition)
                x += 1
                #print(x)
            
            else:
                #print(np.shape(FineFocusArray))
                #FocusArray = np.delete(FocusArray,0)
                np.savetxt(ImageName, FinerFocusArray, delimiter = ',')
                MaxVariance = np.max(FinerVarianceArray)
                MaxFocusPos = np.where(MaxVariance == FinerVarianceArray)
                #print(MaxVariance)
                #print(FocusArray[1])
                #for X in FineFocusArray:
                #    if X[1] == MaxVariance:
                #        MaxFocus = X[0]
                
                MaxFocus = FinerFocusArray[MaxFocusPos]

                #print(FocusArray[0],FocusArray[1])
                #MaxFocusPosition = np.where(MaxVariance == FocusArray[1])
                #MaxFocus = FocusArray[0,MaxFocusPosition]
                print(MaxFocus[0])
                self.mmc.setZPosition(MaxFocus[0])
                Unfocused = False

        ## Finer Focusing End

        Time = GetFormattedHSAP()
        FineImage, MD = self.CaptureImage("Return") ## Before Picture
        PositionString = "%s_%s_%s_" % (np.round(self.Position[0],3),np.round(self.Position[1],3),np.round(self.Position[2],3))
        ImageName = self.DataFolderPath + "/AutoFocus" + Iterator + "/FinestFocus.tiff"
        tifffile.imwrite(ImageName, FineImage)

        Offset = [-Offset[0],-Offset[1],0] ## Dont change Z once focused
        self.MoveStage(Offset) ## Returning to the originial posisition
        ReturnedPosition = self.GetStagePosition()
        DeltaArray = []
        Iterator = 0
        for Coords in PreviousPosition:
            DPos = PreviousPosition[Iterator] - ReturnedPosition[Iterator]            
            DeltaArray.append(DPos)
            Iterator += 1


        
        #FocusArray.to_csv(ImageName, index = False)

        return DeltaArray, FocusArray, FineFocusArray, FinerFocusArray



## Add a program for stitching images together
## Add a program for auto focusing
## Add a progam for automatically correcting offsets
## Add a program for collectin movies
## Add a program for automatically scanning in X,Y,Z Dimensions
## Add a program for stitching Z-frames into a single TIFF stack
## Just make the metadata go to a corresponding txt for simplicity
## Add functions for scanning across XY plane collecing Q movies
## Add functions for scanning through Z plane
## Create "Generate Save Position" Function
## Create separate folders for metadata, autofocusing information, etc.
## Learn String Formatting and create a way to pass a formatted dict into the class to et up automated scanning of movies
## Make focusing faster
## Add position correction to auto focus
## See if I can get moving to be more precise and accurate ## Maybe set discrete steps
## Add metadata as XML files
## Add a description/read me txt to each project folder
## Add a go to random position function
## Collect a grid of focus vs position on slide
## Add a gui (use HTML, supposidly it is easier)
## add a live FFT option
## Add a function that determines a "Quality score for the data based on the goodness of gaussian fits within the data"
## Add a printout to each function (i.e., on autofocusing, print (Starting Autofocus). print( Ending Autofocust: %s (Focus) %s (Time to focus)))
## Add a function to map some property (i.e. focus position) against position on the slide
## Add a tag function to the image and movie collection functions, for collecting and identifying multiple sequential movies of the same data
## Add XML metadata file for each capture
## Camera property "Scan Mode" controls the capture mode, ultraquiet, quiet, normal, etc...
## Add a function that finds the FWHM of beads from a movie over time
## Add autocropping and calibrating movies
## Add a wrapup function that closes out devices, saves debug output, finishes the program, etc.


## For integration with labview
def GetMicroscopeCameraDevice(Microscope):
    Camera = str(Microscope.GetCameraDevice())
    return Camera

def InitializeMicroscope(Path, ProjectName, Debug):
    MicroscopeObject = Microscope(Path,ProjectName,Debug)
    return MicroscopeObject

def KillMicroscope(Microscope):
    Microscope.Unload()
    
def SnapMicroscopeImage(MicroscopeObject,*args):
    Image, MetaData = MicroscopeObject.CaptureImage("Return",*args)
    return Image

def TestFunction(MicroscopeObject):
    ReturnString = MicroscopeObject.TestReturnString()
    return ReturnString

def ExecClassFunction(Class, Function, Parameters):
    ParametersList = ""
    for parameters in Parameters:
        ParametersList += str(parameters) + ","
    ParametersList = ParametersList[:-1]

    ExecutionString = "ValueToReturn = Class.%s(%s)" % (Function, ParametersList)
    Exec_Locals = {}    ## Values to extract from exec
    Exec_Globals = {"Class": Class} ## Values to pass into exec

    exec(ExecutionString, Exec_Globals, Exec_Locals)


def ExecClassFunctionReturn1(Class, Function, Parameters):
    ParametersList = ""
    for parameters in Parameters:
        ParametersList += str(parameters) + ","
    ParametersList = ParametersList[:-1]

    ExecutionString = "ValueToReturn = Class.%s(%s)" % (Function, ParametersList)
    Exec_Locals = {}    ## Values to extract from exec
    Exec_Globals = {"Class": Class} ## Values to pass into exec

    exec(ExecutionString, Exec_Globals, Exec_Locals)

    return Exec_Locals.get("ValueToReturn")

def ExecClassFunctionReturn2(Class, Function, Parameters):

    ParametersList = ""
    for parameters in Parameters:
        if type(parameters) != str:
            ParametersList += str(parameters) + ","
        else:
            ParametersList += parameters + ","
    ParametersList = ParametersList[:-1]

    ExecutionString = "ValueToReturn1, ValuesToReturn2 = Class.%s(%s)" % (Function, ParametersList)
    Exec_Locals = {}    ## Values to extract from exec
    Exec_Globals = {"Class": Class} ## Values to pass into exec
    
    #System.CaptureImage("Save")

    exec(ExecutionString, Exec_Globals, Exec_Locals)

    return [Exec_Locals.get("ValueToReturn1"), Exec_Locals.get("ValueToReturn2")]