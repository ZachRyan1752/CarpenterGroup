## Import Library
import pymmcore_plus as pymmcore
import tifffile
import os
from datetime import datetime
from dicttoxml import dicttoxml
import numpy as np
import useq
#import pycromanager as PM
import time
import pandas as pd
import statsmodels.api as sm
import scipy
import random
import ImageProcessing as IP

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
    def __init__(self):
        self.Position = self.GetStagePosition() ## Getting the position of the stage


    def GetStagePosition(self): ## Saves the new position internally, then outputs a list containing the X, Y, and Z positions of the stage
        XPos = self.mmc.getXPosition()
        YPos = self.mmc.getYPosition()
        ZPos = self.mmc.getZPosition()
        self.Position = [XPos,YPos,ZPos]
        return self.Position

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