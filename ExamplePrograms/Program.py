import MicroscopeController as MC
import ImageProcessing as IP
import GeneralLibrary as GL
import numpy as np
import time
import random
import tifffile
## Hints
## 1) If you arent capturing as long of movies as you want check that the size of you movie isnt at or near 16 Gb in size
## ^ Larger or similarly sized movies will overwhelm the circular buffer <-- This will be fixed, eventually

## Dont collect for too long --> Better to split one big collection into multiple small collections for imaging purposes
## The largest a tifffile can be is ~ 4 Gb before it stops working




def Collection(Iterator,Vars,System,FocusArray):
    print("Here")
    System.MoveInPattern("set", Vars, 0, Autofocus = True)
    print("Here 2")
    NumberOfPositions = Vars[0][2]*Vars[1][2]*Vars[2][2]

    
    Status = "Incomplete"
    while Status != "Complete":
        print("Auto-Focusing")
        Position = System.GetStagePosition()

        print("X:{} um, Y:{} um, Z:{}".format(Position[0],Position[1],Position[2]))
        Status, Focus = System.MoveInPattern("move", Vars, 0, Autofocus = True, Method = "HighestVariance")


        
        #FocusStart = time.time()
        #Focus = System.FocusInStagesByOffset([0.5,0.5,0],Iterator)
        #FocusEnd = time.time()
        #print("Focus Time: %s" % (FocusEnd-FocusStart))
        FocusArray.append(Focus)

        print("Capturing Photo")
        
        System.CaptureImage("Save", Settings = "Internal")
        #System.CaptureMovie("Save",30,1/30,5,Iterator) ## 10 seconds --> 30 fps --> 300 frames, setting this higher than 30 fps changes cameramode automatically
        
        print("Movie Collected, Moving to Position: %s out of %s" % (Iterator,NumberOfPositions))
        
        Iterator += 1


    return FocusArray

if __name__ == "__main__":
    MicroscopeDeviceAdapterPath = ["C://Program Files//Micro-Manager-2.0//"]
    ExtractedFolderPath = "RawData/"
    ConfigurationPath = "C:/Users/ryan.1752/Desktop/UManagerFiles/Configuration/AbleTable_TrainingScope_HamCam_PIStage_06_05_2025.cfg"

    ## Experiment Details
    ## Limits of stage are 0-200 um in x 0-180 um in Y for whatever reason
    ## Scanning from 5 - 195 microns in X
    ## Scanning from 5 - 175 microns in Y
    ## No Z - scan (except for autofocus)

    #System1 = MC.Microscope(Paths,"AutoFocusOverSlideTest40umSteps",False) ##uManagerMovieCollectionTesting
    System1 = MC.ImagingSystem(MicroscopeDeviceAdapterPath, ExtractedFolderPath, ConfigurationPath, 0, CenterStageOnStart = True)
    #Vars1 = [[5,185,10],[5,165,9],[75,100,1]] ## [X,Y,Z]
    StartX = 10
    StopX = 130 
    StepX = 60
    StartY = 10
    StopY = 130
    StepY = 60



    System1.SetInternalCollectionSettings(ExposureTime = 250, FrameSize = [1152, 1152], CropType = "Centered")

    Vars1 = [[StartX,StopX,StepY],[StartY,StopY,StepY],[0,0,0]] ## [X,Y,Z] 
    Iterator1 = 1
    FocusArrayXValuesHeader = [*np.linspace(5,185,2)]
    FocusArrayYValuesHeader = [*np.linspace(5,165,2)]
    #time.sleep(5)
    FocusArray1 = [] ## Just a 2d Array of Z-Values 


    FocusArray1 = Collection(Iterator1,Vars1,System1,FocusArray1)

    FocusArrayXIterator = 0
    FocusArrayStitched = []
    StepsX = int((StopX-StartX)/StepX)
    StepsY = int((StopY-StartY)/StepY)
    while FocusArrayXIterator < StepsX:
        FocusArrayYIterator = 0
        FocusColumn = []

        while FocusArrayYIterator < StepsY:
            
            FocusColumn.append(FocusArray1[FocusArrayXIterator*StepsY+FocusArrayYIterator])
            FocusArrayYIterator += 1

        FocusArrayStitched.append(FocusColumn)
        FocusArrayXIterator += 1

    FocusArrayStitched = np.array(FocusArrayStitched).astype(np.float32)
    print(np.shape(FocusArrayStitched))
    tifffile.imwrite("ExtractedData/FocusArray.tif",FocusArrayStitched)

    #print(FocusArray)
    print("Program Done!")

    exit