import MicroscopeController as MC
import ImageProcessing as IP
import GeneralLibrary as GL
import numpy as np
import time
import random
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
        Status, Focus = System.MoveInPattern("move", Vars, 0, Autofocus = True)


        
        #FocusStart = time.time()
        #Focus = System.FocusInStagesByOffset([0.5,0.5,0],Iterator)
        #FocusEnd = time.time()
        #print("Focus Time: %s" % (FocusEnd-FocusStart))
        FocusArray.append(Focus)

        print("Capturing Photo")
        
        System.CaptureImage("Save")
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
    StepX = 20
    StartY = 10
    StopY = 130
    StepY = 20


    Vars1 = [[StartX,StopX,StepY],[StartY,StopY,StepY],[0,0,0]] ## [X,Y,Z] 
    Iterator1 = 1
    FocusArrayXValuesHeader = [*np.linspace(5,185,2)]
    FocusArrayYValuesHeader = [*np.linspace(5,165,2)]
    #time.sleep(5)
    FocusArray1 = [] ## Just a 2d Array of Z-Values 


    FocusArray1 = Collection(Iterator1,Vars1,System1,FocusArray1)

    FocusArray = np.array(FocusArray1).reshape((int((StopX-StartX)/StepX),int((StopY-StartY)/StepY)))

    print(FocusArray)
    print("Program Done!")

    exit
#GL.WriteCsv2D_Data(FocusArray,)
#FocusArrayXValuesHeader.insert(0,np.float64(0))
#XLabels = FocusArrayXValuesHeader # Adding a corner/vertex datapoint
#YLabels = np.array(FocusArrayYValuesHeader).T


#XLabelsText = []
#for X in XLabels:
#    Text = str(X)
#    XLabelsText.append(Text)

#YIterator = 0
#for Y in YLabels:
#    FocusArray1[YIterator].insert(0,Y)

#    YIterator += 1


#IP.WriteCsv2D_Data(FocusArray1,"D:/Microscope Data/AutoFocusOverSlideTest40umSteps/","FociArray.csv",XLabelsText)


## NOTE:
## THE SNAPIMAGE OPERATION PRODUCES A FLOAT-64 ARRAY THAT IMAGEJ CANNOT READ --> CONVERT TO FLOAT32 OR INT32