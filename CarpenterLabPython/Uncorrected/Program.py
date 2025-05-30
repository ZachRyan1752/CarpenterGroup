import MicroscopeController as MC
import ImageProcessing as IP
import numpy as np
import time
import random
## Hints
## 1) If you arent capturing as long of movies as you want check that the size of you movie isnt at or near 16 Gb in size
## ^ Larger or similarly sized movies will overwhelm the circular buffer <-- This will be fixed, eventually

## Dont collect for too long --> Better to split one big collection into multiple small collections for imaging purposes
## The largest a tifffile can be is ~ 4 Gb before it stops working

Paths = (
    "C:/Users/ryan.1752/Desktop/UManagerFiles/Configuration/AbleTable_TrainingScope_HamCap_NanoDrive_07_05_2025.cfg",
    "D:/Microscope Data/", 
    "C://Program Files//Micro-Manager-2.0//"
)
## 1) Configuration File
## 2) Where Data Folders Will Be Created
## 3) Micro Manager Directory

## Experiment Details
## Limits of stage are 0-200 um in x 0-180 um in Y for whatever reason
## Scanning from 5 - 195 microns in X
## Scanning from 5 - 175 microns in Y
## No Z - scan (except for autofocus)

System1 = MC.Microscope(Paths,"AutoFocusOverSlideTest40umSteps",False) ##uManagerMovieCollectionTesting
#Vars1 = [[5,185,10],[5,165,9],[75,100,1]] ## [X,Y,Z]
Vars1 = [[5,185,2],[5,165,2],[75,100,1]] ## [X,Y,Z] 
Iterator1 = 1
FocusArrayXValuesHeader = [*np.linspace(5,185,2)]
FocusArrayYValuesHeader = [*np.linspace(5,165,2)]
#time.sleep(5)
FocusArray1 = [] ## Just a 2d Array of Z-Values 
XStack1 = []

def Collection(Iterator,Vars,System,FocusArray,XStack):
    System.MoveInPattern("set",Vars,0)
    NumberOfPositions = Vars[0][2]*Vars[1][2]*Vars[2][2]

    YStep = 0
    while System.MoveInPattern("move",Vars,0) != "Complete":
        Position = System.GetStagePosition()
        

        print("X:{} um, Y:{} um, Z:{}".format(Position[0],Position[1],Position[2]))

        print("Auto-Focusing")
        FocusStart = time.time()
        Focus = System.FocusInStagesByOffset([0.5,0.5,0],Iterator)
        FocusEnd = time.time()
        print("Focus Time: %s" % (FocusEnd-FocusStart))

        XStack.append(Focus)

        print("Capturing Photo")
        
        System.CaptureImage("Save")
        
        #System.CaptureMovie("Save",30,1/30,5,Iterator) ## 10 seconds --> 30 fps --> 300 frames, setting this higher than 30 fps changes cameramode automatically
        
        print("Movie Collected, Moving to Position: %s out of %s" % (Iterator,NumberOfPositions))
        
        YStep += 1
        Iterator += 1

        if YStep == 2: ## Y-step 
            YStep = -1
            FocusArray.append([XStack])
            XStack = []

    return FocusArray

FocusArray1 = Collection(Iterator1,Vars1,System1,FocusArray1,XStack1)
FocusArrayXValuesHeader.insert(0,np.float64(0))
XLabels = FocusArrayXValuesHeader # Adding a corner/vertex datapoint
YLabels = np.array(FocusArrayYValuesHeader).T


XLabelsText = []
for X in XLabels:
    Text = str(X)
    XLabelsText.append(Text)

YIterator = 0
for Y in YLabels:
    FocusArray1[YIterator].insert(0,Y)

    YIterator += 1


IP.WriteCsv2D_Data(FocusArray1,"D:/Microscope Data/AutoFocusOverSlideTest40umSteps/","FociArray.csv",XLabelsText)
print("Program Done!")

exit

## NOTE:
## THE SNAPIMAGE OPERATION PRODUCES A FLOAT-64 ARRAY THAT IMAGEJ CANNOT READ --> CONVERT TO FLOAT32 OR INT32