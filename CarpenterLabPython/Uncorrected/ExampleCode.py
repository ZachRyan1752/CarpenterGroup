def InitializeMicroscope():
    import MicroscopeController as MC
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
        "C:/Program Files/Micro-Manager-2.0/"
    )
    ## 1) Configuration File
    ## 2) Where Data Folders Will Be Created
    ## 3) Micro Manager Directory

    System1 = MC.Microscope(Paths,"RonchyRullingTest",False) ##uManagerMovieCollectionTesting


    return System1

def StartMicroscopeProgram(System):
    Vars = [[35,65,1],[35,65,1],[75,100,1]] ## [X,Y,Z]
    Iterator = 1

    System.MoveInPattern("set",Vars,5)
    NumberOfPositions = Vars[0][2]*Vars[1][2]*Vars[2][2]


    while System.MoveInPattern("move",Vars,10) != "Complete":
        Position = System.GetStagePosition()
        
        print("X:{} um, Y:{} um, Z:{}".format(Position[0],Position[1],Position[2]))
        print("Auto-Focusing")
        System.FocusInStagesByOffset([0.5,0.5,0],Iterator) 
        print("Taking Movie")
        System.CaptureImage("Save")
        #System.CaptureMovie("Save",30,1/30,5,Iterator) ## 10 seconds --> 30 fps --> 300 frames, setting this higher than 30 fps changes cameramode automatically
        print("Movie Collected, Moving to Position: %s out of %s" % (Iterator,NumberOfPositions))
        Iterator += 1

    print("Program Done!")


def UnloadMicroscope(System):
    
    System.Unload()
    exit

## NOTE:
## THE SNAPIMAGE OPERATION PRODUCES A FLOAT-64 ARRAY THAT IMAGEJ CANNOT READ --> CONVERT TO FLOAT32 OR INT32

System = InitializeMicroscope()
StartMicroscopeProgram(System)