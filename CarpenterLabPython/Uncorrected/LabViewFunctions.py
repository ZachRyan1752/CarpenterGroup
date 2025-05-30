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
    "C://Program Files//Micro-Manager-2.0//"
)
## 1) Configuration File
## 2) Where Data Folders Will Be Created
## 3) Micro Manager Directory


System = MC.Microscope(Paths,"AutoFocusLengthvsHeighTest",False) ##uManagerMovieCollectionTesting
Vars = [[35,65,1],[35,65,1],[75,100,1]] ## [X,Y,Z]
Iterator = 2


Image, MetaData = MC.ExecClassFunctionReturn2(System,'CaptureImage',["'Save'"])