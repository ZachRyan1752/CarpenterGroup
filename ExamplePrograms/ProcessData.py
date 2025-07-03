from bs4 import BeautifulSoup
import tifffile
import GeneralLibrary as GL
import ImageProcessing as IP
import time
import numpy as np

### New Settings ###
PeakPickingSettings = (2, 20, 6) ## Gaussian Smooth radius, minimum grad hieght, minimum distance
PeakGroupingSettings = (8,0) ## How close can peaks be across frames before they are considered the same peak
PeakFittingSettings = (13,0.65,0.8) ## PSF distance (will be fit over this 2*distance), min R2, minimum fraction of time that is fitted, recommended to be ~ 1 for beads, lower thresholds for SM samples
### New Settings ###

### Settings ###
SettingsPeakPicking = (5, 200, 20) ## Minimum Distance (Pixels), ## Minimum Peak Height, ## Minimum Distance From the Edge
SettingsFitting = (4, 0.80) ## Viewing Window for Peaks (Choose so that function --> 0), ## Threshold R2
SettingsFit = (10,10,1.7,1.7,0.001,100) ## X, Y, Sigma X, Sigma Y, Theta, offset
SettingsPeakGrouping = (7,7) ## Maximum Distance, ## Area for peak comparisons
SettingsConversionFactors = (0,20) ## Photoelectrons per incoming photon, ## Framerate of data collection

### Main Code ###
def Main(MoviePath, MovieNumber):
    MovieFullPath = DataPath + MoviePath
    MovieOutputPath = MoviePath.replace(".tif", "/")
    #MovieObject = IP.Movie("Movie%s" % MovieNumber,ExtractedDataPath+"Movie%s" % MovieNumber,MovieFullPath,SettingsPeakPicking,SettingsFitting,SettingsFit,SettingsPeakGrouping,SettingsConversionFactors)

    MovieData = tifffile.imread(MovieFullPath)
    if MovieNumber == 0:
        MovieDebug = True

    else:
        MovieDebug = False

    MovieObject = IP.Movie29_05_2025(MovieData, np.empty((0,2304,2304)), ExtractedDataPath + MovieOutputPath, PeakPickingSettings, PeakGroupingSettings, PeakFittingSettings, DebugPicking = MovieDebug)

    t0 = time.time()
    MovieObject.RunGetPeaksByGradient()
    t1 = time.time()
    print("Done Picking Peaks: %s Seconds" % (t1-t0))
    
    t2 = time.time()
    MovieObject.RunOrganizePeaksThroughTime()
    MovieObject.TruncateOrganizedPeaksList(20)
    t3 = time.time()
    print("Done Filtering Peaks: %s Seconds" % (t3-t2))

    t4 = time.time()
    MovieObject.RunFitPeaksOverMovie()
    t5 = time.time()
    print("Done Fitting Peaks: %s Seconds" % (t5-t4))

    
    t14 = time.time()
    print("Time to Finish Processing: %s" % (t14-t0))
    


if __name__ == "__main__":
    t_start = time.time()

    MovieNumber = 0

    DataPath = GL.AskForDirectory(InputMessage = "Grabbing data from: ")
    ExtractedDataPath = GL.AskForDirectory(InputMessage = "Extracting data to: ")

    MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif")

    for MoviePath in MoviePaths:
        Main(MoviePath, MovieNumber)
        MovieNumber += 1


    t_end = time.time()
    print("Time to Process All Movies: %s" % (t_end-t_start))
    ### Main Code ###