from bs4 import BeautifulSoup
import tifffile
import GeneralLibrary as GL
import ImageProcessing as IP
import time

### New Settings ###
PeakPickingSettings = (100, 13) ## Minimum peak height, minimum distance between peaks
PeakGroupingSettings = (4,0) ## How close can peaks be across frames before they are considered the same peak
PeakFittingSettings = (13,0) ## PSF distance (will be fit over this 2*distance)
### New Settings ###

### Settings ###
SettingsPeakPicking = (5, 200, 20) ## Minimum Distance (Pixels), ## Minimum Peak Height, ## Minimum Distance From the Edge
SettingsFitting = (4, 0.80) ## Viewing Window for Peaks (Choose so that function --> 0), ## Threshold R2
SettingsFit = (10,10,1.7,1.7,0.001,100) ## X, Y, Sigma X, Sigma Y, Theta, offset
SettingsPeakGrouping = (7,7) ## Maximum Distance, ## Area for peak comparisons
SettingsConversionFactors = (0,10) ## Photoelectrons per incoming photon, ## Framerate of data collection
### Settings ###

ProjectInformationFile = "ProjectInformation.xml"
with open(ProjectInformationFile) as PIFFile:
    PIF = PIFFile.read()
    ProjectInformationXML = BeautifulSoup(PIF, "xml")

ProjectFolderData = ProjectInformationXML.find_all('Directory')
for Information in ProjectFolderData:
    ProjectFolderPath = str(Information.text)

ProjectFolderPath = GL.FormatXMLString(ProjectFolderPath)

DataPath = ProjectFolderPath + "CorrectedCroppedData/"
ExtractedDataPath = ProjectFolderPath + "ExtractedData/"
MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif")
BlankFrameFolder = ProjectFolderPath + "Calibration/CalibrationData/BlankAverageGainCorrected.tif"
BlankData = tifffile.imread(BlankFrameFolder)
MovieNumber = 0

### Main Code ###
def Main(MoviePath):
    MovieFullPath = DataPath + MoviePath
    
    #MovieObject = IP.Movie("Movie%s" % MovieNumber,ExtractedDataPath+"Movie%s" % MovieNumber,MovieFullPath,SettingsPeakPicking,SettingsFitting,SettingsFit,SettingsPeakGrouping,SettingsConversionFactors)

    MovieData = tifffile.imread(MovieFullPath)
    MovieObject = IP.Movie29_05_2025(MovieData, BlankData, ExtractedDataPath+"Movie%s/" % MovieNumber, PeakPickingSettings, PeakGroupingSettings, PeakFittingSettings)

    t0 = time.time()
    MovieObject.GetPeaksByGradient()
    t1 = time.time()
    print("Done Picking Peaks: %s Seconds" % (t1-t0))
    
    t2 = time.time()
    #MovieObject.OrganizePeaksThroughTime()
    t3 = time.time()
    print("Done Filtering Peaks: %s Seconds" % (t3-t2))

    t4 = time.time()
    #MovieObject.FitPeaksOverMovie()
    t5 = time.time()
    print("Done Fitting Peaks: %s Seconds" % (t5-t4))

    
    t14 = time.time()
    print("Time to Finish Processing: %s" % (t14-t0))
    


if __name__ == "__main__":
    t_start = time.time()


    for MoviePath in MoviePaths:
        Main(MoviePath)


    t_end = time.time()
    print("Time to Process All Movies: %s" % (t_end-t_start))
    ### Main Code ###