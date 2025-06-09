#from carpenterlabpython 
import GeneratesCMOSCalibrationData as GCCD
from bs4 import BeautifulSoup

ProjectInformationFile = "ProjectInformation.xml"
with open(ProjectInformationFile) as PIFFile:
    PIF = PIFFile.read()
    ProjectInformationXML = BeautifulSoup(PIF, "xml")

CalibrationInformation = ProjectInformationXML.find_all('Directory')
for Information in CalibrationInformation:
    CalibrationFolderPath = Information.string

## Data collection for calibration:
## Generally:
## 1) Movies should not have any saturated pixels
## 2) Movies should be 100 frames, taken at 100 ms exposures, ultraquiet mode (or another mode if assigned)
## 3) Gain calibration movies should be at least 10,000 frames
## 3) Movies should have the same or greater ROI then the movies you intend to take
## 4) Ronchi rulings should be single image
## 5) All collections should be done with the lights off
## 6) Carpet sample should be taken at an OD that does not saturate the camera (probably 2.5)
## 7) PVAVA (blank) and SM samples should be taken at the same OD.
## Collect these movies:
## 1) Various wide field illumination and no illumination movies with no objective in (or the objective at a different positon)
## --> Move to calibration/gaincalibration (or just save them directly there)
## 2) Collect images of ruanchy rulings in one orientation.
## --> Move to calibration/resolutioncalibration
## 3) Carpet sample with laser (Taken day of experiment)
## 3a) Get the laser power at the sample
## --> Move to calibration/beamprofile
## 4) PVAVA blank (Taken day of experiment)
### --> Move to calibration/background

## To verify calibration is accurate:
## Gain Calibration
## Check R2, offsets, and gain coefficients to ensure they make sense
## Space calibration
## Check the pixels/nm resolution and verify it against the optical path
## Beam profile
## Check R2 of fit


## This may take a while to complete

## Gain calibration
if __name__ == '__main__': ## If this is removed concurrent.futures.ProcessPool executor will execute this entire library instead of the code of interest
    #GCCD.GainCalibrationFromMovies(CalibrationFolderPath, 3)
    ## Folder, lines per mm of ronchi ruling, camera pixel size
    #GCCD.SpatialCalibrationFromRuling(CalibrationFolderPath, 200, 6.5)
    ## Folder, beam power in mW
    #GCCD.GetBeamProfile(CalibrationFolderPath, 0.435)
    ## Averageing the blank
    #GCCD.BlankAverage(CalibrationFolderPath)
    ## Apply the calibrations
    GCCD.ApplyCalibrationsToMovies(CalibrationFolderPath)