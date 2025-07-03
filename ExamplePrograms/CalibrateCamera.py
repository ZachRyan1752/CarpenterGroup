#from carpenterlabpython 
import GeneratesCMOSCalibrationData as GCCD

## This may take a while to complete
## VERIFY THE RESULTS OF CALIBRATION BEFORE YOU CONTINUE TO PROCESS DATA!!

## Gain calibration
if __name__ == '__main__': ## If this is removed concurrent.futures.ProcessPool executor will execute this entire library instead of the code of interest
    ## Run the gain calibration from movies, pass along the # of illumination levels
    GCCD.RunGainCalibrationFromMovies(5)
    ## Calculate the spatial calibration of the system from LPM of a ronchi ruling and camera pixel size (in um)
    GCCD.SpatialCalibrationFromRuling(200, 6.5)
    ## Calculate the beam profile and peak power, providing a laser intensity at sample in mW
    GCCD.GetBeamProfile(0.435)
    ## Apply the calibrations, pass along which calibrations if you would like to chooose from them
    GCCD.ApplyCalibrationsToMovies()