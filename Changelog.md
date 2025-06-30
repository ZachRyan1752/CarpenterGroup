# To-do List
## Documentation
### GeneralLibrary.py
### GeneratesCMOSCalibrationData.py
### ImageProcessing.py
### MicroscopeController.py
### Fits.py
### Graphing.py

## Code clean up
### GeneralLibrary.py
### GeneratesCMOSCalibrationData.py
  - Remove extraneous, commented out code
  - Clean up some of the debug outputs
  - Remove the old hard-wired based directory search method with the GL based file explorer search as default behaviour for all functions
### ImageProcessing.py
### MicroscopeController.py
  - Clean up some debug outputs
  - Remove extransous, commented out code
### Fits.py
### Graphing.py

## Implementation
### General
  - Create a default directory folder, where the calibrations files, raw data files, etc. directories are stored, so they dont have to be manually inputted every time
### GeneralLibrary.py
  - Add a preserve root structure option to a function that scans for files or directories
    - Scan for files you want to process in a head folder, then recreate that subfile format in the data deposition directory
   
  - Add a format for file name function
    - Takes some input variable or string, then it removes characters not suitable for file names and returns the string
### GeneratesCMOSCalibrationData.py
  - Add automatic output of histograms like manufacturers provide for gain calibration code
  - Remove OPENCV2 from spatial calibration code
  - Extract frequencies directly from numpy, rather than from peak picking

### ImageProcessing.py
### MicroscopeController.py
  - Add an events based microscope controller system
  - Add autodetection of all connected devices, and subsequent autoconfiguration of the config file, for when you dont want to use EVERY device all at once
### Fits.py
### Graphing.py

# Version 2.2.0 (2025_06_24)
## GeneratesCMOSCalibrationData.py:
  - Added some further documentation of the code
  - Made changes to some variable names for clarity
  - Add the ability to select data through a file explorer window, rather than manually setting the locations or using default directories
  

## MicroscopeController.py:
  - Added USEQ.MDASequence support: [Read-Here](https://pymmcore-plus.github.io/pymmcore-plus/guides/mda_engine/)
  - Added a movie collection function, capable of handling frame rates of up to 600 fps at 144 x 144 movies, 150 fps at 1152 x 1152 movies
  - Added a joint USEQ.MDASequence/movie collection function, for MDA based image acquisition
  - Implemented low latency movement before collection, to measure movement jitter in stages

## Added StageRinging.py:
  - Example script
  - Takes 600 Hz videos for 2 seconds in a particular pattern and saves the data for later processing

# Version 2.1.0 (2025_06_11)
## GeneratesCMOSCalibrationData.py:
  - Implemented further optimizations, tested using a supercomputer
  - Removed extraneous code

  - ### GainCalibrationFromMovies:
    - Removed extraneous code
    - Moved away from np.vstack to creating a list of arrays, then stitching them into a numpy movie to help save memory
    - Implemented bug fixes
    - Removed max_workers from concurrent.futures.ProcessThreadPool, as the default should be sufficient, and *too* many workers can actually slow the system from overhead.
    - Fixed accidentally subtracting illumination frames from themselves, instead of the dark frame from itself (to remove the offset)
    - Added read noise frame generation for comparison of our calibration steps to manufacturer steps
    - Added a wrapper function to prevent recursive calling of parallel execution functions.
   
    
  - ### ApplyCalibrationsToMovies:
    - Now outputs a file that explicetly states where each movie comes from after the calibration is done.


## Added SinglePixelFit.py
  - Goal: Crop data by any arbitrary ammount to produce smaller or differently shaped videos


## Added PlotHistogram.py
  - Goal: Plot a histogram to compare our gain calibration data against Hamamatsu's


## General:
  - Removed broken glass


# Version 2.0.0 (2025_06_06)
## GeneralLibrary.py:
- ### Added a ImportAllCsvFromDirectory function
  - Takes a directory, returns a compound array marked with the folder, path, and data, and returns a big array
- ### Added a FormatXMLString function
  - Takes a string from the XML reader and removes the extraneous \n and spaces from it


## GeneratesCMOSCalibrationData.py:
  - Optimized the GainCalibration function to use less ram (hopefully)
  - Added multiprocessing to multiple functions


## General:
- Added example functions for how to use the code to the GitHub
- Added a SM data simulation file
- Added the GenerateProject.py program for generating file directories
- Removed broken glass

  
Markdown cheatsheets: 

https://github.com/im-luka/markdown-cheatsheet
https://github.com/im-luka/markdown-cheatsheet?tab=readme-ov-file#user-content-fn-1-3cd6da4860b03d7ee151fb7ff7186453
