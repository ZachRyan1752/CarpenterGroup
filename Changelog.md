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
    - Now outputs a file that explicetely states where each movie comes from after the calibration is done.
    

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
