# Change log
## Version 3 (2025_07_03)
### GeneratesCMOSCalibrationData.py (Updated)
  - Removed unused multiprocessing and gc libraries
  - Removed extraneous commented out code

### ImageProcessing.py (Added)
  - Allows the basic processing of SM images and movie data
  - Containes various functions, such as:
  - Finding peaks by frame maximum
  - Finding peaks via gradient
  - Organizing peaks that overlap through time
  - Fitting identified peaks throughout the movie
  - A class that allows for data to be processed, and integrated concurrent.futures parallelization

### pyproject.toml (Added)
  - Allows the libraries to be access natively by python by using the command:
  - pip -e install .

### __init__.py (Added)
  - Allows the importation of the libraries from subfolders in the same parent directory:
  - Read this for usage: https://www.geeksforgeeks.org/python/what-is-__init__-py-file-in-python/

## Version 2
### GeneratesCMOSCalibrationData.py (Updated)
  - Made the opening of file explorer for selection of data default behaviour for all functions

## Version 1
### GeneratesCMOSCalibrationData.py (Added)
  - Takes calibration data and processes it into a gain calibration frame, a beam profile, a beam profile fit, and spatial calibration
  - Enables the auto application of said calibrations to data

### GeneralLibrary.py (Added)
  - Adds various functions, such as:
  - Scanning for files in a path via a tag
  - Prompting the user with file explorer menu for either a directory or files
  - Creating a path if it does not exists
  - Getting the distance between two points
  - Getting the R2 of a data and fitted data
  - Making an .XML file with inputted data
  - Creating a 2D FFT
  - Writeing 2D data to a .csv file
  - Importing and stitching all .csv files from a directory
  - Formatting an XML string to remove extraneous characters
