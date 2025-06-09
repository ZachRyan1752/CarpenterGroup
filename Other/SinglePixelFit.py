import GeneralLibrary as GL
import tifffile
import numpy as np

DataPath = "D:\Microscope Data\ZDRSMTrainingZDR_Rho110\Calibration\GainCalibration"
MoviePaths = GL.ScanForFilesInPathByTag(DataPath,".tif")

def ExtractRegionFromImage(Image,Region):
    ImageCropped = Image[Region[0]:Region[1],Region[2]:Region[3]]

    return ImageCropped

def ProcessMovie(MoviePath, Region):
    Movie = tifffile.imread(MoviePath)

    MovieCropped = []
    for Frames in Movie:
        CroppedImage = ExtractRegionFromImage(Frames, Region)
        MovieCropped.append(CroppedImage)

    MovieCropped = np.array(MovieCropped)

    return MovieCropped

MovieNumber = 0
Region = [0,576,0,576]

GL.CreateIfNotExist(DataPath + "/CroppedMovies/")

for Paths in MoviePaths:

    MovieCropped = ProcessMovie(Paths, Region)
    tifffile.imwrite(DataPath + "/CroppedMovies/Movie%s.tif" % MovieNumber, MovieCropped)
    print(MovieNumber)
    MovieNumber += 1