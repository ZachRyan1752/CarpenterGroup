import numpy as np
import tifffile
import random

def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset): # Taken from: https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
    x, y = xy
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

ImageDimensions = [576,576]
FrameCount = 100
BlinkingProbability = 0.185 ## Lower means more likely active, less means more likely blanked
QuenchingProbability = 0.065 ## Higher means more likely to quench

SMPeakCount = 50 ## Get more than you think you will need

SMPeakParameters = [800, 0, 0, 1.67, 1.67, 0, 0] ## 11000, 0, 0, 1.67, 1.67, 0, 0
BackgroundMax = 200
BackgroundSigmaX = 85
BackgroundSigmaY = 85
BackgroundTheta = 0
BackgroundOffset = 0

BeamProfileHieght = 19189
BPSigmaX = 55.5
BPSigmaY = 55.5
BPTheta = 0
BPOffset = 128

xaxis = np.linspace(0, ImageDimensions[0], ImageDimensions[0])
yaxis = np.linspace(0, ImageDimensions[1], ImageDimensions[1])
xaxis, yaxis = np.meshgrid(xaxis, yaxis)

Background = twoD_Gaussian((xaxis, yaxis), BackgroundMax, ImageDimensions[0] / 2, ImageDimensions[1] / 2, BackgroundSigmaX, BackgroundSigmaY, BackgroundTheta, BackgroundOffset).reshape((ImageDimensions[0],ImageDimensions[1]))
BeamProfile = twoD_Gaussian((xaxis, yaxis), BeamProfileHieght, ImageDimensions[0] / 2, ImageDimensions[1] / 2, BPSigmaX, BPSigmaY, BPTheta, BPOffset).reshape((ImageDimensions[0],ImageDimensions[1]))
BeamProfileNormalized = (BeamProfile - np.min(BeamProfile)) / (np.max(BeamProfile) - np.min(BeamProfile))

Peaks = np.zeros([FrameCount, *ImageDimensions])
while SMPeakCount > 0:
    FrameCountGenerator = 0
    FrameCountGenerator += FrameCount
    RandomX = random.randint(ImageDimensions[0] / 2 - BPSigmaX * 2, ImageDimensions[0] / 2 + BPSigmaX * 2)
    RandomY = random.randint(ImageDimensions[1] / 2 - BPSigmaY * 2, ImageDimensions[1] / 2 + BPSigmaY * 2)
    
    IntensityFactor = BeamProfileNormalized[RandomX,RandomY]
    
    PeakMax = SMPeakParameters[0] * IntensityFactor
    SigmaX = SMPeakParameters[3] * random.uniform(0.9,1.1)
    SigmaY = SMPeakParameters[3] * random.uniform(0.9,1.1)
    
    Theta = random.uniform(-0.2,0.2)
    
    Peak = twoD_Gaussian((xaxis, yaxis), PeakMax, RandomX, RandomY, SigmaX, SigmaY, Theta, 0).reshape((ImageDimensions[0],ImageDimensions[1]))
    Peak = np.expand_dims(Peak, axis=0)

    PeakOverTime = np.empty([0,*ImageDimensions])
    while FrameCountGenerator > 0:
        FramePeak = 0
        Probability = random.uniform(0,1)
        if Probability > BlinkingProbability:
            FramePeak += Peak * 1## Blinking
        if Probability < BlinkingProbability:
            FramePeak += Peak * 0
        if Probability > QuenchingProbability:
            pass
        if Probability < QuenchingProbability:
            Peak = Peak * 0
        
        PeakOverTime = np.vstack((PeakOverTime,FramePeak))
        
        FrameCountGenerator -= 1
    
    Peaks += PeakOverTime
    print("%s Peaks Remaining" % SMPeakCount)
    SMPeakCount -= 1

Noise = np.empty([0, *ImageDimensions])

FrameCountGenerator = 0
FrameCountGenerator += FrameCount
while FrameCountGenerator > 0:


    FrameNoise = np.zeros(ImageDimensions) + np.random.normal(115,12,ImageDimensions)
    FrameNoise = np.expand_dims(FrameNoise, axis=0)
    
    Noise = np.vstack((Noise,FrameNoise))

    FrameCountGenerator -= 1
    
    



SimulatedImage = np.zeros([FrameCount, *ImageDimensions]) + Peaks + Background #+ Noise


tifffile.imwrite("SimulatedSMMovie2.tif",SimulatedImage.astype(np.int16))