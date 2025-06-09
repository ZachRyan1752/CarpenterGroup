## This library contains general functions that are useful for all libraries, but do not necessarily belong in any particular library


# Import Library
import datetime
import os
import numpy as np
import xml.etree.cElementTree as ET
import lxml.builder
import lxml.etree
from xml.dom import minidom
import cv2
import math
from bs4 import BeautifulSoup
import csv

# Defs
def GetFormattedYMDHSAP():
    ## Returns the date and time in string format of Year_Month_Day_Hour_Minute_AM/PM
    DateAndTime = datetime.datetime.now()
    DateAndTimeFormatted = datetime.datetime.strftime(DateAndTime,"20%y_%m_%d_%I_%M_%p")
    return DateAndTimeFormatted


def CreateIfNotExist(Path):
    ## Checks if a folder exists, if it does not, it creates it
    if os.path.exists(Path):
        pass
    else:
        os.mkdir(Path)


def ScanForFilesInPathByTag(Path, Tag):
    ## Returns the paths of files within a directory with a certain tag present
    TaggedFiles = []
    for Files in os.listdir(Path):
        if Tag in Files:
            TaggedFiles.append(Files)

    return TaggedFiles

def Distance(Coordinates1,Coordinates2):
    distance = math.sqrt((Coordinates2[0]-Coordinates1[0])**2+(Coordinates2[1]-Coordinates1[1])**2)
    return distance

def getR2(rawdata, fitteddata):
    residuals = rawdata - fitteddata
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((rawdata-np.mean(rawdata))**2)
    R2 = 1 - ( ss_res / ss_tot )
    return R2

def MakeXMLWithHeader(Header, SubHeaders, SubHeadersDataNames, Data, FilePath, FileName):
    ## Creates a XML of arbitrary data
    ## Taken from: https://www.geeksforgeeks.org/create-xml-documents-using-python/

    XMLRoot = minidom.Document()
    XML = XMLRoot.createElement(Header)
    XMLRoot.appendChild(XML)

    XMLIterator = 0
    for SubHeader in SubHeaders:
        productChild = XMLRoot.createElement(SubHeader)
        productChild.setAttribute(SubHeadersDataNames[XMLIterator], str(Data[XMLIterator]))
        XML.appendChild(productChild)
    
        XMLIterator += 1

    XMLStringToFile = XMLRoot.toprettyxml(indent = "\t")
    with open(FilePath+FileName, "w") as File:
        File.write(XMLStringToFile)

def MakeXMLWithHeaderbs4(Header, SubHeaders, Data, FilePath, FileName):
    ## Creates a XML of arbitrary data
    ## Taken from: https://www.geeksforgeeks.org/create-xml-documents-using-python/
    XMLData = BeautifulSoup(features = "xml")

    Root = XMLData.new_tag(Header)
    XMLData.append(Root)

    XMLIterator = 0
    for SubHeader in SubHeaders:
        XMLTag = XMLData.new_tag(SubHeader)
        XMLTag.string = str(Data[XMLIterator])
        Root.append(XMLTag)
        
        

        XMLIterator += 1

    with open(FilePath+FileName, "w", encoding = "utf-8") as File:
        File.write(XMLData.prettify())


def twoD_FFT(Image):
    FFT = cv2.dft(np.float32(Image), flags = cv2.DFT_COMPLEX_OUTPUT)
    FFTCentered = np.fft.fftshift(FFT)

    FFTMagnitude = 20*np.log(cv2.magnitude(FFTCentered[:,:,0],FFTCentered[:,:,1]))

    FFTNormalized = cv2.normalize(FFTMagnitude, None, 0, 255, cv2.NORM_MINMAX, cv2.CV_8UC1)
    return FFTNormalized

def WriteCsv2D_Data(Data,Path,Filename,Header):
    #if not filename in os.listdir(path):
    with open(Path+Filename,'w') as CsvFile:
        CsvFile.write("")
            
    CsvFile.close()
     
    with open(Path+Filename,'a') as CsvFile:
        for Rows in Header:
            CsvFile.write(str(Rows)+",")
        
        CsvFile.write("\n")   
            
        for Rows in Data:
            for Columns in Rows:
                CsvFile.write(str(Columns)+",")
            CsvFile.write("\n")        
    
    CsvFile.close()

def ImportAllCsvFromDirectory(DataFolder):
    ## Import all data from a directory contained in .csv files and return them as a singular array
    ## Ordered as [[Peak 1 Data], [Peak 2 Data], ...]
    ## Where [Peak 1 Data] is ordered as [Frame #, Peak Height, Peak X, Peak Y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset, R2, Area sum 5x5, Photoelectrons Per Second]
    CSVFiles = ScanForFilesInPathByTag(DataFolder, ".csv")
    DataArrayOutput = []
    for DataPath in CSVFiles:
        with open(DataFolder + DataPath, mode = 'r') as File:
            DataArray = csv.reader(File)
            for Lines in DataArray:
                LinesOut = [DataFolder, DataPath, *Lines]
                DataArrayOutput.append(LinesOut)

    return DataArrayOutput

def FormatXMLString(String):
    String = String.replace("\n  ","")
    String = String.replace("\n ","")

    return String