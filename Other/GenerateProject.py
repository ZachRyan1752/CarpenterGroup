## This library is for the generation of project data in the form of an XML file, standardizing the structure of data collection for easier processing

# Import Library
import os
import time
import lxml.etree
import lxml.builder
import xml.etree.cElementTree as ET
import datetime
from bs4 import BeautifulSoup

# Defs
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
        File.write(XMLData.prettify()) # 

def CreateIfNotExist(Path):
    ## Checks if a folder exists, if it does not, it creates it
    if os.path.exists(Path):
        pass
    else:
        os.mkdir(Path)
        

def GetFormattedYMDHSAP():
    ## Returns the date and time in string format of Year_Month_Day_Hour_Minute_AM/PM
    DateAndTime = datetime.datetime.now()
    DateAndTimeFormatted = datetime.datetime.strftime(DateAndTime,"20%y_%m_%d_%I_%M_%p")
    return DateAndTimeFormatted


def GenerateNewProjectDirectory(Initials, ProjectName, DataDirectory):
    ## Generates a directory for storing project information, and creates an associated XML file containing essential project info
    ## Some basic syntax for generating XML files can be found here:
    ## https://stackoverflow.com/questions/3605680/creating-a-simple-xml-file-using-python

    MainFolder = DataDirectory + Initials + ProjectName + "/"

    CreateIfNotExist(MainFolder)

    CalibrationFolder = MainFolder + "Calibration/"
    GainCalibrationFolder = CalibrationFolder + "GainCalibration/"
    BeamProfileFolder = CalibrationFolder + "BeamProfile/"
    CalibrationDataFolder = CalibrationFolder + "CalibrationData/"
    ResolutionCalibrationFolder = CalibrationFolder + "Resolution/"
    BackgroundCalibrationFolder = CalibrationFolder + "Background/"

    CreateIfNotExist(CalibrationFolder)
    CreateIfNotExist(GainCalibrationFolder)
    CreateIfNotExist(BeamProfileFolder)
    CreateIfNotExist(CalibrationDataFolder)
    CreateIfNotExist(ResolutionCalibrationFolder)
    CreateIfNotExist(BackgroundCalibrationFolder)

    RawDataFolder = MainFolder + "RawData/"
    CorrectedCroppedFolder = MainFolder + "CorrectedCroppedData/"

    CreateIfNotExist(RawDataFolder)
    CreateIfNotExist(CorrectedCroppedFolder)

    ExtractedData = MainFolder + "ExtractedData/"

    CreateIfNotExist(ExtractedData)

    #LibrariesFolder = MainFolder + "Libs/"

    #CreateIfNotExist(LibrariesFolder)

    DateAndTime = GetFormattedYMDHSAP()

    print(MainFolder)
    Header = "Project_Information"
    SubHeaders = ["Project_Name","Initials","Date_and_Time_of_Creation", "Directory"]
    Data = [ProjectName, Initials, DateAndTime, MainFolder]
    MakeXMLWithHeaderbs4(Header, SubHeaders, Data, MainFolder, "ProjectInformation.xml")


Initials = "WMC"
ProjectName = "Z_Stacks"
DataDirectory = "D:\Microscope Data/"

GenerateNewProjectDirectory(Initials, ProjectName, DataDirectory)