import GeneralLibrary as GL
import time
import pandas as pd
import Graphing
import numpy as np

Folder = "D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData/Deviation/"
CSVPaths = GL.ScanForFilesInPathByTag(Folder,".csv")

SubDirectories = ["0um","0_001um","0_01um","0_1um","1um","10um"]
NameArray = ["Home","PosX1","PosY1","NegX1","NegY1","PosY2","PosX2","NegY2","NegX2","PosXY","NegXY"]
GroupedArray = [[],[],[],[],[],[],[],[],[],[],[]]



for SubDirs in SubDirectories[::-1]:
    SubDirMovies = GL.ScanForFilesInPathByTag(Folder,SubDirs)
    for Names in NameArray:
        Movies = GL.ScanForFilesInPathByTag(Folder,Names)
        DoubleCheckedPaths = []
        for MoviePath in Movies:
            if MoviePath in SubDirMovies and MoviePath == SubDirs + "Movie" + Names + ".csv":
                AppendIndex = NameArray.index(Names)
                GroupedArray[AppendIndex].append(MoviePath)
                #DoubleCheckedPaths.append(MoviePath)
            
DoubleGroupedArray = [[],[],[],[],[],[],[]]

Index2 = 0
PosX = []
PosY = []
NegX = []
NegY = []

for Items in GroupedArray:
    if Index2 == 1:
        for Item in Items:
            PosX.append(Item)
    if Index2 == 6:
        for Item in Items:
            PosX.append(Item)    
        for Item in PosX:
            DoubleGroupedArray[1].append(Item)

    if Index2 == 2:
        for Item in Items:
            PosY.append(Item)    
    if Index2 == 5:
        for Item in Items:
            PosY.append(Item)   
        for Item in PosY:
            DoubleGroupedArray[2].append(Item) 

    if Index2 == 3:
        for Item in Items:
            NegX.append(Item)   
    if Index2 == 8:
        for Item in Items:
            NegX.append(Item)   
        for Item in NegX:
            DoubleGroupedArray[3].append(Item)

    if Index2 == 4:
        for Item in Items:
            NegY.append(Item) 
    if Index2 == 7:
        for Item in Items:
            NegY.append(Item) 
        for Item in NegY:
            DoubleGroupedArray[4].append(Item)    

    if Index2 == 0:
        for Item in Items:
            DoubleGroupedArray[0].append(Item)

    if Index2 == 9:
        for Item in Items:
            DoubleGroupedArray[5].append(Item)

    if Index2 == 10:    
        for Item in Items:
            DoubleGroupedArray[6].append(Item)
    Index2 += 1


window_size = 1
FrameLimit = 150
window = np.ones(window_size) / window_size

SkipListArray = []
for items in NameArray:
    string = "0umMovie%s.csv" % items
    string2 = "0_001umMovie%s.csv" % items
    SkipListArray.append(string)
    #SkipListArray.append(string2)

NameArray2 = ["Home","PosX","PosY","NegX","NegY","PosXY","NegXY"]
BigIterator = 0
for Items in DoubleGroupedArray:
    XDataArray = [[],[]]
    YDataArray = [[],[]]
    Labels = [[],[]]
    SkipList = []
    Title = "Dynamic Stability"
    AxisLabels = ("Time (s)", "Deviation from Average (nm)")
    print(Items)
    print(Items[::-1])
    Iterator = 0
    Mult = 0
    for SubItem1 in Items:
        if BigIterator != -1:
            for items234 in SkipListArray:
                try:
                    Items.remove(items234)
                except:
                    pass
        for SubItem2 in Items[::-1]:
            if SubItem1 != SubItem2 and SubItem1[0:7] == SubItem2[0:7]:
                print("Yup: %s, %s" % (SubItem1, SubItem2))
                print(SubItem1, SubItem2)
                CSVData1 = pd.read_csv(Folder+SubItem1)
                CSVData2 = pd.read_csv(Folder+SubItem2)
        
                XAxis = CSVData1['Frames'] * 1 / 600
                YAxis = CSVData1['Frames'] * 1 / 600
                
                XAxis = XAxis[0:FrameLimit]
                YAxis = YAxis[0:FrameLimit]

                XData1 = CSVData1['X_DEVIATION']
                YData1 = CSVData1['Y_DEVIATION']

                XData2 = CSVData2['X_DEVIATION']
                YData2 = CSVData2['Y_DEVIATION']

                XData = np.average((XData1, XData2), axis = 0) + Iterator * Mult
                YData = np.average((YData1, YData2), axis = 0) + Iterator * Mult

                XData = np.convolve(XData, window, mode = "same")
                YData = np.convolve(YData, window, mode = "same")

                XData = XData[0:FrameLimit]
                YData = YData[0:FrameLimit]

                XDataArray[0].append(XAxis)
                XDataArray[1].append(YAxis)

                YDataArray[0].append(XData)
                YDataArray[1].append(YData)

                Labels[0].append("X-Axis Response: %s" % SubItem1.replace(".csv",""))
                Labels[1].append("Y-Axis Response: %s" % SubItem1.replace(".csv",""))
                SkipList.append(SubItem1)
                SkipList.append(SubItem2)
                Items.remove(SubItem2)

                Iterator += 1


        if SubItem1 not in SkipList:
            print("Here2")
            print(SubItem1, SubItem2, SkipList)
            CSVData = pd.read_csv(Folder+SubItem1)

            XAxis = CSVData['Frames'] * 1 / 600
            YAxis = CSVData['Frames'] * 1 / 600

            XAxis = XAxis[0:FrameLimit]
            YAxis = YAxis[0:FrameLimit]

            XData = CSVData['X_DEVIATION'] + Iterator * Mult
            YData = CSVData['Y_DEVIATION'] + Iterator * Mult

            XData = XData[0:FrameLimit]
            YData = YData[0:FrameLimit]
            
            XDataArray[0].append(XAxis)
            XDataArray[1].append(YAxis)

            XData = np.convolve(XData, window, mode = "same")
            YData = np.convolve(YData, window, mode = "same")

            YDataArray[0].append(XData)
            YDataArray[1].append(YData)

            Labels[0].append("X-Axis Response: %s" % SubItem1.replace(".csv",""))
            Labels[1].append("Y-Axis Response: %s" % SubItem1.replace(".csv",""))

            SkipList.append(SubItem1)
            #SkipList.append(SubItem2)

            Iterator += 1


    if BigIterator == 946:
        #YInitial = YDataArray
        #XDataArray2 = [[],[]]
        #XDataArray2[0].append(XDataArray[0][0])
        #XDataArray2[1].append(XDataArray[1][0])
        
        #YDataArray2 = [[],[]]
        #YDataArray2[0].append(np.average(YDataArray[0], axis = 0))
        #YDataArray2[1].append(np.average(YDataArray[1], axis = 0))

        #YDataArray = YDataArray2
        #XDataArray = XDataArray2
        #Labels = [["No movement"],["No movement"]]
        pass
    else:
        Index = 0
        Mult2 = 5
        for elements in YDataArray[0]:
            YDataArray[0][Index] -= - Index * Mult2
            YDataArray[1][Index] -= - Index * Mult2

            #YDataArray[0][Index] = np.convolve(YDataArray[0][Index], window, mode = "valid")
            #YDataArray[1][Index] = np.convolve(YDataArray[1][Index], window, mode = "valid")
            Index += 1



    Graphing.StackedDoubleXY(XDataArray, YDataArray, Labels, AxisLabels, Title, ShareY = False, FigureName = "%s.png" % NameArray2[BigIterator])
    BigIterator += 1