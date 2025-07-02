import matplotlib.pyplot as plt
import numpy as np
import time
import matplotlib as mpl


## Do NOT Modify these EVER
### XY Plot Settings
def XY(XData, YData, Label, AxisLabels, Title):
    plt.rcParams["figure.figsize"] = (4,3)
    mpl.rcParams['figure.dpi'] = 300

    Fig, Ax = plt.subplots()

    Ax.set(xlabel = AxisLabels[0], ylabel = AxisLabels[1])
    Ax.plot(XData, YData, label = Label)

    Ax.set_title(Title)
    
    plt.show()

def DoubleXY(XData, YData, Labels, AxisLabels, Title, **kwargs):
    ## https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
    #plt.rcParams["figure.figsize"] = (4,2)
    mpl.rcParams['figure.dpi'] = 600
    plt.rc('axes', titlesize = 12, labelsize = 8)
    #plt.rc('figure', titlesize = 4)
    #plt.rc('legend', fontsize = 3)
    plt.rc('xtick', labelsize = 8)
    plt.rc('ytick', labelsize = 8)

    AxisLims = kwargs.get("AxisLimits", False)
    ShareY = kwargs.get("ShareY", True)
    ShareX = kwargs.get("ShareX", True)

    Fig, Axs = plt.subplots(figsize = (6, 3), ncols = 2, sharey = ShareY, sharex = ShareX)




    
    if AxisLims != False:
        XLimits = AxisLims[0]
        YLimits = AxisLims[1]

        AxisIndex = 0
        for Ax in Axs.flat:
            Ax.set(xlim = XLimits[AxisIndex], ylim = YLimits[AxisIndex])

            AxisIndex += 1
    
    for Ax in Axs.flat:
        Ax.set(xlabel = AxisLabels[0], ylabel = AxisLabels[1])
        #Ax.label_outer()

    #Axs[0].set(xlabel = AxisLabels[0], ylabel = AxisLabels[1])
    #Axs[1].set(xlabel = AxisLabels[0], ylabel = AxisLabels[1])
    
    Axs[0].plot(XData[0],YData[0], label = Labels[0])
    Axs[1].plot(XData[1],YData[1], label = Labels[1])

    
    Axs[0].legend(loc = "upper right")
    Axs[1].legend(loc = "upper right")

    Fig.suptitle(Title, y = .95)
    #Fig.set_figwidth = 8
    #Fig.set_figheight = 4
    
    Fig.tight_layout()

    #plt.subplots_adjust(left = 0.15, right = 1, top = 0.8, bottom = 0.25)

    plt.show()

def StackedDoubleXY(XData, YData, Labels, AxisLabels, Title, **kwargs):
    ## https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
    #plt.rcParams["figure.figsize"] = (4,2)
    mpl.rcParams['figure.dpi'] = 300
    plt.rc('axes', titlesize = 12, labelsize = 8)
    #plt.rc('figure', titlesize = 4)
    plt.rc('legend', fontsize = 4)
    plt.rc('xtick', labelsize = 8)
    plt.rc('ytick', labelsize = 8)

    AxisLims = kwargs.get("AxisLimits", False)
    ShareY = kwargs.get("ShareY", True)
    ShareX = kwargs.get("ShareX", True)
    FigureName = kwargs.get("FigureName", "Defaualt")

    Fig, Axs = plt.subplots(figsize = (6, 3), ncols = 2, sharey = ShareY, sharex = ShareX)




    
    if AxisLims != False:
        XLimits = AxisLims[0]
        YLimits = AxisLims[1]

        AxisIndex = 0
        for Ax in Axs.flat:
            Ax.set(xlim = XLimits[AxisIndex], ylim = YLimits[AxisIndex])

            AxisIndex += 1
    
    for Ax in Axs.flat:
        Ax.set(xlabel = AxisLabels[0], ylabel = AxisLabels[1])
        #Ax.label_outer()

    #Axs[0].set(xlabel = AxisLabels[0], ylabel = AxisLabels[1])
    #Axs[1].set(xlabel = AxisLabels[0], ylabel = AxisLabels[1])
    
    DataIndex = 0
    for data in YData[0]:
        #print(len(XData[0]))
        #print(XData[0][DataIndex])
        #print(YData[0][DataIndex])
        #print(Labels[0][DataIndex])
        Axs[0].plot(XData[0][DataIndex],YData[0][DataIndex], label = Labels[0][DataIndex], linewidth = 1)
        Axs[1].plot(XData[1][DataIndex],YData[1][DataIndex], label = Labels[1][DataIndex], linewidth = 1)
        DataIndex += 1

    Axs[0].legend(bbox_to_anchor = (0.5, 1.2), loc = "upper center", ncol = 2)
    Axs[1].legend(bbox_to_anchor = (0.5, 1.2), loc = "upper center", ncol = 2)
    #Axs[0].legend(loc = "upper left")
    #Axs[1].legend(loc = "upper left")

    Fig.suptitle(Title, y = .95)
    #Fig.set_figwidth = 8
    #Fig.set_figheight = 4
    
    Fig.tight_layout()

    #plt.subplots_adjust(left = 0.15, right = 1, top = 0.8, bottom = 0.25)
    plt.savefig("D:\Microscope Data\ZDRPI_XYZ_Characterization\ExtractedData\Figures/%s" % FigureName)
    plt.show()
