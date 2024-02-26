# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:56:12 2024

@author: Clement

Made to analyze the results of the simulations to compare PnET and
Biomass Succession.
"""

#%% IMPORTING MODULES

import sys, os, glob
import time
import pandas as pd
from osgeo import gdal
import numpy as np
import math
import pickle
import matplotlib.pyplot as plt
# import geopandas as geo
# import rasterio
# from rasterio.features import shapes
from tqdm import tqdm
# import numpy.ma as ma
# import matplotlib.patches as mpatches

import glob, os, re, math
from osgeo import gdal
from collections import Counter
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter
import pandas as pd
import statistics
from tqdm import tqdm
import textwrap

#%% DEFINING FUNCTIONS

def getRasterData(path):
    raster = gdal.Open(path)
    rasterData = raster.GetRasterBand(1)
    rasterData = rasterData.ReadAsArray()
    rasterData - rasterData.astype('float64')
    return(np.array(rasterData))

def getRasterDataAsList(path):
    return(getRasterData(path).tolist())

def writeNewRasterData(rasterDataArray, pathOfTemplateRaster, pathOfOutput):
    # Saves a raster in int16 with a nodata value of 0
    # Inspired from https://gis.stackexchange.com/questions/164853/reading-modifying-and-writing-a-geotiff-with-gdal-in-python
    # Loading template raster
    template = gdal.Open(pathOfTemplateRaster)
    driver = gdal.GetDriverByName("GTiff")
    [rows, cols] = template.GetRasterBand(1).ReadAsArray().shape
    outputRaster = driver.Create(pathOfOutput, cols, rows, 1, gdal.GDT_Int16)
    outputRaster.SetGeoTransform(template.GetGeoTransform())##sets same geotransform as input
    outputRaster.SetProjection(template.GetProjection())##sets same projection as input
    outputRaster.GetRasterBand(1).WriteArray(rasterDataArray)
    outputRaster.GetRasterBand(1).SetNoDataValue(0)##if you want these values transparent
    outputRaster.FlushCache() ##saves to disk!!
    outputRaster = None
    
def writeNewRasterDataFloat32(rasterDataArray, pathOfTemplateRaster, pathOfOutput):
    # Saves a raster in Float32 with a nodata value of 0.0
    # Inspired from https://gis.stackexchange.com/questions/164853/reading-modifying-and-writing-a-geotiff-with-gdal-in-python
    # Loading template raster
    template = gdal.Open(pathOfTemplateRaster)
    driver = gdal.GetDriverByName("GTiff")
    [rows, cols] = template.GetRasterBand(1).ReadAsArray().shape
    outputRaster = driver.Create(pathOfOutput, cols, rows, 1, gdal.GDT_Float32)
    outputRaster.SetGeoTransform(template.GetGeoTransform())##sets same geotransform as input
    outputRaster.SetProjection(template.GetProjection())##sets same projection as input
    outputRaster.GetRasterBand(1).WriteArray(rasterDataArray)
    outputRaster.GetRasterBand(1).SetNoDataValue(-1)##if you want these values transparent
    outputRaster.FlushCache() ##saves to disk!!
    outputRaster = None

def writeExistingRasterData(rasterDataArray, pathOfRasterToEdit):
    # Edits the data of an existing raster
    rasterToEdit = gdal.Open(pathOfRasterToEdit, gdal.GF_Write)
    rasterToEdit.GetRasterBand(1).WriteArray(rasterDataArray)
    rasterToEdit.FlushCache() ##saves to disk!!
    rasterToEdit = None

# From https://pynative.com/python-write-list-to-file/
# write list to binary file
def write_list(a_list, filePath):
    # store list in binary file so 'wb' mode
    with open(filePath, 'wb') as fp:
        pickle.dump(a_list, fp)
        print('List saved at path :' + str(filePath))

# Read list to memory
def read_list(filePath):
    # for reading also binary mode is important
    with open(filePath, 'rb') as fp:
        n_list = pickle.load(fp)
        return n_list

def weightedStandardDeviation(listOfValues, listOfWeights):
    """Taken from https://stackoverflow.com/a/65434084, and validated
    against results from  DescrStatsW(values, weights=weights).std (see https://stackoverflow.com/a/36464881)"""
    return math.sqrt(np.average((listOfValues - np.average(listOfValues, weights=listOfWeights))**2, weights=listOfWeights))

def getFamilyZoneName(rasterName):
    familyZoneName = rasterName.split("/")[-1].split("\\")[-1][0:-4]
    return(familyZoneName)


def ComputeMooseHQI(youngForestSurface,
                            coniferForestSurface,
                            mixedLeafyForestSurface,
                            WelandSurface,
                            nonForestSurface):
    """
From Koitzsch, K. B. (2002). APPLICATION OF A MOOSE HABITAT SUITABILITY INDEX
MODEL TO VERMONT WILDLIFE MANAGEMENT UNITS. Alces: A Journal Devoted to the
Biology and Management of Moose, 38, 89‑107.

Formulas are interpretations from the figure of the article via WebPlotDigitizer
and finding the equations of each parts of the curves by hand (see note 04-05122023).
    """
    totalSurface = (youngForestSurface + coniferForestSurface +
                    mixedLeafyForestSurface + WelandSurface + nonForestSurface)
    percentageYoung = (youngForestSurface / totalSurface) * 100
    percentageConifer = (coniferForestSurface / totalSurface) * 100
    percentageMixedLeafy = (mixedLeafyForestSurface / totalSurface) * 100
    percentageWetland = (WelandSurface / totalSurface) * 100
    
    if percentageYoung < 40:
        youngHQI = 0.025 * percentageYoung
    elif percentageYoung < 50:
        youngHQI = 1 
    else:
        youngHQI = -0.02*percentageYoung + 2
        
    if percentageConifer < 4.5:
        coniferHQI = 0.2222223 * percentageConifer
    elif percentageConifer < 15:
        coniferHQI = 1 
    else:
        coniferHQI = -0.0117*percentageConifer + 1.17
    
    if percentageMixedLeafy < 35:
        mixedLeafyHQI = 0.0285 * percentageMixedLeafy
    elif percentageMixedLeafy < 55:
        mixedLeafyHQI = 1 
    else:
        mixedLeafyHQI = -0.023*percentageMixedLeafy + 2.23
        
    if percentageWetland < 4.5:
        wetlandHQI = 0.222223 * percentageWetland
    elif percentageWetland < 10:
        wetlandHQI = 1 
    else:
        wetlandHQI = -0.011111*percentageWetland + 1.11111
        
    return((youngHQI * coniferHQI * mixedLeafyHQI * wetlandHQI)**(0.25))

def ComputeMooseHQIDussault2006(dictForestTypeSurface, coverFoodEdges, Percentile70thCoverFoodEdgeDensity):
    """
From Dussault et al. (2006). A habitat suitability index model to assess moose
habitat selection at multiple spatial scales. Can J For Res 36:1097–1107.
https://doi.org/10.1139/x05-310


Formulas are from Material and methods.
    """
    SIfood = ((dictForestTypeSurface["Mi10"] + dictForestTypeSurface["Dt50"] + dictForestTypeSurface["Mt50"])
              + (dictForestTypeSurface["Di50"] + dictForestTypeSurface["Mi30"]) * 0.5 
              + dictForestTypeSurface["Mi50"] * 0.4 
              + dictForestTypeSurface["C10"] * 0.3
              + dictForestTypeSurface["CF30"] * 0.15
              + dictForestTypeSurface["IMP"] * 0.1
              + dictForestTypeSurface["CS30"] * 0.05)
    
    betweenStandEdgeIndex = coverFoodEdges / Percentile70thCoverFoodEdgeDensity
    betweenStandEdge = (1 - dictForestTypeSurface["Mi50"]) * betweenStandEdgeIndex
    withinStandEdge = dictForestTypeSurface["Mi50"]
    SIedge = withinStandEdge + betweenStandEdge
        
    return(SIfood * 0.45 + SIedge * 0.55)

#%% READING DATA

os.chdir(r"D:\OneDrive - UQAM\1 - Projets\Thèse - Comparaison PnET Biomass Succession\scenarios")

# We want a dictionnary containing all of the rasters of biomass for every species
# simulated; we'll use them to make our figures.

biomassArraysResults = dict()
for extension in ["PnET", "BiomassSuccession"]:
    biomassArraysResults[extension] = dict()
    for species in ["ABIE.BAL", "BETU.PAP", "LARI.LAR", "PICE.GLA", "PICE.MAR",
                    "PINU.BAN", "POPU.TRE", "BETU.ALL"]:
        biomassArraysResults[extension][species] = dict()
        for timestep in range(0, 105, 5):
            rasterDataPath = ("./" + str(extension)
                              + "/output/biomass/bio-"
                              + str(species) + "-" + str(timestep) + ".img")
            biomassArraysResults[extension][species][timestep] = getRasterData(rasterDataPath)

#%% TRANSFORMING DATA

# Structure we want : dict[variable][first factor][second factor] and ["Mean"] or ["SD"] at the end
# First factor = extension, second factor = the row (which corresponds to the number of species in the pixels)



dictVariablesValues = dict()

# First : Evolution of total biomass

dictVariablesValues["TotalBiomass"] = dict()

for extension in ["PnET", "BiomassSuccession"]:
    dictVariablesValues["TotalBiomass"][extension] = dict()
    for numberOfSpecies in range(1, 6):
        dictVariablesValues["TotalBiomass"][extension][numberOfSpecies] = dict()
        dictVariablesValues["TotalBiomass"][extension][numberOfSpecies]["Mean"] = list()
        dictVariablesValues["TotalBiomass"][extension][numberOfSpecies]["SD"] = list()
        for timestep in range(0, 105, 5):
            biomassAllSpecies = np.zeros_like(biomassArraysResults[extension]["ABIE.BAL"][timestep])
            for species in ["ABIE.BAL", "BETU.PAP", "LARI.LAR", "PICE.GLA", "PICE.MAR",
                            "PINU.BAN", "POPU.TRE", "BETU.ALL"]:
                biomassAllSpecies += biomassArraysResults[extension][species][timestep]
                
            meanSumTimestepForRow = np.mean(biomassAllSpecies, axis = 1)[numberOfSpecies - 1]
            sdForPixelsOfRow = np.std(biomassAllSpecies, axis = 1)[numberOfSpecies - 1]
            
            dictVariablesValues["TotalBiomass"][extension][numberOfSpecies]["Mean"].append(meanSumTimestepForRow)
            dictVariablesValues["TotalBiomass"][extension][numberOfSpecies]["SD"].append(sdForPixelsOfRow)

# Second : Biomass for each function group
# It's simply the biomass per species, since we have only 5 species 
# in the initial communities
speciesList = ["ABIE.BAL", "PICE.GLA", "PINU.BAN", "BETU.ALL", "POPU.TRE"]

for species in speciesList:
    dictVariablesValues["Biomass " + str(species)] = dict()

    for extension in ["PnET", "BiomassSuccession"]:
        dictVariablesValues["Biomass " + str(species)][extension] = dict()
        for numberOfSpecies in range(1, 6):
            dictVariablesValues["Biomass " + str(species)][extension][numberOfSpecies] = dict()
            dictVariablesValues["Biomass " + str(species)][extension][numberOfSpecies]["Mean"] = list()
            dictVariablesValues["Biomass " + str(species)][extension][numberOfSpecies]["SD"] = list()
            for timestep in range(0, 105, 5):
                meanSumTimestepForRow = np.mean(biomassArraysResults[extension][species][timestep], axis = 1)[numberOfSpecies - 1]
                sdForPixelsOfRow = np.std(biomassArraysResults[extension][species][timestep], axis = 1)[numberOfSpecies - 1]
                dictVariablesValues["Biomass " + str(species)][extension][numberOfSpecies]["Mean"].append(meanSumTimestepForRow)
                dictVariablesValues["Biomass " + str(species)][extension][numberOfSpecies]["SD"].append(sdForPixelsOfRow)

# Finally : FDiv (exponential of the shannon diversity index)
# Computed for each pixels, and then averaged
# Computation based on what's done at line 527 of pythonAnalysis_chapter3_v1.2.py

dictVariablesValues["ShannonFDIV"] = dict()

for extension in ["PnET", "BiomassSuccession"]:
    dictVariablesValues["ShannonFDIV"][extension] = dict()
    for numberOfSpecies in range(1, 6):
        dictVariablesValues["ShannonFDIV"][extension][numberOfSpecies] = dict()
        dictVariablesValues["ShannonFDIV"][extension][numberOfSpecies]["Mean"] = list()
        dictVariablesValues["ShannonFDIV"][extension][numberOfSpecies]["SD"] = list()
        for timestep in range(0, 105, 5):
            listOfPixelShannonValues = list()
            for pixel in range(0, len(biomassArraysResults[extension][species][timestep][numberOfSpecies - 1])):
                biomassOfSpeciesInPixel = list()
                for species in speciesList:
                    biomassOfSpeciesInPixel.append(biomassArraysResults[extension][species][timestep][numberOfSpecies - 1][pixel])
                biomassOfSpeciesInPixel = np.array([x for x in biomassOfSpeciesInPixel if x > 0],dtype=np.int64)
                if sum(biomassOfSpeciesInPixel) == 0:
                    shannonFdivInPixel = 0
                else:
                    relativeAbundanceOfGroups = [x/sum(biomassOfSpeciesInPixel) for x in biomassOfSpeciesInPixel]
                    shannonFdivInPixel = math.exp(-sum([(x*math.log(x)) for x in relativeAbundanceOfGroups if x > 0]))
                
                listOfPixelShannonValues.append(shannonFdivInPixel)
            dictVariablesValues["ShannonFDIV"][extension][numberOfSpecies]["Mean"].append(np.mean(listOfPixelShannonValues))
            dictVariablesValues["ShannonFDIV"][extension][numberOfSpecies]["SD"].append(np.std(listOfPixelShannonValues))
                
#%% MAKING FIGURES

# Here, we want a figure for one variable for each graph
# On each graph : as many plots as for the first factor
# Each curve is the mean value for the second factor,
# with SD represented as transparent buffers around the curves

# The code comes from Script_figures_chapitre3_BasicAveraged_v3.0.py
# from my third chapter. I didn't took the time to clean it.

variables = dictVariablesValues.keys()
Linux = False

def GenerateGraphsOfResults(dictVariablesValues, 
                            variables,
                            pathToSave,
                            firstFactorName,
                            levelsOfFirstFactor,
                            secondFactorName,
                            levelsOfsecondFactor,
                            legendOff = False, 
                            saveTheGraphs = False, 
                            noXAxis = False,
                            variableSelection = [],
                            antiVariableSelection = []):
    """Generates the graphs for the density and cost results.
    Allow the graphs to be with or without a legend, with or
    without x axis, and to be displayed or saved at a given path."""
    for variable in variables:
        doThisVariable = True
        if len(variableSelection) > 0:
            for selectionItem in variableSelection:
                if selectionItem not in variable:
                    doThisVariable = False
        if len(antiVariableSelection)> 0:
            for antiSelectionItem in antiVariableSelection:
                if antiSelectionItem in variable:
                    doThisVariable = False
        if doThisVariable:
            # We create a figure with a number of subplots based on the number
            # of levels of the first factor
            numberOfColumns = len(levelsOfFirstFactor)
            numberOfLines = 1
            i = 0
            fig, axes = plt.subplots(numberOfLines, numberOfColumns, sharey= True)
            
            # We resize the figure and graphs
            ResizeFigure(fig, legendOff)
            
            # We do a graph for each level of the first factor
            legendAlreadyDrawn = False
            for levelOfFirstFactor in levelsOfFirstFactor:
                legendItems = list()
                # We get the color palette of the curves for the second factor
                colorPalette = ChooseColorPaletteAccordingToFactor(secondFactorName, levelsOfsecondFactor)
                # We plot the curves
                DrawCurvesForLevelsOfFactor(axes[i], dictVariablesValues, levelOfFirstFactor,
                                            levelsOfsecondFactor, variable, colorPalette, legendItems)
                # We customize the limits of the y-axis
                CustomGraphLimits(axes[i], dictVariablesValues, variable)
                # We make the graph more beautiful
                MakeGraphPrettier(axes[i], variable, dictVariablesValues, levelOfFirstFactor,
                                  levelsOfsecondFactor, (legendOff or legendAlreadyDrawn), legendItems,  (i == 0))
                legendAlreadyDrawn = True
                # We show the plot
                plt.show()
                # We increment the graph number
                i = i + 1
            # When all graphes have been made, we make the figure prettier.
            MakeFigurePrettier(axes, levelsOfFirstFactor, noXAxis)
            # If needed, we save the figures
            if saveTheGraphs:
                fig.savefig(pathToSave + variable + ".svg", format='svg', dpi=600, transparent=False)
                plt.close()

def FindMaximumValueForVariable(dictVariablesValues, variable):
    """Return the maximal value of a given variable in
    all of the data frames of a result dictionary."""
    listOfMaximalValues = list()
    for firstFactorLevel in dictVariablesValues[variable].keys():
        for secondFactorLevel in dictVariablesValues[variable][firstFactorLevel].keys():
            listOfMaximalValues.append(max(np.array(dictVariablesValues[variable][firstFactorLevel][secondFactorLevel]["Mean"])
                                           + np.array(dictVariablesValues[variable][firstFactorLevel][secondFactorLevel]["SD"])))
    return(max(listOfMaximalValues))

def MakeFactorBeautiful(factor):
    """Return a more beautiful version of the name of the factor."""
    if factor == "Climate":
        return "Climate"
    elif factor == "Management":
        return "Forest management"
    elif factor == "MegaFire":
        return "🔥 Mega Fire"
    elif factor == "MegaDrought":
        return "💧 Mega Drought"
    elif factor == "MountainPineBeetle":
        return "🐞 Mountain Pine Beetle"
    elif factor == "NoCatastrophy":
        return "✅ No Catastrophy"
    else:
        return factor
    
def MakeVariableBeautiful(variable):
    # To know how to write mathematical equation style letters in Matplotlib :
    # https://matplotlib.org/3.3.1/tutorials/text/mathtext.html
    activated = False
    if activated:
        beautifullVariable = ""
        if "TOTAL MATURE BIOMASS" in variable.upper():
            beautifullVariable += "Harvestable timber - total" + " (Tons)"
            if "-" in variable.upper():
                beautifullVariable += " - " + variable.upper()[-8:len(variable.upper())]
        elif "BIOMASS ECONOMIC SPECIES" in variable.upper():
            beautifullVariable += "Harvestable timber - economic species\n(Tons)"
        elif "MEAN FUNCTIONAL DIVERSITY" in variable.upper():
            beautifullVariable += "Mean Functional Diversity"
        else:
            beautifullVariable += variable.upper()
    else:
        beautifullVariable = variable
    return beautifullVariable
    
def ResizeFigure(figureObject, legendOff, specialFigure7 = False):
    """Resizes the graphs and the figure made by the script.""" 
    if specialFigure7:
        width = 3.5
        leftSubplot = 0.230
    else:
        width = 12.5
        leftSubplot = 0.080
        
    if legendOff :
        figureObject.set_size_inches(width, 3.5)
    else:
        figureObject.set_size_inches(width, 4.2)
    if legendOff:
        plt.subplots_adjust(wspace=0.15,
            top = 0.95,
            bottom = 0.25,
            left = leftSubplot,
            right = 0.990)
    else:
        plt.subplots_adjust(wspace=0.15,
                top = 0.838,
                bottom = 0.22,
                left = leftSubplot,
                right = 0.990)

def ChooseColorPaletteAccordingToFactor(secondFactorName, levelsOfsecondFactor):
        """Return the color palette needed for the curves of a given factor."""
        if secondFactorName == "Management":
            # colorsCustom = ["#4D2A29", "#F4BA2D", "#7CA982", "#33658A"] # First version by Clement
            # colorsCustom = ["#AD0328", "#55ABC6", "#DF841A", "#8EA604"] # From GeoDataViz on Github
            colorsCustom = ["#4d2a29", "#FFBA49", "#d95d39", "#33658a"] # Another version by Coolors
            
            colorPalette=iter(colorsCustom)
        else:
            colorPalette=iter(plt.cm.inferno(np.linspace(1,0,len(levelsOfsecondFactor))))
            
        return(colorPalette)
    
def findHatchType(labelForCurve):
    hatchType = "-"
    if labelForCurve == "TRIAD+":
        hatchType = "---"
    elif labelForCurve == "normal-TRIAD":
        hatchType = "|||"
    elif labelForCurve == "BAU-PlantFunct":
        hatchType = "///"
    elif labelForCurve == "BAU-NoPlant":
        hatchType = "\\\\\\"
    return(hatchType)
    
def findLinestyle(labelForCurve):
    lineStyle = "solid"
    if labelForCurve == "TRIAD+":
        lineStyle = ((0, (2, 2)))
    elif labelForCurve == "normal-TRIAD":
        lineStyle = "solid"
    elif labelForCurve == "BAU-PlantFunct":
        lineStyle = ((2, (2, 2)))
    elif labelForCurve == "BAU-NoPlant":
        lineStyle = "solid"
    return(lineStyle)

def findLinewidth(labelForCurve):
    lineWidth = 1
    if labelForCurve == "TRIAD+":
        lineWidth = 2
    elif labelForCurve == "normal-TRIAD":
        lineWidth = 2
    elif labelForCurve == "BAU-PlantFunct":
        lineWidth = 3.5
    elif labelForCurve == "BAU-NoPlant":
        lineWidth = 3.5
    return(lineWidth)

def DrawCurvesForLevelsOfFactor(axis, dictVariablesValues, firstFactorLevel,
                                levelsOfsecondFactor, variable, colorPalette, legendItems):
        """Draws the curves for a variable, a factor and levels of factor on
        a given axis."""
        # We get the right data frame for the factor/factor level, then the
        # list of timesteps, and then we plot.
        for levelOfsecondFactor in levelsOfsecondFactor:
            meanResults = dictVariablesValues[variable][firstFactorLevel][levelOfsecondFactor]["Mean"]
            SDResults = dictVariablesValues[variable][firstFactorLevel][levelOfsecondFactor]["SD"]
            sequenceOfYears = range(0, 105, 5)
            colorOfCurve = next(colorPalette)
            # We get the label for the curve, which is used to make the legend
            # EDIT : We change TRIAD+++ to TRIAD+ here, because data already uses TRIAD+++.
            labelForCurve = levelOfsecondFactor
            # axis.errorbar(sequenceOfYears, meanResults, SDResults, linewidth=1, marker='o', color = colorOfCurve)
            # We plot the standard deviation around the curve with same color,
            # but less opacity
            hatchType = findHatchType(labelForCurve)
            lineStyle = findLinestyle(labelForCurve)
            lineWidth = findLinewidth(labelForCurve)
            f = axis.fill_between(sequenceOfYears, list(np.array(meanResults) + np.array(SDResults)),
                     list(np.array(meanResults) - np.array(SDResults)), color=colorOfCurve,
                     alpha=0.3, label = labelForCurve, edgecolor=(0,0,0,0), linewidth=0.0)
            # f2 = axis.fill_between(sequenceOfYears, list(np.array(meanResults) + np.array(SDResults)),
            #          list(np.array(meanResults) - np.array(SDResults)), color=colorOfCurve,
            #          facecolor="none", hatch= hatchType, edgecolor=(0,0,0,0.06), linewidth=0.0)
            # We plot the curve
            p = axis.plot(sequenceOfYears, meanResults,
                     label = labelForCurve, color = colorOfCurve, linewidth = lineWidth, linestyle = lineStyle)
            axis.set_xticks([0,25,50,75,100],
                            ["0","25","50","75","100"])
            legendItems.append((f, p[0]))

def CustomGraphLimits(axis, dictVariablesValues, variable, specialFigure7 = False):
    """Given the axis of a graph and a variable, sets the right limits
    on the y axis to better display the variation of the variable."""
    if "clumpy" in variable: # Varies from -1 to 1; but we start him at 0
        axis.set_ylim(0,1)
    elif "tca" in variable: # Custom
        if specialFigure7:
            axis.set_ylim(0,100)
        else:
            axis.set_ylim(0,50)
    # If no custom limits, the limits are just from zero to the maximum value
    # for the variable, with a margin of 5%
    else:
        axis.set_ylim(bottom = 0, top = FindMaximumValueForVariable(dictVariablesValues, variable) * 1.05)
        
def MakeGraphPrettier(axis, variable, dictVariablesValues, levelOfFirstFactor,
                      levelsOfsecondFactor, legendOff, legendItems, leftGraph = False):
    """Changes the font size and various aspects
    of a graph to make it look better.
    Also remove the y-axis ticks on graphs that are not the first
    on the left of the figure, to avoid repeating the axis."""
    # Title format and presence if legend is activated
    # if not legendOff:
    titleToPut = MakeFactorBeautiful(levelOfFirstFactor)
    axis.set_title(titleToPut, fontsize = 16, y = -0.32)
    # Global font size
    plt.rcParams.update({'font.size': 12}) # Change the size of the fonts, except the title
    # Tick labels
    axis.tick_params(reset = True, axis='both', which='major', labelsize= 13, right = False, top = False)
    # Ewriting the timestep
    axis.set_xlabel("Time step", fontsize= 15)
    # Legend position and size
    if not legendOff:
        # Number of columsn in legend depend on number of levels of second factor
        numberOflevelsOfFactor = len(levelsOfsecondFactor)
        columnsNumber = numberOflevelsOfFactor
        # if numberOflevelsOfFactor > 2:
        #     columnsNumber = 2
        axis.legend(legendItems, levelsOfsecondFactor,
                    prop={'size': 12}, loc = "upper center",  bbox_to_anchor = (0, -0.01, 1, 1), bbox_transform = plt.gcf().transFigure,
                    ncol=columnsNumber, handleheight = 2, handlelength = 3)
    # No labels for graphs not at the left of the figure
    if leftGraph:
        axis.set_ylabel(MakeVariableBeautiful(variable), fontsize= 15, wrap=True)
        # If the variable goes above 1000 => scientific notation, 
        # unless it's the cost in dollars
        if FindMaximumValueForVariable(dictVariablesValues, variable) >= 1000 and "Estimated cost" not in variable:
            plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))     
            axis.yaxis.offsetText.set_visible(True)
    # Remove offset text on the other graphs
    elif FindMaximumValueForVariable(dictVariablesValues, variable) >= 1000 and not leftGraph:
        axis.yaxis.offsetText.set_visible(False)
    # We limit the number of ticks on the y-axis
    axis.yaxis.set_major_locator(ticker.MaxNLocator(4))
    
def MakeFigurePrettier(axes, levelsOfFirstFactor, noXAxis = False, specialFigure7 = False):
    """Removes y-axis tick labels on graphs that are not left of the figure,
    the top and right border of the graphs, and also the x axis if needed."""
    # Remove the unwanted ticks labels on graphs that are not left
    if not specialFigure7:
        for i in range(1, len(levelsOfFirstFactor)):
            plt.setp(axes[i].get_yticklabels(), visible=False)
            axes[i].spines['left'].set_visible(False)
            axes[i].yaxis.set_ticks_position('none')
    # Removing top and right border of the graphs to avoid "box" effect
    for i in range(0, len(levelsOfFirstFactor)):
        axes[i].spines['right'].set_visible(False)
        axes[i].spines['top'].set_visible(False)
        # We add the grid for better visibility
        axes[i].grid(color = "gainsboro", linestyle='-', linewidth=0.5) 
        axes[i].set_facecolor('#FCFCFC')
        axes[i].set_axisbelow(True)
        # Removing the X axis if needed
        if noXAxis == True:
            axes[i].spines['bottom'].set_visible(False)
            axes[i].xaxis.set_major_formatter(NullFormatter())
            axes[i].xaxis.set_ticklabels([])
            axes[i].xaxis.set_ticks_position('none')

def GetCustomLevelsOfFactors(resultDict, factor, variable):
    """Returns an list of customized levels of factors and
    a customized variable, to adapt to the custom figure for fragmentation
    indices."""
    variableToUse = variable
    redisgnedFactor = factor
    if "Percentage of biomass harvested with Uneven-aged methods" in factor and "North" in factor:
        redisgnedFactor = "Percentage of biomass harvested with Uneven-aged methods"
        if "NORTH" not in variable.upper() and "SOUTH" not in variable.upper() and "ALLFOREST" not in variable.upper():
            variableToUse = variable + "_north"
    elif "Percentage of biomass harvested with Uneven-aged methods" in factor and "South" in factor:
        redisgnedFactor = "Percentage of biomass harvested with Uneven-aged methods"
        if "SOUTH" not in variable.upper() and "NORTH" not in variable.upper() and "ALLFOREST" not in variable.upper():
            variableToUse = variable + "_south"
    elif "Aggregation of the cuts" in factor and "North" in factor:
        redisgnedFactor = "Aggregation of the cuts"
        if "SOUTH" not in variable.upper() and "NORTH" not in variable.upper() and "ALLFOREST" not in variable.upper():
            variableToUse = variable + "_north"
    elif "Aggregation of the cuts" in factor and "South" in factor:
        redisgnedFactor = "Aggregation of the cuts"
        if "SOUTH" not in variable.upper() and "NORTH" not in variable.upper() and "ALLFOREST" not in variable.upper():
            variableToUse = variable + "_south"
    return(redisgnedFactor, resultDict[redisgnedFactor], variableToUse)
    


pathToSave = r"D:\OneDrive - UQAM\1 - Projets\Thèse - Comparaison PnET Biomass Succession\FiguresResults\\"

GenerateGraphsOfResults(dictVariablesValues, 
                        variables,
                        pathToSave,
                        "Extension",
                        ["PnET", "BiomassSuccession"],
                        "Number Of Species Cohabitating",
                        range(1, 6),
                        legendOff = False, 
                        saveTheGraphs = True, 
                        noXAxis = False,
                        variableSelection = [],
                        antiVariableSelection = [])