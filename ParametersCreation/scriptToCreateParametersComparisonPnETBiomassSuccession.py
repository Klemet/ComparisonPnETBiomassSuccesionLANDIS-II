# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 11:09:13 2024

@author: Clement

Script used to create the few parameters files needed to be made "by hand" for
the simulations to compare the dynamic of vegetation with PnET and Biomass Succession.
"""

#%% IMPORTING MODULES

import sys, os, csv, json
import pandas as pd
from osgeo import gdal
from osgeo import ogr
import numpy as np
from tqdm import tqdm
# import math
import statistics
# from collections import Counter
import random
# import itertools
import shutil
import pickle

#%% FUNCTIONS

def getRasterData(path):
    raster = gdal.Open(path)
    rasterData = raster.GetRasterBand(1)
    rasterData = rasterData.ReadAsArray()
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
    outputRaster.GetRasterBand(1).SetNoDataValue(0)##if you want these values transparent
    outputRaster.FlushCache() ##saves to disk!!
    outputRaster = None
    
def writeNewRasterDataWithoutTemplate(rasterDataArray, pathOfOutput):
    # Saves a raster in int16 with a nodata value of 0
    # Inspired from https://gis.stackexchange.com/questions/164853/reading-modifying-and-writing-a-geotiff-with-gdal-in-python
    # Loading template raster
    driver = gdal.GetDriverByName("GTiff")
    [rows, cols] = rasterDataArray.shape
    outputRaster = driver.Create(pathOfOutput, cols, rows, 1, gdal.GDT_Int16)
    outputRaster.GetRasterBand(1).WriteArray(rasterDataArray)
    outputRaster.GetRasterBand(1).SetNoDataValue(0)##if you want these values transparent
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

#%% READING NEEDED FILES

os.chdir(r"D:\OneDrive - UQAM\1 - Projets\Thèse - Comparaison PnET Biomass Succession\ParametersCreation")

#%% CREATING RASTERS

# These rasters don't need projections - use writeNewRasterDataWithoutTemplate



### ECOREGION RASTER

# Just a square with the same number in all pixels.

# We want the rasters to have 5 rows and 500 columns.

# For PnET, we put the ID 102 everywhere
ecoregionsPnETArray = np.full((5, 500), 102)

# For BiomassSuccession, we put the ID 263 everywhere
ecoregionsBioSuccessArray = np.full((5, 500), 236)

# We save the rasters
writeNewRasterDataWithoutTemplate(ecoregionsPnETArray,
                                  "./ecoregionsPnET.tif")
writeNewRasterDataWithoutTemplate(ecoregionsBioSuccessArray,
                                  "./ecoregionsBioSuccession.tif")

### INITIAL COMMUNITIES RASTER

# A difference code for every pixel
initialCommunitiesArray = np.zeros((5, 500))
for row in range(0, 5):
    rowID = row * 1000
    for column in range (0, 500):
        initialCommunitiesArray[row][column] = rowID + column + 1
        
writeNewRasterDataWithoutTemplate(initialCommunitiesArray,
                                  "./initialCommunitiesPnetAndBioSuccession.tif")

#%% CREATING TEXTE FILES

### INITIAL COMMUNITIES FILES

# The 5 species we're going to use are : 
speciesToUse = ["ABIE.BAL", "PICE.GLA", "PINU.BAN", "BETU.ALL", "POPU.TRE"]

# Making the permutations for each map code of each row.
# First, we get the list of mapcodes for each rows
# We separate by rows as we'll want to change the number of species planted
mapCodesByRows = list()
for row in range(0, len(initialCommunitiesArray)):
    mapCodesByRows.append(range(int(initialCommunitiesArray[row][0]), int(initialCommunitiesArray[row][-1] + 1)))

# We fill a dictionnary
mapCodeSpeciesDict = dict()

for row in range(0, len(mapCodesByRows)):
    nbOfSpeciesToPlant = row + 1
    for pixel in mapCodesByRows[row]:
        mapCodeSpeciesDict[pixel] = random.sample(speciesToUse, nbOfSpeciesToPlant)
        
# Then, we write the text file
linesOfInitialCommunitiesTextFile = ["LandisData   \"Initial Communities\"\n\n"]
for pixel in mapCodeSpeciesDict.keys():
    linesOfInitialCommunitiesTextFile.append("MapCode " + str(pixel) +"\n")
    for species in mapCodeSpeciesDict[pixel]:
        linesOfInitialCommunitiesTextFile.append(str(species) + " 10\n")
    linesOfInitialCommunitiesTextFile.append("\n")

with open("./initialCommunitiesComparisonPnETBioSuccession.txt", 'w') as file:
    file.writelines(linesOfInitialCommunitiesTextFile)