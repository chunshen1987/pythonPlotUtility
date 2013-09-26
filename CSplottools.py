#! /usr/bin/env python

from numpy import *

def getPlotElements(idx):
   """
   return line style, marker type, and color for curve. It loops over the whole 
   lists of options.
   """
   lineStyleList = ['-', '--', '-.', ':']
   colorList = ['k', 'r', '#009900', 'b', '#990099', '#CC6600', '#009999', '#FF8000']
   shadowColorList = ['#C0C0C0', '#FFA8A8', '#B3FFB8', '#A4CCFF', '#EF9BFF', '#FFDE80', '#99FFFF', '#FFB266']
   MarkerList = ['s', 'o', '^', 'v', 'D', 'p', '*', 'H', '.', ',', '<', 
                 '>', '1', '2', '3', '4', 'h', '+', 'x', 'd', '|', '_']
   return lineStyleList[idx%len(lineStyleList)], MarkerList[idx%len(MarkerList)], colorList[idx%len(colorList)], shadowColorList[idx%len(shadowColorList)]

def getBinnedAveragedDatawithErrorbars(dataMatrix, nbin, bincol=0):
   """
   Return the binned average of data, together with the count of the number 
   of data and normalized probability distribution. It returns another matrix
   for the statistical error bars for each binned quantities
   """
   if dataMatrix.ndim == 1:
      ncol = 1
      dataMatrix = dataMatrix.reshape(len(dataMatrix), 1)
   else:
      ncol = len(dataMatrix[0, :])
   binColumn = dataMatrix[:, bincol]
   ntotal = len(binColumn)
   
   binMin = min(binColumn); binMax = max(binColumn)+1e-8
   binEdges = linspace(binMin, binMax, nbin+1)
   
   binnedData = zeros([nbin, ncol+2])
   binnedData_err = zeros([nbin, ncol])

   for idxbin in range(nbin):
      binWidth = binEdges[idxbin+1] - binEdges[idxbin]
      binmid = (binEdges[idxbin+1] + binEdges[idxbin])/2.
      idxdata = logical_and(binColumn >= binEdges[idxbin], binColumn < binEdges[idxbin+1])
      nsamples = len(dataMatrix[idxdata, bincol])
      for icol in range(ncol):
         if nsamples == 0: 
            binnedData[idxbin, icol] = 0
            binnedData_err[idxbin, icol] = 0
         else:
            binnedData[idxbin, icol] = mean(dataMatrix[idxdata, icol])
            binnedData_err[idxbin, icol] = std(dataMatrix[idxdata, icol])/sqrt(nsamples)
      if nsamples == 0:
         binnedData[idxbin, bincol] = binmid
      binnedData[idxbin, ncol] = nsamples
      binnedData[idxbin, ncol+1] = nsamples/binWidth/ntotal
   return binnedData, binnedData_err

