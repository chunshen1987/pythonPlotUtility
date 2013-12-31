#! /usr/bin/env python

from sys import argv, exit
from os import path, remove
from DBR import SqliteDB
from numpy import *
from minbiasEccReader import minbiasEccReader

#define colors
purple = "\033[95m"
green = "\033[92m"
red = "\033[91m"
normal = "\033[0m"

def printHelpMessageandQuit():
    print "Usage : "
    print "minbias_CentralityCut databaseName -type cutType"
    print "Usage of minbias_CentralityCut command line arguments: "
    print "-type   the type of quantity used to cut centrality: " + green + " Npart, Ncoll, b, total_entropy" + normal
    print "-h | -help   This message"
    exit(0)

if __name__ == "__main__":
    multiplicityFactor = 1.0
    if len(argv) <= 2:
        printHelpMessageandQuit()
    dbName = str(argv[1])
    while len(argv) > 2:
        option = argv[2]; del argv[2]
        if option == '-type':
            cutType = str(argv[2]); del argv[2]
            if not cutType in ['Npart', 'Ncoll', 'b', 'total_entropy']:
                print argv[0], ": invalid cutType", red + cutType + normal
                printHelpMessageandQuit()
        elif option == '-mult':
            multiplicityFactor = float(argv[2]); del argv[2]
        elif option == '-h' or option == '-help':
            printHelpMessageandQuit()
        else:
            print argv[0], ': invalid option', option
            printHelpMessageandQuit()
    reader = minbiasEccReader(dbName)
    reader.cutCentralitieswitheccStatistics(cutType, multiplicityFactor)
