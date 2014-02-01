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
    print "minbias_distribution outputs" + purple + " X, P(X), P(X)_err, dP(X)/dX, dP(X)/dX_err" + normal
    print "Usage : "
    print "minbias_distribution %s -disType disType -nbin nbin -cutType cutType -cen centralityBoundary" % (purple + 'databaseName' + normal)
    print "Usage of minbias_CentralityCut command line arguments: "
    print "-disType   the type of quantity of the output distribution: " + green + " Npart, Ncoll, b, total_entropy, ecc_n, deformed_cosTheta1" + normal
    print "-nbin      set number of bins "
    print "-cutType   the type of quantity used to cut centrality: " + green + " Npart, Ncoll, b, total_entropy" + normal
    print "-cen       the boundary of centrality cut: " + purple + " e.g. 0-5" + normal
    print "-h | -help   This message"
    exit(0)

if __name__ == "__main__":
    # set default values
    cutType = 'total_entropy'
    cen = [0, 100]
    nbin = 30
    if len(argv) <= 2:
        printHelpMessageandQuit()
    dbName = str(argv[1])
    while len(argv) > 2:
        option = argv[2]; del argv[2]
        if option == '-disType':
            disType = str(argv[2]); del argv[2]
            if not disType in ['Npart', 'Ncoll', 'b', 'total_entropy'] and not 'ecc' in disType and not 'deformed' in disType:
                print argv[0], ": invalid disType", red + disType + normal
                printHelpMessageandQuit()
        elif option == '-cutType':
            cutType = str(argv[2]); del argv[2]
            if not cutType in ['Npart', 'Ncoll', 'b', 'total_entropy']:
                print argv[0], ": invalid cutType", red + cutType + normal
                printHelpMessageandQuit()
        elif option == '-cen':
            cen = str(argv[2]).split('-'); del argv[2]
            cen = [float(x) for x in cen]
        elif option == '-nbin':
            nbin = int(argv[2]); del argv[2]
        elif option == '-h' or option == '-help':
            printHelpMessageandQuit()
        else:
            print argv[0], ': invalid option', option
            printHelpMessageandQuit()
    reader = minbiasEccReader(dbName)
    print(reader.getDistribution(disType, nbin, cutType, centralityBound = cen))
