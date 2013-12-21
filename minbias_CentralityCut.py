#! /usr/bin/env python

from sys import argv, exit
from os import path, remove
from DBR import SqliteDB
from numpy import *

#define centrality boundaries
centralityBoundaries = [(0,0.2),(0,1),(0,5),(5,10),(10,20),(20,30),
                       (30,40),(40,50),(50,60),(60,70),(70,80)]

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

def cutCentralities(dbName, cutType, multiplicityFactor = 1.0):
    if not path.exists(path.abspath(dbName)):
        print "Error: can not find database file: %s" % (red + dbName + normal)
        exit(1)
    database = SqliteDB(path.abspath(dbName))
    nevent = database._executeSQL("select count() from collisionParameters").fetchall()[0][0]

    for icen in range(len(centralityBoundaries)):
        lowerbound = centralityBoundaries[icen][0]
        upperbound = centralityBoundaries[icen][1]
        nsample = int(nevent*(upperbound - lowerbound)/100)-1
        noffset = int(nevent*lowerbound/100)
        fetchedData = array(database._executeSQL("select Npart, b, total_entropy from collisionParameters order by -%s limit %d offset %d" % (cutType, nsample, noffset)).fetchall())
        cenCentral = (upperbound + lowerbound)/2.
        Npartmean = mean(fetchedData[:,0])
        Npartmin = min(fetchedData[:,0]); Npartmax = max(fetchedData[:,0])
        bmin = min(fetchedData[:,1]); bmax = max(fetchedData[:,1])
        dSdymin = min(fetchedData[:,2])/multiplicityFactor; dSdymax = max(fetchedData[:,2])/multiplicityFactor
        print(cenCentral, Npartmean, Npartmin, Npartmax, dSdymin, dSdymax, bmin, bmax)
        fetchedData = array(database._executeSQL("select ecc_id, n, ecc_real, ecc_imag from eccentricities where event_id in (select event_id from collisionParameters order by -collisionParameters.%s limit %d offset %d)" % (cutType, nsample, noffset)).fetchall())
        for eccType in range(1,3):
            tempidx = (fetchedData[:,0] == eccType)
            tempdata = fetchedData[tempidx, :]
            for iorder in range(1,10):
                idx = (tempdata[:,1] == iorder)
                eccn2 = sqrt(mean(tempdata[idx,2]**2 + tempdata[idx,3]**2))
                eccn2err = std(tempdata[idx,2]**2 + tempdata[idx,3]**2)/(2.*eccn2)/sqrt(nsample)
                print(cenCentral, eccType, iorder, eccn2, eccn2err)

    
    database.closeConnection()

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
    print 'cutting centralities from database %s according to %s ....' % (purple + dbName + normal, green + cutType + normal)
    cutCentralities(dbName, cutType, multiplicityFactor)
