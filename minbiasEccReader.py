#! /usr/bin/env python

from sys import argv, exit
from os import path, remove
from DBR import SqliteDB
from numpy import *

#define colors
purple = "\033[95m"
green = "\033[92m"
red = "\033[91m"
normal = "\033[0m"

class minbiasEccReader(object):
    """
        This class contains functions to perform statistical analysis on initial eccentricities
        from minimum bias events generated from superMC
    """
    def __init__(self, databaseName):
        """
            Register a sqlite database
        """
        # setup database
        if isinstance(databaseName, str):
            if path.exists(databaseName):
                database = SqliteDB(databaseName)
                self.dbName = databaseName
            else:
                raise ValueError("EbeDBReader.__init__: the input argument must be an existing database file. %s can not be found" % (red + databaseName + normal))
        if isinstance(database, SqliteDB):
            self.db = database
        else:
            raise TypeError("EbeDBReader.__init__: the input argument must be a string or a SqliteDB database.")

        #define centrality boundaries
        self.centralityBoundaries = [(0,0.2),(0,1),(0,5),(5,10),(10,20),(20,30),
                                     (30,40),(40,50),(50,60),(60,70),(70,80)]
      
    def cutCentralitieswitheccStatistics(self, cutType, multiplicityFactor = 1.0):
        """
            this function cut the centralities and also output event averaged ecc_n with 
            statistical error
        """
        print 'cutting centralities from database %s according to %s ....' % (purple + self.dbName + normal, green + cutType + normal)
        centralityOutput = open('centralityCut_%s.dat' % cutType, 'w')
        eccnStatedOutput = open('eccnStatistics_ed_%s.dat' % cutType, 'w')
        eccnStatsdOutput = open('eccnStatistics_sd_%s.dat' % cutType, 'w')
        nevent = self.db._executeSQL("select count() from collisionParameters").fetchall()[0][0]

        for icen in range(len(self.centralityBoundaries)):
            lowerbound = self.centralityBoundaries[icen][0]
            upperbound = self.centralityBoundaries[icen][1]
            nsample = int(nevent*(upperbound - lowerbound)/100)-1
            noffset = int(nevent*lowerbound/100)
            fetchedData = array(self.db._executeSQL("select Npart, b, total_entropy from collisionParameters order by -%s limit %d offset %d" % (cutType, nsample, noffset)).fetchall())
            cenCentral = (upperbound + lowerbound)/2.
            Npartmean = mean(fetchedData[:,0])
            Npartmin = min(fetchedData[:,0]); Npartmax = max(fetchedData[:,0])
            bmin = min(fetchedData[:,1]); bmax = max(fetchedData[:,1])
            dSdymin = min(fetchedData[:,2])/multiplicityFactor; dSdymax = max(fetchedData[:,2])/multiplicityFactor
            centralityOutput.write("%6.4f  %d  %d  %d  %18.8e  %18.8e  %18.8e  %18.8e \n" % (cenCentral, Npartmean, Npartmin, Npartmax, dSdymin, dSdymax, bmin, bmax))
            fetchedData = array(self.db._executeSQL("select ecc_id, n, ecc_real, ecc_imag from eccentricities where event_id in (select event_id from collisionParameters order by -collisionParameters.%s limit %d offset %d)" % (cutType, nsample, noffset)).fetchall())
            for eccType in range(1,3):
                tempidx = (fetchedData[:,0] == eccType)
                tempdata = fetchedData[tempidx, :]
                eccOutput = []
                for iorder in range(1,10):
                    idx = (tempdata[:,1] == iorder)
                    eccn2 = sqrt(mean(tempdata[idx,2]**2 + tempdata[idx,3]**2))
                    eccn2err = std(tempdata[idx,2]**2 + tempdata[idx,3]**2)/(2.*eccn2)/sqrt(nsample)
                    eccOutput += [eccn2, eccn2err]
                if eccType == 1:
                   eccnStatsdOutput.write("%6.4f  " % cenCentral)
                   for tempecc in eccOutput:
                      eccnStatsdOutput.write("%18.8e  " % tempecc)
                   eccnStatsdOutput.write("\n")
                elif eccType == 2:
                   eccnStatedOutput.write("%6.4f  " % cenCentral)
                   for tempecc in eccOutput:
                      eccnStatedOutput.write("%18.8e  " % tempecc)
                   eccnStatedOutput.write("\n")

        centralityOutput.close()
        eccnStatsdOutput.close()
        eccnStatedOutput.close()

if __name__ == "__main__":
    reader = minbiasEccReader(str(argv[1]))
    reader.cutCentralitieswitheccStatistics('total_entropy', 1.0)
