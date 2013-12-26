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

class particleReader(object):
    """
        This class is used to perform statistical analysis on particle database from UrQMD.
    """
    def __init__(self, database):
        """
            Register a sqlite database
        """
        # setup database
        if isinstance(database, str):
            if path.exists(database):
                database = SqliteDB(database)
            else:
                raise ValueError("EbeDBReader.__init__: the input argument must be an existing database file.")
        if isinstance(database, SqliteDB):
            self.db = database
        else:
            raise TypeError("EbeDBReader.__init__: the input argument must be a string or a SqliteDB database.")
        # setup lookup tables
        self.pid_lookup = dict(self.db.selectFromTable("pid_lookup"))

        # set self.hasInitializedStringSubstitution to none for lazy initialization in evaluateExpression function
        self.hasInitializedStringSubstitution = False

    def getNumberOfHydroEvents(self):
        """
            return total number hydro events stored in the database
        """
        Nev = self.db._executeSQL("select count() from (select distinct hydroEvent_id from particle_list)").fetchall()[0][0]
        return(Nev)

    def getNumberOfUrQMDEvents(self, hydroEventid):
        """
            return number of UrQMD events for the given hydro event
        """
        Nev = self.db._executeSQL("select count() from (select distinct UrQMDEvent_id from particle_list where hydroEvent_id = %d)" % hydroEventid).fetchall()[0][0]
        return(Nev)

    def getNumberOftotalEvents(self):
        """
            return total number of events stored in the database
        """
        hydroEventid = self.db._executeSQL("select distinct hydroEvent_id from particle_list").fetchall()
        Nev = 0
        for iev in range(len(hydroEventid)):
           Nev += self.getNumberOfUrQMDEvents(hydroEventid[iev][0])
        return(Nev)


    def collectParticleSpectrum(self, particleName="pion_p", pT_range=linspace(0, 3, 20)):
        """
            return event averaged particle spectrum (pT, dN/(dydpT))
            event loop over all the hydro + UrQMD events
        """
        pid = self.pid_lookup[particleName]
        pTcenter = []
        dNdydpT = []
        dNdydpTerr = []
        Nev = self.getNumberOftotalEvents()
        for ipT in range(len(pT_range)-1):
            pTlow = pT_range[ipT]
            pThigh = pT_range[ipT+1]
            dpT = pThigh - pTlow
            pTcenter.append((pTlow + pThigh)/2.)
            whereClause = "pid = %d and %g<=pT and pT<=%g" % (pid, pTlow, pThigh)
            deltaN = self.db.selectFromTable("particle_list", "count()", whereClause=whereClause)[0][0]
            dNdydpT.append(deltaN/dpT)
            dNdydpTerr.append(sqrt(deltaN/dpT)/sqrt(Nev))

        return(array([pTcenter, dNdydpT, dNdydpTerr]).transpose())
     
    def getParticleSpectrum(self, particleName = 'pion_p', pT_range = linspace(0,3,20)):
        """
            check particle spectrum is pre-collected or not. If not, collect particle spectrum 
            into the database. It returns event averaged particle spectrum with pT_range
            required by the users.
        """
        #check particle spectrum table is exist or not
        if self.db.createTableIfNotExists("particleSpectra", (("pid","integer"), ("pT","real"), ("dNdydpT", "real"), ("dNdydpTerr", "real"))):
            for aParticle in ['pion_p', 'kaon_p', 'proton']:
                pid = self.pid_lookup[aParticle]
                dNdata = self.collectParticleSpectrum(aParticle)
                print(pid, dNdata)
                for idx in range(len(dNdata[:,0])):
                    self.db.insertIntoTable("particleSpectra", (pid, dNdata[idx,0], dNdata[idx,1], dNdata[idx,2]))
        
        # close connection to commit changes
        self.db.closeConnection()
        

def printHelpMessageandQuit():
    print "Usage : "
    print "particleReader.py databaseName"
    exit(0)

if __name__ == "__main__":
    if len(argv) < 2:
        printHelpMessageandQuit()
    test = particleReader(str(argv[1]))
    test.getParticleSpectrum()

