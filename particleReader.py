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

        # define all charged hadrons
        self.chargedHadronList = ["pion_p", "pion_m", 
                                  "kaon_p", "kaon_m", 
                                  "proton", "anti_proton",
                                  "sigma_p", "sigma_m", "anti_sigma_p", "anti_sigma_m",
                                  "xi_m", "anti_xi_m"]

    def getNumberOfHydroEvents(self):
        """
            return total number hydro events stored in the database
        """
        Nev = self.db._executeSQL("select count(*) from (select distinct hydroEvent_id from particle_list)").fetchall()[0][0]
        return(Nev)

    def getNumberOfUrQMDEvents(self, hydroEventid):
        """
            return number of UrQMD events for the given hydro event
        """
        Nev = self.db._executeSQL("select count(*) from (select distinct UrQMDEvent_id from particle_list where hydroEvent_id = %d)" % hydroEventid).fetchall()[0][0]
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

    def collectParticleSpectrum(self, particleName="pion_p", rapidity_range = [-0.5, 0.5], pseudorap_range = [-0.5, 0.5]):
        """
            return event averaged particle spectrum (pT, dN/(dydpT), dN/(detadpT))
            event loop over all the hydro + UrQMD events
        """
        pT_range = linspace(0, 3, 31)
        pid = self.pid_lookup[particleName]
        if pid == 1:    # all charged hadrons
            pidString = ""
            for aPart in self.chargedHadronList:
                pidString += "pid = %d or " % self.pid_lookup[aPart]
            pidString += "pid = 1"
        else:
            pidString = "pid = %d" % pid
        pTavg = []
        dNdydpT = []
        dNdydpTerr = []
        dNdetadpT = []
        dNdetadpTerr = []
        Nev = self.getNumberOftotalEvents()
        for ipT in range(len(pT_range)-1):
            pTlow = pT_range[ipT]
            pThigh = pT_range[ipT+1]
            dpT = pThigh - pTlow
            #fetch dN/dydpT data
            whereClause = "(%s)" % pidString
            whereClause += " and (%g<=pT and pT<=%g) and (%g<=rapidity and rapidity<=%g)" % (pTlow, pThigh, rapidity_range[0], rapidity_range[1])
            tempdata = self.db.selectFromTable("particle_list", "count(*), avg(pT)", whereClause=whereClause)
            deltaN = tempdata[0][0]
            if tempdata[0][1] == None:
                pTavg.append((pTlow + pThigh)/2.)
            else:
                pTavg.append(tempdata[0][1])
            dNdydpT.append(deltaN/dpT/Nev)
            dNdydpTerr.append(sqrt(deltaN/dpT/Nev)/sqrt(Nev))
            #fetch dN/detadpT data
            whereClause = "(%s)" % pidString
            whereClause += " and (%g<=pT and pT<=%g) and (%g<=pseudorapidity and pseudorapidity<=%g)" % (pTlow, pThigh, pseudorap_range[0], pseudorap_range[1])
            tempdata = self.db.selectFromTable("particle_list", "count(*), avg(pT)", whereClause=whereClause)
            deltaN = tempdata[0][0]
            dNdetadpT.append(deltaN/dpT/Nev)
            dNdetadpTerr.append(sqrt(deltaN/dpT/Nev)/sqrt(Nev))

        return(array([pTavg, dNdydpT, dNdydpTerr, dNdetadpT, dNdetadpTerr]).transpose())
    
    def collectBasicParticleSpectra(self):
        """
            collect particle spectra into database for commonly interested hadrons
        """
        BasicParticleList = ['pion_p', 'kaon_p', 'proton']
        for aParticle in BasicParticleList:
            pid = self.pid_lookup[aParticle]
            dNdata = self.collectParticleSpectrum(aParticle)
            for idx in range(len(dNdata[:,0])):
                self.db.insertIntoTable("particleSpectra", (pid, dNdata[idx,0], dNdata[idx,1], dNdata[idx,2], dNdata[idx,3], dNdata[idx,4]))
            self.db._dbCon.commit()  # commit changes
     
    def getParticleSpectrum(self, particleName = 'pion_p', pT_range = linspace(0,3,20), rapidity_range=[-0.5, 0.5]):
        """
            check particle spectrum is pre-collected or not. If not, collect particle spectrum 
            into the database. It returns event averaged particle spectrum with pT_range and 
            eta_range required by the users.
        """
        eps = 1e-15
        pid = self.pid_lookup[particleName]
        #check particle spectrum table is exist or not
        if self.db.createTableIfNotExists("particleSpectra", (("pid","integer"), ("pT","real"), ("dNdydpT", "real"), ("dNdydpTerr", "real"), ("dNdetadpT", "real"), ("dNdetadpTerr", "real") )):
            self.collectBasicParticleSpectra()
        dNdata = array(self.db._executeSQL("select pT, dNdydpT, dNdydpTerr, dNdetadpT, dNdetadpTerr from particleSpectra where pid = %d " % pid).fetchall())
        if(dNdata.size == 0):
            dNdata = self.collectParticleSpectrum(particleName)
            if (sum(abs(dNdata[:,1])) <  1e-15):
                print "There is no record of particle: %s in the database" % particleName
                return None
            else:
                for idx in range(len(dNdata[:,0])):
                    self.db.insertIntoTable("particleSpectra", (pid, dNdata[idx,0], dNdata[idx,1], dNdata[idx,2], dNdata[idx,3], dNdata[idx,4]))
                self.db._dbCon.commit()  # commit changes

        #interpolate results to desired pT range
        dNdyinterp = exp(interp(pT_range, dNdata[:,0], log(dNdata[:,1]+eps)))
        dNdyinterp_err = exp(interp(pT_range, dNdata[:,0], log(dNdata[:,2]+eps)))
        dNdetainterp = exp(interp(pT_range, dNdata[:,0], log(dNdata[:,3]+eps)))
        dNdetainterp_err = exp(interp(pT_range, dNdata[:,0], log(dNdata[:,4]+eps)))
        results = array([pT_range, dNdyinterp, dNdyinterp_err, dNdetainterp, dNdetainterp_err])
        return(transpose(results))

    def collectParticleYield(self, particleName = 'pion_p'):
        """
            return event averaged particle yield (y or eta, dN/dy, dN/(deta))
            event loop over all the hydro + UrQMD events
        """
        eta_range = linspace(-3, 3, 61)
        pid = self.pid_lookup[particleName]
        if pid == 1:    # all charged hadrons
            pidString = ""
            for aPart in self.chargedHadronList:
                pidString += "pid = %d or " % self.pid_lookup[aPart]
            pidString += "pid = 1"
        else:
            pidString = "pid = %d" % pid
        etaavg = []
        dNdy = []
        dNdyerr = []
        dNdeta = []
        dNdetaerr = []
        Nev = self.getNumberOftotalEvents()
        for ieta in range(len(eta_range)-1):
            etalow = eta_range[ieta]
            etahigh = eta_range[ieta+1]
            deta = etahigh - etalow
            #fetch dN/dy data
            whereClause = "(%s)" % pidString
            whereClause += " and (%g<=rapidity and rapidity<=%g)" % (etalow, etahigh)
            tempdata = self.db.selectFromTable("particle_list", "count(*), avg(rapidity)", whereClause=whereClause)
            deltaN = tempdata[0][0]
            if tempdata[0][1] == None:
                etaavg.append((etalow + etahigh)/2.)
            else:
                etaavg.append(tempdata[0][1])
            dNdy.append(deltaN/deta/Nev)
            dNdyerr.append(sqrt(deltaN/deta/Nev)/sqrt(Nev))
            #fetch dN/deta data
            whereClause = "(%s)" % pidString
            whereClause += " and (%g<=pseudorapidity and pseudorapidity<=%g)" % (etalow, etahigh)
            tempdata = self.db.selectFromTable("particle_list", "count(*), avg(pseudorapidity)", whereClause=whereClause)
            deltaN = tempdata[0][0]
            dNdeta.append(deltaN/deta/Nev)
            dNdetaerr.append(sqrt(deltaN/deta/Nev)/sqrt(Nev))

        return(array([etaavg, dNdy, dNdyerr, dNdeta, dNdetaerr]).transpose())
    
    def collectBasicParticleYield(self):
        """
            collect particle yield into database for commonly interested hadrons
        """
        BasicParticleList = ['pion_p', 'kaon_p', 'proton']
        for aParticle in BasicParticleList:
            pid = self.pid_lookup[aParticle]
            dNdata = self.collectParticleYield(aParticle)
            for idx in range(len(dNdata[:,0])):
                self.db.insertIntoTable("particleYield", (pid, dNdata[idx,0], dNdata[idx,1], dNdata[idx,2], dNdata[idx,3], dNdata[idx,4]))
            self.db._dbCon.commit()  # commit changes

    def getParticleYieldvsrap(self, particleName = 'pion_p', rap_range = linspace(-2.5, 2.5, 51)):
        """
            check particle yield is pre-collected or not. If not, collect particle yield
            into the database. It returns event averaged particle yield with rapidity_range
            or pseudorap_range required by the users.
        """
        eps = 1e-15
        pid = self.pid_lookup[particleName]
        #check particle spectrum table is exist or not
        if self.db.createTableIfNotExists("particleYield", (("pid","integer"), ("eta","real"), ("dNdy", "real"), ("dNdyerr", "real"), ("dNdeta", "real"), ("dNdetaerr", "real") )):
            self.collectBasicParticleYield()
        dNdata = array(self.db._executeSQL("select eta, dNdy, dNdyerr, dNdeta, dNdetaerr from particleYield where pid = %d " % pid).fetchall())
        if(dNdata.size == 0):
            dNdata = self.collectParticleYield(particleName)
            if (sum(abs(dNdata[:,1])) <  1e-15):
                print "There is no record of particle: %s in the database" % particleName
                return None
            else:
                for idx in range(len(dNdata[:,0])):
                    self.db.insertIntoTable("particleYield", (pid, dNdata[idx,0], dNdata[idx,1], dNdata[idx,2], dNdata[idx,3], dNdata[idx,4]))
                self.db._dbCon.commit()  # commit changes

        #interpolate results to desired pT range
        dNdyinterp = interp(rap_range, dNdata[:,0], dNdata[:,1])
        dNdyinterp_err = interp(rap_range, dNdata[:,0], dNdata[:,2])
        dNdetainterp = interp(rap_range, dNdata[:,0], dNdata[:,3])
        dNdetainterp_err = interp(rap_range, dNdata[:,0], dNdata[:,4])
        results = array([rap_range, dNdyinterp, dNdyinterp_err, dNdetainterp, dNdetainterp_err])
        return(transpose(results))

    def getParticleYield(self, particleName = 'pion_p', rapRange = [-0.5, 0.5], pseudorapRange = [-0.5, 0.5]):
        """
            return the particle yield of particle species "particleName" within given 
            rapidity or pseudorapidity range by users
        """
        npoint = 50
        Nev = self.getNumberOftotalEvents()
        
        rapArray = linspace(rapRange[0], rapRange[1], npoint)
        dy = rapArray[1] - rapArray[0]
        dNdydata = self.getParticleYieldvsrap(particleName, rap_range = rapArray)
        dNdy = sum(dNdydata[:,1])*dy
        dNdyerr = sqrt(dNdy)/sqrt(Nev)
        
        pseudorapArray = linspace(pseudorapRange[0], pseudorapRange[1], npoint)
        deta = pseudorapArray[1] - pseudorapArray[0]
        dNdetadata = self.getParticleYieldvsrap(particleName, rap_range = pseudorapArray)
        dNdeta = sum(dNdetadata[:,3])*deta
        dNdetaerr = sqrt(dNdeta)/sqrt(Nev)

        return(dNdy, dNdyerr, dNdeta, dNdetaerr)


def printHelpMessageandQuit():
    print "Usage : "
    print "particleReader.py databaseName"
    exit(0)

if __name__ == "__main__":
    if len(argv) < 2:
        printHelpMessageandQuit()
    test = particleReader(str(argv[1]))
    print(test.getParticleSpectrum('charged', pT_range = linspace(0,3,31)))
    print(test.getParticleYieldvsrap('charged', rap_range = linspace(-2,2,41)))
    print(test.getParticleYield('charged'))

