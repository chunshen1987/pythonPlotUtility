#! /usr/bin/env python

from sys import argv, exit
from os import path, remove
from DBR import SqliteDB
from numpy import *
from random import shuffle
import scipy.special
import time

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

        # get number of events
        # pre-collect it to shorten the initial loading time of the database
        if self.db.createTableIfNotExists("number_of_events", (("Nev_tot", "integer"), ("Nev_hydro", "integer"))):
            self.totNev = self.getNumberOftotalEvents()
            self.hydroNev = self.getNumberOfHydroEvents()
            self.db.insertIntoTable("number_of_events", (int(self.totNev), int(self.hydroNev)))
            self.db._dbCon.commit()  # commit changes
        else:
            self.totNev = self.db.executeSQLquery("select Nev_tot from number_of_events").fetchall()[0][0]
            self.hydroNev = self.db.executeSQLquery("select Nev_hydro from number_of_events").fetchall()[0][0]

        if self.db.createTableIfNotExists("UrQMD_NevList", (("hydroEventId", "integer"), ("Number_of_UrQMDevents", "integer"))):
            for hydroEventId in range(1, self.hydroNev+1):
                UrQMDNev = self.getNumberOfUrQMDEvents(hydroEventId)
                self.db.insertIntoTable("UrQMD_NevList", (int(self.hydroNev), int(UrQMDNev)))
            self.db._dbCon.commit()  # commit changes

    ################################################################################
    # functions to get number of events
    ################################################################################ 
    def getNumberOfHydroEvents(self):
        """
            return total number hydro events stored in the database
        """
        Nev = self.db.executeSQLquery("select count(*) from (select distinct hydroEvent_id from particle_list)").fetchall()[0][0]
        return(Nev)

    def getNumberOfUrQMDEvents(self, hydroEventid):
        """
            return number of UrQMD events for the given hydro event
        """
        Nev = self.db.executeSQLquery("select count(*) from (select distinct UrQMDEvent_id from particle_list where hydroEvent_id = %d)" % hydroEventid).fetchall()[0][0]
        return(Nev)

    def getNumberOftotalEvents(self):
        """
            return total number of events stored in the database
        """
        hydroEventid = self.db.executeSQLquery("select distinct hydroEvent_id from particle_list").fetchall()
        Nev = 0
        for iev in range(len(hydroEventid)):
            Nev += self.getNumberOfUrQMDEvents(hydroEventid[iev][0])
        return(Nev)

    def getPidString(self, particleName = 'pion_p'):
        if isinstance(particleName, list):
            pidList = []
            for aPart in particleName:
                pidList.append(str(self.pid_lookup[aPart]))
            pidString = " or ".join(map(lambda(x): 'pid = ' + x, pidList))
        else:
            pid = self.pid_lookup[particleName]
            if pid == 1:    # all charged hadrons
                pidString = self.getPidString(self.chargedHadronList)
            else:
                pidString = "pid = %d" % pid
        return(pidString)
    
    ################################################################################
    # functions to collect particle spectra and yields
    ################################################################################ 
    def collectParticleSpectrum(self, particleName="pion_p", rapidity_range = [-0.5, 0.5], pseudorap_range = [-0.5, 0.5]):
        """
            return event averaged particle spectrum (pT, dN/(dydpT), dN/(detadpT))
            event loop over all the hydro + UrQMD events
        """
        print("collect particle spectra of %s ..." % particleName)
        pT_range = linspace(0, 3, 31)
        pidString = self.getPidString(particleName)
        pTavg = []
        dNdydpT = []
        dNdydpTerr = []
        pTetaavg = []
        dNdetadpT = []
        dNdetadpTerr = []
        Nev = self.totNev
        for ipT in range(len(pT_range)-1):
            pTlow = pT_range[ipT]
            pThigh = pT_range[ipT+1]
            dpT = pThigh - pTlow
            #fetch dN/dydpT data
            tempdata = self.db.executeSQLquery("select count(*), avg(pT) from particle_list where (%s) and (%g <= pT and pT < %g) and (%g <= rapidity and rapidity <= %g)" % (pidString, pTlow, pThigh, rapidity_range[0], rapidity_range[1])).fetchall()
            deltaN = tempdata[0][0]
            if tempdata[0][1] == None:
                pTavg.append((pTlow + pThigh)/2.)
            else:
                pTavg.append(tempdata[0][1])
            dNdydpT.append(deltaN/dpT/Nev)
            dNdydpTerr.append(sqrt(deltaN/dpT/Nev)/sqrt(Nev))
            #fetch dN/detadpT data
            tempdata = self.db.executeSQLquery("select count(*), avg(pT) from particle_list where (%s) and (%g <= pT and pT < %g) and (%g <= pseudorapidity and pseudorapidity <= %g)" % (pidString, pTlow, pThigh, pseudorap_range[0], pseudorap_range[1])).fetchall()
            if tempdata[0][1] == None:
                pTetaavg.append((pTlow + pThigh)/2.)
            else:
                pTetaavg.append(tempdata[0][1])
            deltaN = tempdata[0][0]
            dNdetadpT.append(deltaN/dpT/Nev)
            dNdetadpTerr.append(sqrt(deltaN/dpT/Nev)/sqrt(Nev))

        return(array([pTavg, dNdydpT, dNdydpTerr, pTetaavg, dNdetadpT, dNdetadpTerr]).transpose())
    
    def collectBasicParticleSpectra(self):
        """
            collect particle spectra into database for commonly interested hadrons
        """
        BasicParticleList = ['pion_p', 'kaon_p', 'proton']
        for aParticle in BasicParticleList:
            pid = self.pid_lookup[aParticle]
            dNdata = self.collectParticleSpectrum(aParticle)
            for idx in range(len(dNdata[:,0])):
                self.db.insertIntoTable("particleSpectra", (pid, dNdata[idx,0], dNdata[idx,1], dNdata[idx,2], dNdata[idx,3], dNdata[idx,4], dNdata[idx,5]))
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
        if self.db.createTableIfNotExists("particleSpectra", (("pid","integer"), ("pT","real"), ("dNdydpT", "real"), ("dNdydpTerr", "real"), ("pTeta", "real"), ("dNdetadpT", "real"), ("dNdetadpTerr", "real") )):
            self.collectBasicParticleSpectra()
        dNdata = array(self.db.executeSQLquery("select pT, dNdydpT, dNdydpTerr, pTeta, dNdetadpT, dNdetadpTerr from particleSpectra where pid = %d " % pid).fetchall())
        if(dNdata.size == 0):
            dNdata = self.collectParticleSpectrum(particleName)
            if (sum(abs(dNdata[:,1])) <  1e-15):
                print "There is no record of particle: %s in the database" % particleName
                return None
            else:
                for idx in range(len(dNdata[:,0])):
                    self.db.insertIntoTable("particleSpectra", (pid, dNdata[idx,0], dNdata[idx,1], dNdata[idx,2], dNdata[idx,3], dNdata[idx,4], dNdata[idx,5]))
                self.db._dbCon.commit()  # commit changes

        #interpolate results to desired pT range
        dNdyinterp = exp(interp(pT_range, dNdata[:,0], log(dNdata[:,1]+eps)))
        dNdyinterp_err = exp(interp(pT_range, dNdata[:,0], log(dNdata[:,2]+eps)))
        dNdetainterp = exp(interp(pT_range, dNdata[:,3], log(dNdata[:,4]+eps)))
        dNdetainterp_err = exp(interp(pT_range, dNdata[:,3], log(dNdata[:,5]+eps)))
        results = array([pT_range, dNdyinterp, dNdyinterp_err, dNdetainterp, dNdetainterp_err])
        return(transpose(results))

    def collectParticleYield(self, particleName = 'pion_p'):
        """
            return event averaged particle yield (y or eta, dN/dy, dN/(deta))
            event loop over all the hydro + UrQMD events
        """
        print("collect particle yield of %s" % particleName)
        eta_range = linspace(-3, 3, 61)
        pidString = self.getPidString(particleName)
        etaavg = []
        dNdy = []
        dNdyerr = []
        dNdeta = []
        dNdetaerr = []
        Nev = self.totNev
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
        dNdata = array(self.db.executeSQLquery("select eta, dNdy, dNdyerr, dNdeta, dNdetaerr from particleYield where pid = %d " % pid).fetchall())
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
        Nev = self.totNev
       
        dy = (rapRange[1] - rapRange[0])/npoint
        rapArray = linspace(rapRange[0], rapRange[1]-dy, npoint) + dy/2.
        dNdydata = self.getParticleYieldvsrap(particleName, rap_range = rapArray)
        dNdy = sum(dNdydata[:,1])*dy
        dNdyerr = sqrt(dNdy)/sqrt(Nev)
        
        deta = (pseudorapRange[1] - pseudorapRange[0])/npoint
        pseudorapArray = linspace(pseudorapRange[0], pseudorapRange[1] - deta, npoint) + deta/2.
        dNdetadata = self.getParticleYieldvsrap(particleName, rap_range = pseudorapArray)
        dNdeta = sum(dNdetadata[:,3])*deta
        dNdetaerr = sqrt(dNdeta)/sqrt(Nev)

        return(dNdy, dNdyerr, dNdeta, dNdetaerr)

    
    ################################################################################
    # functions to collect particle emission function
    ################################################################################ 
    def collectParticleYieldvsSpatialVariable(self, particleName = 'pion_p', SVtype = 'tau', SV_range = linspace(0, 12, 61), rapType = 'rapidity', rap_range = [-0.5, 0.5]):
        """
            return event averaged particle yield per spacial variable, tau, x, y, or eta_s. 
            Default output is (tau, dN/dtau)
            event loop over all the hydro + UrQMD events
        """
        print("collect particle yield as a function of %s for %s" % (SVtype, particleName))
        pidString = self.getPidString(particleName)
        SV_avg = []
        dNdSV = []
        dNdSVerr = []
        Nev = self.totNev
        for iSV in range(len(SV_range)-1):
            SVlow = SV_range[iSV]
            SVhigh = SV_range[iSV+1]
            dSV = SVhigh - SVlow
            #fetch dN/dSV data
            dNdata = self.db.executeSQLquery("select count(*), avg(%s) from particle_list where (%s) and (%g <= %s and %s < %g) and (%g <= %s and %s <= %g)" % (SVtype, pidString, SVlow, SVtype, SVtype, SVhigh, rap_range[0], rapType, rapType, rap_range[1])).fetchall()
            deltaN = dNdata[0][0]
            if dNdata[0][1] == None:
                SV_avg.append((SVlow + SVhigh)/2.)
            else:
                SV_avg.append(dNdata[0][1])
            dNdSV.append(deltaN/dSV/Nev)
            dNdSVerr.append(sqrt(deltaN/dSV/Nev)/sqrt(Nev))

        return(array([SV_avg, dNdSV, dNdSVerr]).transpose())
    
    def getParticleYieldvsSpatialVariable(self, particleName = 'pion_p', SVtype = 'tau', SV_range = linspace(0, 12, 61), rapType = 'rapidity', rap_range = [-0.5, 0.5]):
        """
            store and retrieve event averaged particle yield per spacial variable, 
            tau, x, y, or eta_s. 
            Default output is (tau, dN/dtau)
            event loop over all the hydro + UrQMD events
        """
        pidString = self.getPidString(particleName)
        pid = self.pid_lookup[particleName]
        Nev = self.totNev
        if self.db.createTableIfNotExists("particleEmission_d%s_%s" % (SVtype, rapType), (("pid","integer"), ("%s" % SVtype, "real"), ("dNdy", "real"), ("dNdyerr", "real") )):
            dNdata = self.collectParticleYieldvsSpatialVariable(particleName, SVtype, SV_range, rapType, rap_range)
            for i in range(len(dNdata[:,0])):
                self.db.insertIntoTable("particleEmission_d%s_%s" % (SVtype, rapType), (pid, dNdata[i,0], dNdata[i,1], dNdata[i,2]))
            self.db._dbCon.commit()  # commit changes
        else:
            dNdata = array(self.db.executeSQLquery("select %s, dNdy, dNdyerr from particleEmission_d%s_%s where pid = %d" % (SVtype, SVtype, rapType, pid)).fetchall())
            if dNdata.size == 0:
                dNdata = self.collectParticleYieldvsSpatialVariable(particleName, SVtype, SV_range, rapType, rap_range)
                if dNdata.size == 0:
                    print("can not find the particle %s" % particleName)
                    return(array([]))
                else:
                    for i in range(len(dNdata[:,0])):
                        self.db.insertIntoTable("particleEmission_d%s_%s" % (SVtype, rapType), (pid, dNdata[i,0], dNdata[i,1], dNdata[i,2]))
                    self.db._dbCon.commit()  # commit changes

        return(dNdata)

    

    ################################################################################
    # functions to collect particle anisotropic flows
    ################################################################################ 
    def collectAvgdiffvnflow(self, particleName = 'pion_p', psiR = 0., rapidity_range = [-0.5, 0.5], pseudorap_range = [-0.5, 0.5]):
        """
            collect <cos(n*(phi_i - psiR))>, which is averaged over for all particles for the given particle 
            species over all events
        """
        print("collect averaged diff vn flow of %s ..." % particleName)
        pT_range = linspace(0, 3, 31)
        pidString = self.getPidString(particleName)
        pid = self.pid_lookup[particleName]
        norder = 6
        Nev = self.totNev
        tableNamesList = ["Avgdiffvnflow", "Avgdiffvnetaflow"]
        rapTypes = ['rapidity', 'pseudorapidity']
        for itable in range(len(tableNamesList)):
            for ipT in range(len(pT_range)-1):
                pTlow = pT_range[ipT]
                pThigh = pT_range[ipT+1]
                dpT = pThigh - pTlow
                data = array(self.db.executeSQLquery("select pT, phi_p from particle_list where (%s) and (%g <= pT and pT <= %g) and (%g <= %s and %s <= %g)" % (pidString, pTlow, pThigh, rapidity_range[0], rapTypes[itable], rapTypes[itable], rapidity_range[1])).fetchall())
                for iorder in range(1, norder+1):
                    if len(data) == 0:
                        pTavg = (pTlow + pThigh)/2.
                        vn_real = 0.0; vn_real_err = 0.0
                        vn_imag = 0.0; vn_imag_err = 0.0
                    else:
                        pTavg = mean(data[:,0])
                        temparray = cos(iorder*(data[:,1] - psiR))
                        vn_real = mean(temparray)
                        vn_real_err = std(temparray)/sqrt(Nev)
                        temparray = sin(iorder*(data[:,1] - psiR))
                        vn_imag = mean(temparray)
                        vn_imag_err = std(temparray)/sqrt(Nev)
                        self.db.insertIntoTable(tableNamesList[itable], (pid, iorder, pTavg, vn_real, vn_real_err, vn_imag, vn_imag_err))
            self.db._dbCon.commit()  # commit changes

    def collectBasicParticleAvgdiffvnflow(self, psi_R = 0.):
        """
            collect particle averaged vn flow into database for commonly interested hadrons
        """
        BasicParticleList = ['pion_p', 'kaon_p', 'proton']
        for aParticle in BasicParticleList:
            self.collectAvgdiffvnflow(particleName = aParticle, psiR = psi_R)
    
    def getAvgdiffvnflow(self, particleName = "pion_p", psiR = 0., order = 2, pT_range = linspace(0.2, 2.5, 20), rapidity_range = [-0.5, 0.5], pseudorap_range = [-0.5, 0.5]):
        """
            retrieve and interpolate the particle average vn flow data with respect
            to event plane angle, psiR. 
            Collect the average vn flow into database if data does not exist
        """
        pid = self.pid_lookup[particleName]
        collectFlag = False
        if self.db.createTableIfNotExists("Avgdiffvnflow", (("pid", "integer"), ("n", "integer"), ("pT", "real"), ("vn_real", "real"), ("vn_real_err", "real"), ("vn_imag", "real"), ("vn_imag_err", "real") )):
            collectFlag = True
        if self.db.createTableIfNotExists("Avgdiffvnetaflow", (("pid", "integer"), ("n", "integer"), ("pT", "real"), ("vn_real", "real"), ("vn_real_err", "real"), ("vn_imag", "real"), ("vn_imag_err", "real") )):
            collectFlag = True
        if collectFlag:
            self.collectBasicParticleAvgdiffvnflow() 
        vndata = array(self.db.executeSQLquery("select pT, vn_real, vn_real_err, vn_imag, vn_imag_err from Avgdiffvnflow where pid = %d and n = %d" % (pid, order)).fetchall())
        vnetadata = array(self.db.executeSQLquery("select pT, vn_real, vn_real_err, vn_imag, vn_imag_err from Avgdiffvnetaflow where pid = %d and n = %d" % (pid, order)).fetchall())
        if(vndata.size == 0):
            self.collectAvgdiffvnflow(particleName)
            vndata = array(self.db.executeSQLquery("select pT, vn_real, vn_real_err, vn_imag, vn_imag_err from Avgdiffvnflow where pid = %d and n = %d" % (pid, order)).fetchall())
            vnetadata = array(self.db.executeSQLquery("select pT, vn_real, vn_real_err, vn_imag, vn_imag_err from Avgdiffvnetaflow where pid = %d and n = %d" % (pid, order)).fetchall())
        
        vnmag = vndata[:,1]*cos(order*psiR) + vndata[:,3]*sin(order*psiR)
        vnmag_err = sqrt((vndata[:,2]*cos(order*psiR))**2. + (vndata[:,4]*sin(order*psiR))**2.)
        vnetamag = vnetadata[:,1]*cos(order*psiR) + vnetadata[:,3]*sin(order*psiR)
        vnetamag_err = sqrt((vnetadata[:,2]*cos(order*psiR))**2. + (vnetadata[:,4]*sin(order*psiR))**2.)
        #interpolate results to desired pT range
        vnpTinterp = interp(pT_range, vndata[:,0], vnmag)
        vnpTinterp_err = interp(pT_range, vndata[:,0], vnmag_err)
        vnpTetainterp = interp(pT_range, vnetadata[:,0], vnetamag)
        vnpTetainterp_err = interp(pT_range, vnetadata[:,0], vnetamag_err)
        results = array([pT_range, vnpTinterp, vnpTinterp_err, vnpTetainterp, vnpTetainterp_err])
        return(transpose(results))

    def collectAvgintevnflow(self, particleName = 'pion_p', psiR = 0.):
        """
            collect pT integrated <cos(n*(phi_i - psiR))> and <cos(n*(phi_i - psiR))> 
            as a function of rapidity or pseudorapidity.
            The ensemble average is averaged over for all particles for the given 
            particle species over all events
        """
        print("collect averaged inte vn flow of %s ..." % particleName)
        pidString = self.getPidString(particleName)
        pid = self.pid_lookup[particleName]
        eta_range = linspace(-3, 3, 61)
        norder = 6
        Nev = self.totNev
        tableNamesList = ["Avgintevnflow", "Avgintevnetaflow"]
        rapTypes = ['rapidity', 'pseudorapidity']
        for itable in range(len(tableNamesList)):
            for ieta in range(len(eta_range)-1):
                etalow = eta_range[ieta]
                etahigh = eta_range[ieta+1]
                deta = etahigh - etalow
                data = array(self.db.executeSQLquery("select %s, phi_p from particle_list where (%s) and (%g <= %s and %s <= %g)" % (rapTypes[itable], pidString, etalow, rapTypes[itable], rapTypes[itable], etahigh)).fetchall())
                for iorder in range(1, norder+1):
                    if len(data) == 0:
                        eta_avg = (etalow + etahigh)/2.
                        vn_real = 0.0; vn_real_err = 0.0
                        vn_imag = 0.0; vn_imag_err = 0.0
                    else:
                        eta_avg = mean(data[:,0])
                        temparray = cos(iorder*(data[:,1] - psiR))
                        vn_real = mean(temparray)
                        vn_real_err = std(temparray)/sqrt(Nev)
                        temparray = sin(iorder*(data[:,1] - psiR))
                        vn_imag = mean(temparray)
                        vn_imag_err = std(temparray)/sqrt(Nev)
                        self.db.insertIntoTable(tableNamesList[itable], (pid, iorder, eta_avg, vn_real, vn_real_err, vn_imag, vn_imag_err))
            self.db._dbCon.commit()  # commit changes

    def collectBasicParticleAvgintevnflow(self, psi_R = 0.):
        """
            collect particle averaged pT integrated vn flow into database for 
            commonly interested hadrons
        """
        BasicParticleList = ['pion_p', 'kaon_p', 'proton']
        for aParticle in BasicParticleList:
            self.collectAvgintevnflow(particleName = aParticle, psiR = psi_R)
    
    def getAvgintevnflowvsrap(self, particleName = "pion_p", psiR = 0., order = 2, rap_range = linspace(-2.5, 2.5, 20)):
        """
            retrieve and interpolate the particle average pT integrated vn flow data
            as a function of rapidity and pseudorapidity. 
            Collect the average vn flow into database if data does not exist
        """
        pid = self.pid_lookup[particleName]
        collectFlag = False
        if self.db.createTableIfNotExists("Avgintevnflow", (("pid", "integer"), ("n", "integer"), ("eta", "real"), ("vn_real", "real"), ("vn_real_err", "real"), ("vn_imag", "real"), ("vn_imag_err", "real") )):
            collectFlag = True
        if self.db.createTableIfNotExists("Avgintevnetaflow", (("pid", "integer"), ("n", "integer"), ("eta", "real"), ("vn_real", "real"), ("vn_real_err", "real"), ("vn_imag", "real"), ("vn_imag_err", "real") )):
            collectFlag = True
        if collectFlag:
            self.collectBasicParticleAvgintevnflow() 
        vndata = array(self.db.executeSQLquery("select eta, vn_real, vn_real_err, vn_imag, vn_imag_err from Avgintevnflow where pid = %d and n = %d" % (pid, order)).fetchall())
        vnetadata = array(self.db.executeSQLquery("select eta, vn_real, vn_real_err, vn_imag, vn_imag_err from Avgintevnetaflow where pid = %d and n = %d" % (pid, order)).fetchall())
        if(vndata.size == 0):
            self.collectAvgintevnflow(particleName)
            vndata = array(self.db.executeSQLquery("select eta, vn_real, vn_real_err, vn_imag, vn_imag_err from Avgintevnflow where pid = %d and n = %d" % (pid, order)).fetchall())
            vnetadata = array(self.db.executeSQLquery("select eta, vn_real, vn_real_err, vn_imag, vn_imag_err from Avgintevnetaflow where pid = %d and n = %d" % (pid, order)).fetchall())
        
        vnmag = vndata[:,1]*cos(order*psiR) + vndata[:,3]*sin(order*psiR)
        vnmag_err = sqrt((vndata[:,2]*cos(order*psiR))**2. + (vndata[:,4]*sin(order*psiR))**2.)
        vnetamag = vnetadata[:,1]*cos(order*psiR) + vnetadata[:,3]*sin(order*psiR)
        vnetamag_err = sqrt((vnetadata[:,2]*cos(order*psiR))**2. + (vnetadata[:,4]*sin(order*psiR))**2.)
        #interpolate results to desired pT range
        vnpTinterp = interp(rap_range, vndata[:,0], vnmag)
        vnpTinterp_err = interp(rap_range, vndata[:,0], vnmag_err)
        vnpTetainterp = interp(rap_range, vnetadata[:,0], vnetamag)
        vnpTetainterp_err = interp(rap_range, vnetadata[:,0], vnetamag_err)
        results = array([rap_range, vnpTinterp, vnpTinterp_err, vnpTetainterp, vnpTetainterp_err])
        return(transpose(results))
    
    def getParticleintevn(self, particleName = 'pion_p', psiR = 0., order = 2, rapRange = [-0.5, 0.5], pseudorapRange = [-0.5, 0.5]):
        """
            return pT-integrated vn of particle species "particleName" within given 
            rapidity or pseudorapidity range by users
        """
        npoint = 50
        Nev = self.totNev
        
        rapArray = linspace(rapRange[0], rapRange[1], npoint)
        dy = rapArray[1] - rapArray[0]
        vndata = self.getAvgintevnflowvsrap(particleName, psiR = psiR, order = order, rap_range = rapArray)
        vn = sum(vndata[:,1])*dy
        vn_err = mean(vndata[:,2])   # need to be improved
        
        pseudorapArray = linspace(pseudorapRange[0], pseudorapRange[1], npoint)
        deta = pseudorapArray[1] - pseudorapArray[0]
        vnetadata = self.getAvgintevnflowvsrap(particleName, psiR = psiR, order = order, rap_range = pseudorapArray)
        vneta = sum(vnetadata[:,3])*deta
        vneta_err = mean(vnetadata[:,4])   # need to be improved

        return(vn, vn_err, vneta, vneta_err)

    ################################################################################
    # functions to collect particle event plane anisotropic flows with finite resolution 
    ################################################################################ 
    def collectGlobalQnvectorforeachEvent(self):
        """
            collect event plane Qn and sub-event plane QnA, QnB vectors for nth order 
            harmonic flow calculated using all charged particles in the event
            Qn = 1/Nparticle * sum_i exp[i*n*phi_i]
        """
        particleName = "charged"
        weightTypes = ['1', 'pT']
        norder = 6
        Nev = self.totNev
        pidString = self.getPidString(particleName)
        
        etaGap = 1
        etaBound1 = etaGap/2.; etaBound2 = - etaGap/2.
        if self.db.createTableIfNotExists("globalQnvector", (("hydroEvent_id","integer"), ("UrQMDEvent_id", "integer"), ("weightType", "text"), ("n", "integer"), ("Nparticle", "integer"), ("Qn", "real"), ("Qn_psi", "real"), ("Nparticle_A", "integer"), ("subQnA", "real"), ("subQnA_psi", "real"), ("Nparticle_B", "integer"), ("subQnB", "real"), ("subQnB_psi", "real") )):
            for hydroId in range(1, self.hydroNev+1):
                UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
                cachedNev = 1000; offset = 0
                for UrQMDId in range(1, UrQMDNev+1):
                    # cache small temporary database to speed up processing
                    if UrQMDId % cachedNev == 1:
                        if offset != 0: tempDB.closeConnection(discardChanges = True)
                        print("Cache temporary database, please wait ...")
                        tempDB = SqliteDB(':memory:')
                        tempDB.createTableIfNotExists("particle_list", (("UrQMDEvent_id", "integer"), ("pT", "real"), ("phi_p", "real"), ("pseudorapidity", "real")))
                        tempDB.insertIntoTable("particle_list", self.db.executeSQLquery("select UrQMDEvent_id, pT, phi_p, pseudorapidity from particle_list where hydroEvent_id = %d and (%d < UrQMDEvent_id and UrQMDEvent_id <= %d) and (%s)" % (hydroId, offset, offset+cachedNev, pidString)).fetchall())
                        offset += cachedNev

                    print("processing event: (%d, %d) " % (hydroId, UrQMDId))
                    particleList = array(tempDB.executeSQLquery("select pT, phi_p, pseudorapidity from particle_list where UrQMDEvent_id = %d" % (UrQMDId)).fetchall())

                    pT = particleList[:,0]
                    phi = particleList[:,1]
                    eta = particleList[:,2]
                    Nparticle = len(pT)
                    Nparticle_A = min(len(pT[eta > etaBound1]), len(pT[eta < etaBound2])); Nparticle_B = Nparticle_A
                    for weightType in weightTypes:
                        if weightType == '1':
                            weight = ones(len(pT))
                        elif weightType == 'pT':
                            weight = pT
                        for order in range(1,norder+1):
                            #calculate subQn vector
                            idx_A = eta > etaBound1
                            subQnA_X = sum(weight[idx_A]*cos(order*phi[idx_A]))
                            subQnA_Y = sum(weight[idx_A]*sin(order*phi[idx_A]))
                            subPsi_nA = arctan2(subQnA_Y, subQnA_X)/order
                            subQnA = sqrt(subQnA_X**2. + subQnA_Y**2.)/Nparticle_A
                            idx_B = eta < etaBound2
                            subQnB_X = sum(weight[idx_B]*cos(order*phi[idx_B]))
                            subQnB_Y = sum(weight[idx_B]*sin(order*phi[idx_B]))
                            subPsi_nB = arctan2(subQnB_Y, subQnB_X)/order
                            subQnB = sqrt(subQnB_X**2. + subQnB_Y**2.)/Nparticle_B
                        
                            #calculate Qn vector
                            Qn_X = sum(weight*cos(order*phi))
                            Qn_Y = sum(weight*sin(order*phi))
                            Psi_n = arctan2(Qn_Y, Qn_X)/order
                            Qn = sqrt(Qn_X**2. + Qn_Y**2.)/Nparticle

                            self.db.insertIntoTable("globalQnvector", (hydroId, UrQMDId, weightType, order, Nparticle, Qn, Psi_n, Nparticle_A, subQnA, subPsi_nA, Nparticle_B, subQnB, subPsi_nB))
            self.db._dbCon.commit()  # commit changes
        else:
            print("Qn vectors from all charged particles are already collected!")
            inputval = raw_input("Do you want to delete the existing one and collect again?")
            if inputval.lower() == 'y' or inputval.lower() == 'yes':
                self.db.executeSQLquery("drop table globalQnvector")
                self.collectGlobalQnvectorforeachEvent()

    def getFullplaneResolutionFactor(self, resolutionFactor_sub, Nfactor):
        """
            use binary search for numerical solution for R(chi_s) = resolutionFactor_sub
            where R(chi) is defined in Eq. (7) in arXiv:0904.2315v3 and calculate 
            resolutionFactor for the full event
        """
        # check
        if(resolutionFactor_sub > 1.0):
            print("error: resolutionFactor_sub = % g,  is larger than 1!" % resolutionFactor_sub)
            exit(-1)
        
        tol = 1e-8 # accuracy
        
        #search boundary
        left = 0.0; right = 2.0 # R(2.0) > 1.0
        mid = (right + left)*0.5
        dis = right - left
        while dis > tol:
            midval = self.resolution_Rfunction(mid)
            diff = resolutionFactor_sub - midval
            if abs(diff) < tol:
                chi_s = mid
                break
            elif diff > tol :
                left = mid
            else:
                right = mid
            dis = right - left
            mid = (right + left)*0.5
        chi_s = mid
        return(self.resolution_Rfunction(chi_s*sqrt(Nfactor)))
        
    def resolution_Rfunction(self, chi):
        """
            R(chi) for calculating resolution factor for the full event
        """
        chisq = chi*chi
        result = sqrt(pi)/2*exp(-chisq/2)*chi*(scipy.special.i0(chisq/2) + scipy.special.i1(chisq/2))
        return result

    def calculcateResolutionfactor(self, particleName, order_range = [1, 6], oversampling = 1, weightType = 'pT', pT_range = [0.0, 3.0], rapType = "rapidity", rap_range = [-4.0, 4.0]):
        """
            calculate nth order full resolution factor using given species of
            particles within pT and rapidity cuts from all events with option 
            of oversampling.
        """
        Nev = self.totNev
        pidString = self.getPidString(particleName)
        num_of_order  = order_range[1] - order_range[0] + 1
        resolutionFactor_sub = zeros(num_of_order)
        for hydroId in range(1, self.hydroNev+1):
            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
            cachedNev = 1000; offset = 0
            for UrQMDId in range(1, UrQMDNev+1):
                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
                # cache small temporary database to speed up processing
                if UrQMDId % cachedNev == 1:
                    if offset != 0: tempDB.closeConnection(discardChanges = True)
                    print("Cache temporary database, please wait ...")
                    tempDB = SqliteDB(':memory:')
                    tempDB.createTableIfNotExists("particle_list", (("UrQMDEvent_id", "integer"), ("pT", "real"), ("phi_p", "real")))
                    tempDB.insertIntoTable("particle_list", self.db.executeSQLquery("select UrQMDEvent_id, pT, phi_p from particle_list where hydroEvent_id = %d and (%d < UrQMDEvent_id and UrQMDEvent_id <= %d) and (%s) and (%g <= pT and pT <= %g) and (%g <= %s and %s <= %g)" % (hydroId, offset, offset+cachedNev, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
                    offset += cachedNev
                particleList = array(tempDB.executeSQLquery("select pT, phi_p from particle_list where UrQMDEvent_id = %d " % (UrQMDId)).fetchall())

                #calculate resolution factor
                pT = particleList[:,0]
                phi = particleList[:,1]
                if weightType == '1':
                    weight = ones(len(pT))
                elif weightType == 'pT':
                    weight = pT
                Nparticle = len(pT)
                for iorder in range(num_of_order):
                    order = order_range[0] + iorder
                    for isample in range(oversampling):
                        idx = range(Nparticle); shuffle(idx)
                        idx_A = idx[0:int(Nparticle/2)]
                        shuffle(idx)
                        idx_B = idx[0:int(Nparticle/2)]

                        subQnA_X = 0.0; subQnA_Y = 0.0
                        for i in idx_A:
                            subQnA_X += weight[i]*cos(order*phi[i])
                            subQnA_Y += weight[i]*sin(order*phi[i])
                        subPsi_nA = arctan2(subQnA_Y, subQnA_X)/order
                        subQnB_X = 0.0; subQnB_Y = 0.0
                        for i in idx_B:
                            subQnB_X += weight[i]*cos(order*phi[i])
                            subQnB_Y += weight[i]*sin(order*phi[i])
                        subPsi_nB = arctan2(subQnB_Y, subQnB_X)/order
                        resolutionFactor_sub[iorder] += cos(order*(subPsi_nA - subPsi_nB))
        resolutionFactor_sub = sqrt(resolutionFactor_sub/Nev/oversampling)
        resolutionFactor_full = []
        for iorder in range(num_of_order):
            resolutionFactor_full.append(self.getFullplaneResolutionFactor(resolutionFactor_sub[iorder], 2.0))
        return(resolutionFactor_full)
                    
    def collectGlobalResolutionFactor(self):
        """
            collect and store the full resolution factors calculated using all charged particles
            from all the events for order n = 1-6
            R_sub = sqrt(<QnA*conj(QnB)/(|QnA||QnB|)>_ev), R_SP_sub = sqrt(<QnA*conj(QnB)>_ev)
        """
        weightTypes = ['1', 'pT']
        norder = 6
        if self.db.createTableIfNotExists("resolutionFactorR", (("weightType", "text"), ("n", "integer"), ("R","real"), ("R_sub", "real"), ("R_SP_sub", "real"))):
            for iorder in range(1, norder+1):
                for weightType in weightTypes:
                    subQn_data = array(self.db.executeSQLquery("select subQnA_psi, subQnB_psi, subQnA, subQnB from globalQnvector where weightType = '%s' and n = %d" % (weightType, iorder)).fetchall())
                    resolutionFactor_sub = sqrt(mean(cos(iorder*(subQn_data[:,0] - subQn_data[:,1]))))
                    resolutionFactor_SP_sub = sqrt(mean(subQn_data[:,2]*subQn_data[:,3]*cos(iorder*(subQn_data[:,0] - subQn_data[:,1]))))
                    resolutionFactor_full = self.getFullplaneResolutionFactor(resolutionFactor_sub, 2.0)
                    self.db.insertIntoTable("resolutionFactorR", (weightType, iorder, resolutionFactor_full, resolutionFactor_sub, resolutionFactor_SP_sub))
            self.db._dbCon.commit()  # commit changes
        else:
            print("Resolution factors from all charged particles are already collected!")
            inputval = raw_input("Do you want to delete the existing one and collect again?")
            if inputval.lower() == 'y' or inputval.lower() == 'yes':
                self.db.executeSQLquery("drop table resolutionFactorR")
                self.collectGlobalResolutionFactor()

    def collectinteEventplaneflow(self, particleName = 'pion_p', pT_range = [0.0, 4.0], rap_range = [-0.5, 0.5], rapType = "rapidity"):
        """
            collect nth order event plane harmonic flow, vn{EP}
            event plane angle is determined by all charged particles in the whole event
        """
        Nev = self.totNev
        weightTypes = ['1', 'pT']
        norder = 6
        pidString = self.getPidString(particleName)
        vn_obs = zeros(norder); vn_obs_pTweight = zeros(norder)
        vn_obs_sq = zeros(norder); vn_obs_pTweight_sq = zeros(norder)
        for hydroId in range(1, self.hydroNev+1):
            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
            cachedNev = 1000; offset = 0
            for UrQMDId in range(1, UrQMDNev+1):
                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
                # cache small temporary database to speed up processing
                if UrQMDId % cachedNev == 1:
                    if offset != 0: tempDB.closeConnection(discardChanges = True)
                    print("Cache temporary database, please wait ...")
                    tempDB = SqliteDB(':memory:')
                    tempDB.createTableIfNotExists("particle_list", (("UrQMDEvent_id", "integer"), ("pT", "real"), ("phi_p", "real")))
                    tempDB.insertIntoTable("particle_list", self.db.executeSQLquery("select UrQMDEvent_id, pT, phi_p from particle_list where hydroEvent_id = %d and (%d < UrQMDEvent_id and UrQMDEvent_id <= %d) and (%s) and (%g <= pT and pT <= %g) and (%g <= %s and %s <= %g)" % (hydroId, offset, offset+cachedNev, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
                    offset += cachedNev
                particleList = array(tempDB.executeSQLquery("select pT, phi_p from particle_list where UrQMDEvent_id = %d " % (UrQMDId)).fetchall())

                pT = particleList[:,0]
                phi = particleList[:,1]
                Nparticle = len(pT)
                for weightType in weightTypes:
                    if weightType == '1':
                        weight = ones(len(pT))
                    elif weightType == 'pT':
                        weight = pT
                    for iorder in range(1, norder + 1):
                        QnVector = array(self.db.executeSQLquery("select Nparticle, Qn, Qn_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = '%s' and n = %d" % (hydroId, UrQMDId, weightType, iorder)).fetchall())
                        Qn_X = QnVector[0,0]*QnVector[0,1]*cos(iorder*QnVector[0,2])
                        Qn_Y = QnVector[0,0]*QnVector[0,1]*sin(iorder*QnVector[0,2])

                        #calculate event-plane vn
                        vntemp = 0.0
                        for ipart in range(Nparticle):
                            # subtract self correlation
                            Xn = Qn_X - weight[ipart]*cos(iorder*phi[ipart])
                            Yn = Qn_Y - weight[ipart]*sin(iorder*phi[ipart])
                            psi_n = arctan2(Yn, Xn)/iorder
                            vntemp += cos(iorder*(phi[ipart] - psi_n))
                        if weightType == '1':
                            vn_obs[iorder-1] += vntemp/Nparticle
                            vn_obs_sq[iorder-1] += (vntemp/Nparticle)**2
                        elif weightType == 'pT':
                            vn_obs_pTweight[iorder-1] += vntemp/Nparticle
                            vn_obs_pTweight_sq[iorder-1] += (vntemp/Nparticle)**2
        vn_obs /= Nev
        vn_obs_err = sqrt(vn_obs_sq/Nev - (vn_obs)**2.)/sqrt(Nev)
        vn_obs_pTweight /= Nev
        vn_obs_pTweight_err = sqrt(vn_obs_pTweight_sq/Nev - (vn_obs_pTweight)**2.)/sqrt(Nev)
        vnEP = zeros(norder); vnEP_pTweight = zeros(norder)
        vnEP_err = zeros(norder); vnEP_pTweight_err = zeros(norder)
        for iorder in range(1, norder + 1):
            resolutionFactor = self.db.executeSQLquery("select R from resolutionFactorR where weightType = '1' and n = %d" % (iorder)).fetchall()[0][0]
            resolutionFactor_pTweight = self.db.executeSQLquery("select R from resolutionFactorR where weightType = 'pT' and n = %d" % (iorder)).fetchall()[0][0]
            vnEP[iorder-1] = vn_obs[iorder-1]/resolutionFactor
            vnEP_err[iorder-1] = vn_obs_err[iorder-1]/resolutionFactor
            vnEP_pTweight[iorder-1] = vn_obs_pTweight[iorder-1]/resolutionFactor_pTweight
            vnEP_pTweight_err[iorder-1] = vn_obs_pTweight_err[iorder-1]/resolutionFactor_pTweight
        
        return(vnEP, vnEP_err, vnEP_pTweight, vnEP_pTweight_err)

    def collectdiffEventplaneflow(self, particleName = 'pion_p', pT_range = [0.0, 3.0], rap_range = [-0.5, 0.5], rapType = "rapidity"):
        """
            collect nth order pT differential event plane and subevent plane harmonic 
            flow, vn{EP}(pT) and vn{subEP}(pT). 
            event plane angle is determined by all charged particles in the whole event
            subevent plane angle is determined by all charged particles in one sub event
            with eta gap
        """
        print("Collecting differential event plane flow vn{EP} and vn{subEP} for %s ..." % particleName)
        Nev = self.totNev
        norder = 6; npT = 30
        weightTypes = ['1', 'pT']
        pidString = self.getPidString(particleName)

        pTBoundary = linspace(pT_range[0], pT_range[1], npT+1)
        pTmean = zeros([npT]); NevpT = zeros([npT])
        vn_obs = zeros([norder, npT]); vn_obs_pTweight = zeros([norder, npT])
        vn_obs_sq = zeros([norder, npT]); vn_obs_pTweight_sq = zeros([norder, npT])
        vn_obs_sub = zeros([norder, npT]); vn_obs_pTweight_sub = zeros([norder, npT])
        vn_obs_sq_sub = zeros([norder, npT]); vn_obs_pTweight_sq_sub = zeros([norder, npT])
        vn_obs_err = zeros([norder, npT]); vn_obs_pTweight_err = zeros([norder, npT])
        vn_obs_sub_err = zeros([norder, npT]); vn_obs_pTweight_sub_err = zeros([norder, npT])
        
        cachedNev = 1000; tempDBFlag = False
        for hydroId in range(1, self.hydroNev+1):
            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
            offset = 0
            for UrQMDId in range(1, UrQMDNev+1):
                # cache small temporary database to speed up processing
                if UrQMDId % cachedNev == 1:
                    if tempDBFlag: tempDB.closeConnection(discardChanges = True)
                    print("Cache temporary database, please wait ...")
                    tempDB = SqliteDB(':memory:'); tempDBFlag = True
                    tempDB.createTableIfNotExists("particle_list", (("UrQMDEvent_id", "integer"), ("pT", "real"), ("phi_p", "real")))
                    tempDB.insertIntoTable("particle_list", self.db.executeSQLquery("select UrQMDEvent_id, pT, phi_p from particle_list where hydroEvent_id = %d and (%d < UrQMDEvent_id and UrQMDEvent_id <= %d) and (%s) and (%g <= pT and pT < %g) and (%g <= %s and %s <= %g)" % (hydroId, offset, offset+cachedNev, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
                    offset += cachedNev
                
                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
                Qn_X = zeros([2, norder]); Qn_Y = zeros([2, norder])
                subQnA_psi = zeros([2, norder])
                for iorder in range(1, norder + 1):
                    QnVector = array(self.db.executeSQLquery("select Nparticle, Qn, Qn_psi, subQnA_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = '1' and n = %d" % (hydroId, UrQMDId, iorder)).fetchall())
                    Qn_X[0, iorder-1] = QnVector[0,0]*QnVector[0,1]*cos(iorder*QnVector[0,2])
                    Qn_Y[0, iorder-1] = QnVector[0,0]*QnVector[0,1]*sin(iorder*QnVector[0,2])
                    subQnA_psi[0, iorder-1] = QnVector[0,3]
                    QnVector = array(self.db.executeSQLquery("select Nparticle, Qn, Qn_psi, subQnA_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = 'pT' and n = %d" % (hydroId, UrQMDId, iorder)).fetchall())
                    Qn_X[1, iorder-1] = QnVector[0,0]*QnVector[0,1]*cos(iorder*QnVector[0,2])
                    Qn_Y[1, iorder-1] = QnVector[0,0]*QnVector[0,1]*sin(iorder*QnVector[0,2])
                    subQnA_psi[1, iorder-1] = QnVector[0,3]
                 
                for ipT in range(npT):
                    pTlow = pTBoundary[ipT]; pThigh = pTBoundary[ipT+1]
                    particleList = array(tempDB.executeSQLquery("select pT, phi_p from particle_list where UrQMDEvent_id = %d and (%g <= pT and pT < %g)" % (UrQMDId, pTlow, pThigh)).fetchall())
                    
                    if(particleList.size == 0):  # no particle in the event
                        pTmean[ipT] = (pTlow + pThigh)/2.
                    else:
                        NevpT[ipT] += 1
                        pT = particleList[:,0]
                        phi = particleList[:,1]
                        pTmean[ipT] = mean(pT)
                        Nparticle = len(pT)
                        for weightType in weightTypes:
                            if weightType == '1':
                                weight = ones(len(pT))
                                iweight = 0
                            elif weightType == 'pT':
                                weight = pT
                                iweight = 1
                            for iorder in range(1, norder + 1):
                                Qn_X_local = Qn_X[iweight, iorder-1]
                                Qn_Y_local = Qn_Y[iweight, iorder-1]
                                subPsi_n = subQnA_psi[iweight, iorder-1]

                                #calculate event-plane vn
                                vntemp = 0.0
                                for ipart in range(Nparticle):
                                    # subtract self correlation
                                    Xn = Qn_X_local - weight[ipart]*cos(iorder*phi[ipart])
                                    Yn = Qn_Y_local - weight[ipart]*sin(iorder*phi[ipart])
                                    psi_n = arctan2(Yn, Xn)/iorder
                                    vntemp += weight[ipart]*cos(iorder*(phi[ipart] - psi_n))
                                vntemp = vntemp/Nparticle
                                subvntemp = sum(weight*cos(iorder*(phi - subPsi_n)))/Nparticle
                                if weightType == '1':
                                    vn_obs[iorder-1][ipT] += vntemp
                                    vn_obs_sq[iorder-1][ipT] += vntemp**2
                                    vn_obs_sub[iorder-1][ipT] += subvntemp
                                    vn_obs_sq_sub[iorder-1][ipT] += subvntemp**2
                                elif weightType == 'pT':
                                    vn_obs_pTweight[iorder-1][ipT] += vntemp
                                    vn_obs_pTweight_sq[iorder-1][ipT] += vntemp**2
                                    vn_obs_pTweight_sub[iorder-1][ipT] += subvntemp
                                    vn_obs_pTweight_sq_sub[iorder-1][ipT] += subvntemp**2

        eps = 1e-15
        vn_obs = vn_obs/(NevpT + eps)
        vn_obs_err = sqrt(vn_obs_sq/(NevpT + eps) - (vn_obs)**2.)/sqrt(NevpT + eps)
        vn_obs_pTweight = vn_obs_pTweight/(NevpT + eps)
        vn_obs_pTweight_err = sqrt(vn_obs_pTweight_sq/(NevpT + eps) - (vn_obs_pTweight)**2.)/sqrt(NevpT + eps)
        vn_obs_sub = vn_obs_sub/(NevpT + eps)
        vn_obs_sub_err = sqrt(vn_obs_sq_sub/(NevpT + eps) - (vn_obs_sub)**2.)/sqrt(NevpT + eps)
        vn_obs_pTweight_sub = vn_obs_pTweight_sub/(NevpT + eps)
        vn_obs_pTweight_sub_err = sqrt(vn_obs_pTweight_sq_sub/(NevpT + eps) - (vn_obs_pTweight_sub)**2.)/sqrt(NevpT + eps)

        vnEP = zeros([norder, npT]); vnEP_pTweight = zeros([norder, npT])
        vnEP_err = zeros([norder, npT]); vnEP_pTweight_err = zeros([norder, npT])
        vnsubEP = zeros([norder, npT]); vnsubEP_pTweight = zeros([norder, npT])
        vnsubEP_err = zeros([norder, npT]); vnsubEP_pTweight_err = zeros([norder, npT])
        for iorder in range(1, norder + 1):
            resolutionFactor = self.db.executeSQLquery("select R from resolutionFactorR where weightType = '1' and n = %d" % (iorder)).fetchall()[0][0]
            resolutionFactor_pTweight = self.db.executeSQLquery("select R from resolutionFactorR where weightType = 'pT' and n = %d" % (iorder)).fetchall()[0][0]
            resolutionFactor_sub = self.db.executeSQLquery("select R_sub from resolutionFactorR where weightType = '1' and n = %d" % (iorder)).fetchall()[0][0]
            resolutionFactor_pTweight_sub = self.db.executeSQLquery("select R_sub from resolutionFactorR where weightType = 'pT' and n = %d" % (iorder)).fetchall()[0][0]
            for ipT in range(npT):
                vnEP[iorder-1, ipT] = vn_obs[iorder-1, ipT]/resolutionFactor
                vnEP_err[iorder-1, ipT] = vn_obs_err[iorder-1, ipT]/resolutionFactor
                vnEP_pTweight[iorder-1, ipT] = vn_obs_pTweight[iorder-1, ipT]/resolutionFactor_pTweight
                vnEP_pTweight_err[iorder-1, ipT] = vn_obs_pTweight_err[iorder-1, ipT]/resolutionFactor_pTweight
                vnsubEP[iorder-1, ipT] = vn_obs_sub[iorder-1, ipT]/resolutionFactor_sub
                vnsubEP_err[iorder-1, ipT] = vn_obs_sub_err[iorder-1, ipT]/resolutionFactor_sub
                vnsubEP_pTweight[iorder-1, ipT] = vn_obs_pTweight_sub[iorder-1, ipT]/resolutionFactor_pTweight_sub
                vnsubEP_pTweight_err[iorder-1, ipT] = vn_obs_pTweight_sub_err[iorder-1, ipT]/resolutionFactor_pTweight_sub
        
        return(pTmean, vnEP, vnEP_err, vnEP_pTweight, vnEP_pTweight_err, vnsubEP, vnsubEP_err, vnsubEP_pTweight, vnsubEP_pTweight_err)

    def getdiffEventplaneflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
        """
            retrieve nth order event plane and subevent plane flow data from database 
            for the given species of particles with given pT_range.
            if no results is found, it will collect vn{EP}(pT) and vn{subEP}(pT) and 
            store them into the database
        """
        pid = self.pid_lookup[particleName]
        collectFlag = False
        if rapType == 'rapidity':
            tableName = "diffvnEP"
        elif rapType == 'pseudorapidity':
            tableName = "diffvnEPeta"
        if self.db.createTableIfNotExists(tableName, (("pid", "integer"), ("weightType", "text"), ("n", "integer"), ("pT", "real"), ("vnEP", "real"), ("vnEP_err", "real"), ("vnsubEP", "real"), ("vnsubEP_err", "real") )):
            collectFlag = True
        else:
            vnEPdata = array(self.db.executeSQLquery("select pT, vnEP, vnEP_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName, pid, weightType, order)).fetchall())
            vnsubEPdata = array(self.db.executeSQLquery("select pT, vnsubEP, vnsubEP_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName, pid, weightType, order)).fetchall())
            if vnEPdata.size == 0: collectFlag = True
            if vnsubEPdata.size == 0: collectFlag = True
        if collectFlag:
            pT, vnEP, vnEP_err, vnEP_pTweight, vnEP_pTweight_err, vnsubEP, vnsubEP_err, vnsubEP_pTweight, vnsubEP_pTweight_err = self.collectdiffEventplaneflow(particleName = particleName, rapType = rapType) 
            for iorder in range(len(vnEP[:,0])):
                for ipT in range(len(pT)):
                    self.db.insertIntoTable(tableName, (pid, '1', iorder+1, pT[ipT], vnEP[iorder, ipT], vnEP_err[iorder, ipT], vnsubEP[iorder, ipT], vnsubEP_err[iorder, ipT]))
                    self.db.insertIntoTable(tableName, (pid, 'pT', iorder+1, pT[ipT], vnEP_pTweight[iorder, ipT], vnEP_pTweight_err[iorder, ipT], vnsubEP_pTweight[iorder, ipT], vnsubEP_pTweight_err[iorder, ipT]) )
            self.db._dbCon.commit()
            if weightType == '1':
                vnEPdata = array([pT, vnEP[order-1,:], vnEP_err[order-1,:]]).transpose()
                vnsubEPdata = array([pT, vnsubEP[order-1,:], vnsubEP_err[order-1,:]]).transpose()
            elif weightType == 'pT':
                vnEPdata = array([pT, vnEP_pTweight[order-1], vnEP_pTweight_err[order-1,:]]).transpose()
                vnsubEPdata = array([pT, vnsubEP_pTweight[order-1], vnsubEP_pTweight_err[order-1,:]]).transpose()
            if vnEPdata.size == 0 or vnsubEPdata.size == 0:
                print("There is no record for different event plane flow vn for %s" % particleName)
                return None
      
        #interpolate results to desired pT range
        vnEPpTinterp = interp(pT_range, vnEPdata[:,0], vnEPdata[:,1])
        vnEPpTinterp_err = interp(pT_range, vnEPdata[:,0], vnEPdata[:,2])
        vnsubEPpTinterp = interp(pT_range, vnsubEPdata[:,0], vnsubEPdata[:,1])
        vnsubEPpTinterp_err = interp(pT_range, vnsubEPdata[:,0], vnsubEPdata[:,2])
        results = array([pT_range, vnEPpTinterp, vnEPpTinterp_err, vnsubEPpTinterp, vnsubEPpTinterp_err])
        return(transpose(results))
                
    def collectdiffScalarProductflow(self, particleName = 'pion_p', pT_range = [0.0, 3.0], rap_range = [-0.5, 0.5], rapType = "rapidity"):
        """
            collect nth order pT differential scalar product harmonic flow, vn{SP}(pT)
            event plane angle is determined by all charged particles in the whole event
        """
        print("Collecting differential scalar product flow vn{SP} for %s ..." % particleName)
        Nev = self.totNev
        norder = 6; npT = 30
        weightTypes = ['1', 'pT']
        pidString = self.getPidString(particleName)

        pTBoundary = linspace(pT_range[0], pT_range[1], npT+1)
        pTmean = zeros(npT)
        vn_obs = zeros([norder, npT]); vn_obs_pTweight = zeros([norder, npT])
        vn_obs_sq = zeros([norder, npT]); vn_obs_pTweight_sq = zeros([norder, npT])

        cachedNev = 1000; tempDBFlag = False
        for hydroId in range(1, self.hydroNev+1):
            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
            offset = 0
            for UrQMDId in range(1, UrQMDNev+1):
                # cache small temporary database to speed up processing
                if UrQMDId % cachedNev == 1:
                    if tempDBFlag: tempDB.closeConnection(discardChanges = True)
                    print("Cache temporary database, please wait ...")
                    tempDB = SqliteDB(':memory:'); tempDBFlag = True
                    tempDB.createTableIfNotExists("particle_list", (("UrQMDEvent_id", "integer"), ("pT", "real"), ("phi_p", "real")))
                    tempDB.insertIntoTable("particle_list", self.db.executeSQLquery("select UrQMDEvent_id, pT, phi_p from particle_list where hydroEvent_id = %d and (%d < UrQMDEvent_id and UrQMDEvent_id <= %d) and (%s) and (%g <= pT and pT < %g) and (%g <= %s and %s <= %g)" % (hydroId, offset, offset+cachedNev, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
                    offset += cachedNev

                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
                QnA = zeros([2,norder]); Qn_psiA = zeros([2,norder])
                for iorder in range(1, norder + 1):
                    subQnVectors = array(self.db.executeSQLquery("select subQnA, subQnA_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = '1' and n = %d" % (hydroId, UrQMDId, iorder)).fetchall())
                    QnA[0,iorder-1] = subQnVectors[0,0]
                    Qn_psiA[0,iorder-1] = subQnVectors[0,1]
                    subQnVectors = array(self.db.executeSQLquery("select subQnA, subQnA_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = 'pT' and n = %d" % (hydroId, UrQMDId, iorder)).fetchall())
                    QnA[1,iorder-1] = subQnVectors[0,0]
                    Qn_psiA[1,iorder-1] = subQnVectors[0,1]
              
                for ipT in range(npT):
                    pTlow = pTBoundary[ipT]; pThigh = pTBoundary[ipT+1]
                    particleList = array(tempDB.executeSQLquery("select pT, phi_p from particle_list where UrQMDEvent_id = %d and (%g <= pT and pT <= %g)" % (UrQMDId, pTlow, pThigh)).fetchall())
                    
                    if(particleList.size == 0):  # no particle in the bin
                        pTmean[ipT] = (pTlow + pThigh)/2.
                    else:
                        pT = particleList[:,0]
                        phi = particleList[:,1]
                        pTmean[ipT] = mean(pT)
                        Nparticle = len(pT)

                        for weightType in weightTypes:
                            if weightType == '1':
                                weight = ones(len(pT))
                                iweight = 0
                            elif weightType == 'pT':
                                weight = pT
                                iweight = 1
                            for iorder in range(1, norder + 1):
                                QnA_local = QnA[iweight, iorder-1]
                                Qn_psiA_local = Qn_psiA[iweight, iorder-1]

                                # Qn vector for the interested particles
                                Qn_X = 0.0; Qn_Y = 0.0
                                for ipart in range(Nparticle):
                                    Qn_X += weight[ipart]*cos(iorder*phi[ipart])
                                    Qn_Y += weight[ipart]*sin(iorder*phi[ipart])
                                Qn = sqrt(Qn_X**2. + Qn_Y**2.)/Nparticle
                                Qn_psi = arctan2(Qn_Y, Qn_X)/iorder

                                temp = Qn*QnA_local*cos(iorder*(Qn_psi - Qn_psiA_local))

                                if weightType == '1':
                                    vn_obs[iorder-1][ipT] += temp
                                    vn_obs_sq[iorder-1][ipT] += temp*temp
                                elif weightType == 'pT':
                                    vn_obs_pTweight[iorder-1][ipT] += temp
                                    vn_obs_pTweight_sq[iorder-1][ipT] += temp*temp
        vn_obs = vn_obs/Nev
        vn_obs_err = sqrt(vn_obs_sq/Nev - (vn_obs)**2.)/sqrt(Nev)
        vn_obs_pTweight = vn_obs_pTweight/Nev
        vn_obs_pTweight_err = sqrt(vn_obs_pTweight_sq/Nev - (vn_obs_pTweight)**2.)/sqrt(Nev)

        vnSP = zeros([norder, npT]); vnSP_pTweight = zeros([norder, npT])
        vnSP_err = zeros([norder, npT]); vnSP_pTweight_err = zeros([norder, npT])

        for iorder in range(1, norder + 1):
            resolutionFactor = self.db.executeSQLquery("select R_SP_sub from resolutionFactorR where weightType = '1' and n = %d" % (iorder)).fetchall()[0][0]
            resolutionFactor_pTweight = self.db.executeSQLquery("select R_SP_sub from resolutionFactorR where weightType = 'pT' and n = %d" % (iorder)).fetchall()[0][0]
            for ipT in range(npT):
                vnSP[iorder-1, ipT] = vn_obs[iorder-1, ipT]/resolutionFactor
                vnSP_err[iorder-1, ipT] = vn_obs_err[iorder-1, ipT]/resolutionFactor
                vnSP_pTweight[iorder-1, ipT] = vn_obs_pTweight[iorder-1, ipT]/resolutionFactor_pTweight
                vnSP_pTweight_err[iorder-1, ipT] = vn_obs_pTweight_err[iorder-1, ipT]/resolutionFactor_pTweight
        
        return(pTmean, vnSP, vnSP_err, vnSP_pTweight, vnSP_pTweight_err)

    def getdiffScalarProductflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
        """
            retrieve nth order scalar product flow data from database for the given species 
            of particles with given pT_range.
            if no results is found, it will collect vn{SP}(pT) and store it into database
        """
        pid = self.pid_lookup[particleName]
        collectFlag = False
        if rapType == 'rapidity':
            tableName = "diffvnSP"
        elif rapType == 'pseudorapidity':
            tableName = "diffvnSPeta"
        if self.db.createTableIfNotExists(tableName, (("pid", "integer"), ("weightType", "text"), ("n", "integer"), ("pT", "real"), ("vn", "real"), ("vn_err", "real") )):
            collectFlag = True
        else:
            vnSPdata = array(self.db.executeSQLquery("select pT, vn, vn_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName, pid, weightType, order)).fetchall())
            if vnSPdata.size == 0:
                collectFlag = True
        if collectFlag:
            pT, vnSP, vnSP_err, vnSP_pTweight, vnSP_pTweight_err = self.collectdiffScalarProductflow(particleName = particleName, rapType = rapType) 
            for iorder in range(len(vnSP[:,0])):
                for ipT in range(len(pT)):
                    self.db.insertIntoTable(tableName, (pid, '1', iorder+1, pT[ipT], vnSP[iorder, ipT], vnSP_err[iorder, ipT]))
                    self.db.insertIntoTable(tableName, (pid, 'pT', iorder+1, pT[ipT], vnSP_pTweight[iorder, ipT], vnSP_pTweight_err[iorder, ipT]))
            self.db._dbCon.commit()
            if weightType == '1':
                vnSPdata = array([pT, vnSP[order-1,:], vnSP_err[order-1,:]]).transpose()
            elif weightType == 'pT':
                vnSPdata = array([pT, vnSP_pTweight[order-1], vnSP_pTweight_err[order-1,:]]).transpose()
            if vnSPdata.size == 0:
                print("There is no record for different event plane flow vn for %s" % particleName)
                return None
      
        #interpolate results to desired pT range
        vnpTinterp = interp(pT_range, vnSPdata[:,0], vnSPdata[:,1])
        vnpTinterp_err = interp(pT_range, vnSPdata[:,0], vnSPdata[:,2])
        results = array([pT_range, vnpTinterp, vnpTinterp_err])
        return(transpose(results))
    
    
    ################################################################################
    # functions to collect two particle correlation
    ################################################################################ 
    def collectTwoparticleCumulantflow(self, particleName = 'pion_p', pT_range = [0.0, 3.0], rap_range = [-0.5, 0.5], rapType = "rapidity"):
        """
            collect nth order pT differential two particle cumulant harmonic flow, vn{2}(pT)
            the two particles are taken from the same bin
        """
        print("Collecting differential two particle cumulant flow vn{2} for %s ..." % particleName)
        norder = 6; npT = 30
        weightTypes = ['1', 'pT']
        pidString = self.getPidString(particleName)

        pTBoundary = linspace(pT_range[0], pT_range[1], npT+1)
        pTmean = zeros([npT]); NevpT = zeros(npT)
        vn_obs = zeros([norder, npT])
        vn_obs_sq = zeros([norder, npT])
        vn_obs_err = zeros([norder, npT])

        for hydroId in range(1, self.hydroNev+1):
            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
            cachedNev = 1000; offset = 0
            for UrQMDId in range(1, UrQMDNev+1):
                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
                # cache small temporary database to speed up processing
                if UrQMDId % cachedNev == 1:
                    if offset != 0: tempDB.closeConnection(discardChanges = True)
                    print("Cache temporary database, please wait ...")
                    tempDB = SqliteDB(':memory:')
                    tempDB.createTableIfNotExists("particle_list", (("UrQMDEvent_id", "integer"), ("pT", "real"), ("phi_p", "real")))
                    tempDB.insertIntoTable("particle_list", self.db.executeSQLquery("select UrQMDEvent_id, pT, phi_p from particle_list where hydroEvent_id = %d and (%d < UrQMDEvent_id and UrQMDEvent_id <= %d) and (%s) and (%g <= pT and pT <= %g) and (%g <= %s and %s <= %g)" % (hydroId, offset, offset+cachedNev, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
                    offset += cachedNev

                for ipT in range(npT):
                    pTlow = pTBoundary[ipT]; pThigh = pTBoundary[ipT+1]
                    particleList = array(tempDB.executeSQLquery("select pT, phi_p from particle_list where UrQMDEvent_id = %d and (%g <= pT and pT < %g)" % (UrQMDId, pTlow, pThigh)).fetchall())
                    
                    if(particleList.size == 0):
                        pTmean[ipT] = (pTlow + pThigh)/2.
                    elif(len(particleList[:,0]) == 1):
                        pTmean[ipT] = particleList[0,0]
                    else:
                        NevpT[ipT] += 1
                        pT = particleList[:,0]
                        pTmean[ipT] = mean(pT)
                        phi = particleList[:,1]
                        Nparticle = len(pT)
                        Npairs = Nparticle*(Nparticle - 1)/2.

                        for iorder in range(1, norder + 1):
                            temp = 0.0
                            for ipart in range(Nparticle-1):
                                for iassopart in range(ipart+1, Nparticle):
                                    temp += cos(iorder*(phi[ipart] - phi[iassopart]))

                            vn_obs[iorder-1, ipT] += temp/Npairs
                            vn_obs_sq[iorder-1, ipT] += temp*temp/Npairs
        eps = 1e-15
        vn_obs = vn_obs/(NevpT + eps)
        vn_obs_err = sqrt(vn_obs_sq/(NevpT + eps) - (vn_obs)**2.)/sqrt(NevpT + eps)
        
        vn_2 = zeros([norder, npT])
        vn_2_err = zeros([norder, npT])
        for iorder in range(norder):
            for ipT in range(npT):
                if abs(vn_obs[iorder, ipT]) > eps:
                    vn_2[iorder, ipT] = sqrt(vn_obs[iorder, ipT])
                    vn_2_err[iorder, ipT] = vn_obs_err[iorder, ipT]/(2.*sqrt(vn_obs[iorder, ipT]))
        return(pTmean, vn_2, vn_2_err)

    def getdiffTwoparticlecumulantflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
        """
            retrieve nth order two particle cumulant flow data from database for the 
            given species of particles with given pT_range.
            if no results is found, it will collect vn{2}(pT) and store it into database
        """
        pid = self.pid_lookup[particleName]
        collectFlag = False
        if rapType == 'rapidity':
            tableName = "diffvn2"
        elif rapType == 'pseudorapidity':
            tableName = "diffvn2eta"
        if self.db.createTableIfNotExists(tableName, (("pid", "integer"), ("n", "integer"), ("pT", "real"), ("vn", "real"), ("vn_err", "real") )):
            collectFlag = True
        else:
            vn2data = array(self.db.executeSQLquery("select pT, vn, vn_err from %s where pid = %d and n = %d" % (tableName, pid, order)).fetchall())
            if vn2data.size == 0:
                collectFlag = True
        if collectFlag:
            pT, vn2, vn2_err = self.collectTwoparticleCumulantflow(particleName = particleName, rapType = rapType) 
            for iorder in range(len(vn2[:,0])):
                for ipT in range(len(pT)):
                    self.db.insertIntoTable(tableName, (pid, iorder+1, pT[ipT], vn2[iorder, ipT], vn2_err[iorder, ipT]))
            self.db._dbCon.commit()
            vn2data = array([pT, vn2[order-1,:], vn2_err[order-1,:]]).transpose()
            if vn2data.size == 0:
                print("There is no record for different event plane flow vn for %s" % particleName)
                return None
      
        #interpolate results to desired pT range
        vnpTinterp = interp(pT_range, vn2data[:,0], vn2data[:,1])
        vnpTinterp_err = interp(pT_range, vn2data[:,0], vn2data[:,2])
        results = array([pT_range, vnpTinterp, vnpTinterp_err])
        return(transpose(results))

    def collectTwoparticleCorrelation(self, particleName = 'pion_p', pT_range = [1.0, 2.0]):
        """
            collect two particle correlation function C(\delta phi, \delta eta) from 
            all the events within given pT range for a given particle species, "particleName"
        """
        Nev = self.totNev
        NevMix = 10
        nphi = 15; neta = 20
        dphi_bound = linspace(0, pi, nphi+1)
        deta_bound = linspace(0, 2, neta+1)
        CorrNum = zeros([nphi, neta])
        CorrDenorm = zeros([nphi, neta])
        pidString = self.getPidString(particleName)
        for hydroId in range(1, self.hydroNev+1):
            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
            cachedNev = 1000; offset = 0
            for UrQMDId in range(1, UrQMDNev+1):
                # cache small temporary database to speed up processing
                if UrQMDId % cachedNev == 1:
                    if offset != 0: tempDB.closeConnection(discardChanges = True)
                    print("Cache temporary database, please wait ...")
                    tempDB = SqliteDB(':memory:')
                    tempDB.createTableIfNotExists("particle_list", (("UrQMDEvent_id", "integer"), ("phi_p", "real"), ("pseudorapidity", "real")))
                    tempDB.insertIntoTable("particle_list", self.db.executeSQLquery("select UrQMDEvent_id, phi_p, pseudorapidity from particle_list where hydroEvent_id = %d and (%d < UrQMDEvent_id and UrQMDEvent_id <= %d) and (%s) and (%g <= pT and pT <= %g) " % (hydroId, offset, offset+cachedNev, pidString, pT_range[0], pT_range[1])).fetchall())
                    offset += cachedNev
                #processing pairs within same event for the numerator
                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
                if tempDB.doesTableExist("temp_particleList"):
                    tempDB.executeSQLquery("drop table temp_particleList")
                tempDB.createTableIfNotExists("temp_particleList", (("phi_p", "real"), ("pseudorapidity", "real")))  # create a temporary table
                particleList = tempDB.executeSQLquery("select phi_p, pseudorapidity from particle_list where UrQMDEvent_id = %d" % (UrQMDId)).fetchall()
                tempDB.insertIntoTable("temp_particleList", particleList)
                for iphi in range(nphi):
                    for ieta in range(neta):
                        Ntemp = 0
                        for ipart in range(len(particleList)):
                            phi_low = particleList[ipart][0] + dphi_bound[iphi]
                            phi_high = particleList[ipart][0] + dphi_bound[iphi+1]
                            eta_low = particleList[ipart][1] + deta_bound[ieta]
                            eta_high = particleList[ipart][1] + deta_bound[ieta+1]
                            Ntemp += tempDB.executeSQLquery("select count(*) from temp_particleList where (%g < phi_p and phi_p <= %g) and (%g < pseudorapidity and pseudorapidity <= %g)" % (phi_low, phi_high, eta_low, eta_high)).fetchall()[0][0] 
                        CorrNum[iphi, ieta] += Ntemp
                #processing pairs from different event for the denominator
                for ievMix in range(NevMix):
                    hydroMixevId = random.randint(self.hydroNev) + 1
                    #while (hydroMixevId == hydroId[0]): # require mixed event to be different a hydro event
                    #    hydroMixevId = random.randint(len(hydroIdList)) + 1
                    UrQMDMixevId = random.randint(offset-cachedNev, min(offset, UrQMDNev)) + 1
                    print("processing mixed events %d: event (%d, %d) vs event (%d, %d)" % (ievMix+1, hydroId, UrQMDId, hydroMixevId, UrQMDMixevId))
                    if tempDB.doesTableExist("temp_particleList_mixed"):
                        tempDB.executeSQLquery("drop table temp_particleList_mixed")
                    tempDB.createTableIfNotExists("temp_particleList_mixed", (("phi_p", "real"), ("pseudorapidity", "real")))  #create a temporary table
                    particleList_mixed = tempDB.executeSQLquery("select phi_p, pseudorapidity from particle_list where UrQMDEvent_id = %d" % (UrQMDMixevId)).fetchall()
                    tempDB.insertIntoTable("temp_particleList_mixed", particleList_mixed)
                    for iphi in range(nphi):
                        for ieta in range(neta):
                            Ntemp = 0
                            for ipart in range(len(particleList)):
                                phi_low = particleList[ipart][0] + dphi_bound[iphi]
                                phi_high = particleList[ipart][0] + dphi_bound[iphi+1]
                                eta_low = particleList[ipart][1] + deta_bound[ieta]
                                eta_high = particleList[ipart][1] + deta_bound[ieta+1]
                                Ntemp += tempDB.executeSQLquery("select count(*) from temp_particleList_mixed where (%g < phi_p and phi_p <= %g) and (%g < pseudorapidity and pseudorapidity <= %g)" % (phi_low, phi_high, eta_low, eta_high)).fetchall()[0][0] 
                            CorrDenorm[iphi, ieta] += Ntemp
                    tempDB.executeSQLquery("drop table temp_particleList_mixed")  # delete temporary table
                tempDB.executeSQLquery("drop table temp_particleList")  # delete temporary table

        CorrMatrix = CorrNum/CorrDenorm
        print(CorrMatrix)

def printHelpMessageandQuit():
    print "Usage : "
    print "particleReader.py databaseName"
    exit(0)

if __name__ == "__main__":
    if len(argv) < 2:
        printHelpMessageandQuit()
    test = particleReader(str(argv[1]))
    #for aPart in ['pion_p', 'kaon_p', 'proton']:
    #    print(test.getParticleYieldvsSpatialVariable(particleName = aPart, SVtype = 'tau', SV_range = linspace(0.0, 15.0, 76)))
    #    print(test.getParticleYieldvsSpatialVariable(particleName = aPart, SVtype = 'x', SV_range = linspace(-13.0, 13.0, 101)))
    #    print(test.getParticleYieldvsSpatialVariable(particleName = aPart, SVtype = 'eta', SV_range = linspace(-5.0, 5.0, 51)))
    #print(test.getAvgdiffvnflow(particleName = "charged", psiR = 0., order = 2, pT_range = linspace(0.0, 2.0, 20)))
    #print(test.getAvgintevnflowvsrap(particleName = "charged", psiR = 0., order = 2, rap_range = linspace(-2.0, 2.0, 20)))
    #print(test.getParticleintevn('charged'))
    #print(test.getParticleSpectrum('charged', pT_range = linspace(0,3,31)))
    #print(test.getParticleYieldvsrap('charged', rap_range = linspace(-2,2,41)))
    #print(test.getParticleYield('charged'))
    #test.collectGlobalQnvectorforeachEvent()
    #test.collectGlobalResolutionFactor()
    #print(test.getdiffEventplaneflow(particleName = 'pion_p'))
    #print(test.getdiffEventplaneflow(particleName = 'kaon_p'))
    #print(test.getdiffEventplaneflow(particleName = 'proton'))
    #print(test.getdiffScalarProductflow(particleName = 'pion_p'))
    #print(test.getdiffScalarProductflow(particleName = 'kaon_p'))
    #print(test.getdiffScalarProductflow(particleName = 'proton'))
    print(test.getdiffTwoparticlecumulantflow(particleName = 'pion_p'))
    print(test.getdiffTwoparticlecumulantflow(particleName = 'kaon_p'))
    print(test.getdiffTwoparticlecumulantflow(particleName = 'proton'))
    #test.collectEventplaneflow('charged', 2)
    #test.collectTwoparticleCorrelation()

