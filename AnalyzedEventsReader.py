#! /usr/bin/env python

from sys import argv, exit
from os import path
from DBR import SqliteDB
from numpy import *
from random import shuffle
import scipy.special

#define colors
purple = "\033[95m"
green = "\033[92m"
red = "\033[91m"
normal = "\033[0m"

class AnalyzedDataReader(object):
    """
        This class is used to perform event average for all the results
        that are collected from ParticleReader to get final experimental
        observables
    """
    def __init__(self, database):
        """
            Register databases
        """
        # setup database for analysis
        if isinstance(database, str):
            if path.exists(database):
                database = SqliteDB(database)
            else:
                raise ValueError(
                    "AnalyzedDataReader.__init__: input database %s can not "
                    "be found!" % (database))
        if isinstance(database, SqliteDB):
            self.db = database
        else:
            raise TypeError(
                "AnalyzedDataReader.__init__: the input database must be "
                "a string or a SqliteDB database.")
        
        # setup lookup tables
        self.pid_lookup = dict(self.db.selectFromTable("pid_lookup"))

        # define all charged hadrons
        self.charged_hadron_list = [
            "pion_p", "pion_m", "kaon_p", "kaon_m", "proton", "anti_proton",
            "sigma_p", "sigma_m", "anti_sigma_p", "anti_sigma_m",
            "xi_m", "anti_xi_m"]

        # create index
        if not self.db.doesIndexExist('index_particle_pT_spectra'):
            print("Create index for particle_pT_spectra table ...")
            self.db.executeSQLquery(
                "CREATE INDEX index_particle_pT_spectra ON "
                "particle_pT_spectra(hydro_event_id, urqmd_event_id, pid)")
        for aTable in ['flow_Qn_vectors', 'flow_Qn_vectors_pTdiff']:
            if not self.db.doesIndexExist('index_%s' % aTable):
                print("Create index for table %s ..." % aTable)
                self.db.executeSQLquery( "CREATE INDEX index_%s ON "
                    "%s(hydro_event_id, urqmd_event_id, pid, weight_type, n)"
                    % (aTable, aTable))

        # get number of events
        # pre-collect it to shorten the initial loading time of the database
        if self.db.createTableIfNotExists(
                "number_of_events", (("Nev_tot", "integer"),
                                     ("Nev_hydro", "integer"))):
            self.tot_nev = self.get_number_of_total_events()
            self.hydro_nev = self.get_number_of_hydro_events()
            self.db.insertIntoTable("number_of_events",
                                    (int(self.tot_nev), int(self.hydro_nev)))
            self.db._dbCon.commit()  # commit changes
        else:
            self.tot_nev = self.db.executeSQLquery(
                "select Nev_tot from number_of_events").fetchall()[0][0]
            self.hydro_nev = self.db.executeSQLquery(
                "select Nev_hydro from number_of_events").fetchall()[0][0]

        temp_flag = self.db.createTableIfNotExists(
            "UrQMD_NevList", (("hydroEventId", "integer"),
                              ("Number_of_UrQMDevents", "integer")))
        urqmd_nev_array = zeros(self.hydro_nev)
        for hydroEventId in range(1, self.hydro_nev+1):
            urqmd_nev = self.get_number_of_urqmd_events(hydroEventId)
            urqmd_nev_array[hydroEventId-1] = urqmd_nev
            if temp_flag:
                self.db.insertIntoTable("UrQMD_NevList",
                                        (int(hydroEventId), int(urqmd_nev)))
        self.db._dbCon.commit()  # commit changes

        self.process_nev = 10000
        #determine retrieve event boundaries
        self.nev_bin = int(self.tot_nev/self.process_nev) + 2
        if self.tot_nev % self.process_nev == 0: self.nev_bin -= 1
        self.event_bound_hydro = ones(self.nev_bin)
        self.event_bound_urqmd = ones(self.nev_bin)
        ihydro_ev = 1
        ibin = 1
        temp = urqmd_nev_array[0]
        while ihydro_ev <= self.hydro_nev:
            if temp > ibin*self.process_nev:
                self.event_bound_hydro[ibin] = ihydro_ev
                ibin += 1
            else:
                ihydro_ev += 1
                if ihydro_ev < self.hydro_nev:
                    temp += urqmd_nev_array[ihydro_ev-1]
                else:
                    self.event_bound_hydro[ibin] = ihydro_ev
        for ibin in range(1,self.nev_bin):
            self.event_bound_urqmd[ibin] = (ibin*self.process_nev 
                - sum(urqmd_nev_array[0:self.event_bound_hydro[ibin]-1]))

    
    ###########################################################################
    # functions to get number of events
    ########################################################################### 
    def get_number_of_hydro_events(self):
        """
            return total number hydro events stored in the database
        """
        nev = self.db.executeSQLquery(
            "select count(*) from (select distinct hydro_event_id "
            "from particle_pT_spectra)").fetchall()[0][0]
        return nev

    def get_number_of_urqmd_events(self, hydroeventid):
        """
            return number of UrQMD events for the given hydro event
        """
        nev = self.db.executeSQLquery(
            "select count(*) from (select distinct urqmd_event_id from "
            "particle_pT_spectra where hydro_event_id = %d)" % hydroeventid
        ).fetchall()[0][0]
        return nev

    def get_number_of_total_events(self):
        """
            return total number of events stored in the database
        """
        hydro_event_id = self.db.executeSQLquery(
            "select distinct hydro_event_id from particle_pT_spectra"
        ).fetchall()
        nev = 0
        for iev in range(len(hydro_event_id)):
            nev += self.get_number_of_urqmd_events(hydro_event_id[iev][0])
        return nev

    def get_pid_string(self, particleName):
        if isinstance(particleName, list):
            pidList = []
            for aPart in particleName:
                pidList.append(str(self.pid_lookup[aPart]))
            pidstring = " or ".join(map(lambda(x): 'pid = ' + x, pidList))
        else:
            pid = self.pid_lookup[particleName]
            if pid == 1:    # all charged hadrons
                pidstring = self.get_pid_string(self.charged_hadron_list)
            else:
                pidstring = "pid = %d" % pid
        return(pidstring)

    ###########################################################################
    # functions to collect particle spectra and yields
    ########################################################################### 
    def get_particle_spectra(self, particle_name, pT_range = linspace(0,3,31), 
                             rap_type = 'rapidity'):
        """
            This function performs event average for particle pT-spectra
            from the database and interpolates to the desired pT values
            specified by the users.
            It returns (pT, dN/(dydpT), dN/(dydpT)_err)
        """
        print("processing particle spectra for %s ... " % particle_name)
        eps = 1e-15
        pid = self.pid_lookup[particle_name]
        if rap_type == 'rapidity':
            analyzed_table_name = 'particle_pT_spectra'
        elif rap_type == 'pseudorapidity':
            analyzed_table_name = 'particle_pT_spectra_eta'
        else:
            raise ValueError("unrecognized rap_type: %s" % rap_type)

        npT = len(self.db.executeSQLquery(
            "select pT from %s where hydro_event_id = %d and "
            "urqmd_event_id = %d and pid = %d" 
            % (analyzed_table_name, 1, 1, pid)).fetchall())

        dN_avg = zeros([npT, 3])
        #fetch data
        for ibin in range(1, self.nev_bin):
            print("processing events %d to %d ..." 
                % ((ibin-1)*self.process_nev, ibin*self.process_nev))
            hydro_ev_bound_low = self.event_bound_hydro[ibin-1]
            hydro_ev_bound_high = self.event_bound_hydro[ibin]
            urqmd_ev_bound_low = self.event_bound_urqmd[ibin-1]
            urqmd_ev_bound_high = self.event_bound_urqmd[ibin]
            if hydro_ev_bound_low == hydro_ev_bound_high:
                temp_data = array(self.db.executeSQLquery(
                    "select pT, dN_dydpT from %s where pid = %d and "
                    "hydro_event_id = %d and (%d <= urqmd_event_id "
                    "and urqmd_event_id < %d)"
                    % (analyzed_table_name, pid, hydro_ev_bound_low, 
                       urqmd_ev_bound_low, urqmd_ev_bound_high)).fetchall())
            else:
                temp_data = array(self.db.executeSQLquery(
                    "select pT, dN_dydpT from %s where pid = %d and "
                    "((hydro_event_id = %d and urqmd_event_id >= %d) "
                    " or (%d < hydro_event_id and hydro_event_id < %d) "
                    " or (hydro_event_id = %d and urqmd_event_id < %d))"
                    % (analyzed_table_name, pid, hydro_ev_bound_low, 
                       urqmd_ev_bound_low, hydro_ev_bound_low, 
                       hydro_ev_bound_high, hydro_ev_bound_high, 
                       urqmd_ev_bound_high)).fetchall())
            for ipT in range(npT):
                dN_avg[ipT,0] += (
                    sum(temp_data[ipT::npT,0]*temp_data[ipT::npT,1]))
                dN_avg[ipT,1] += sum(temp_data[ipT::npT,1])
                dN_avg[ipT,2] += sum(temp_data[ipT::npT,1]**2)
        
        # calculate mean pT, <dN/dydpT>, and <dN/dydpT>_err 
        dN_avg[:,0] = dN_avg[:,0]/dN_avg[:,1]
        dN_avg[:,1] = dN_avg[:,1]/self.tot_nev
        dN_avg[:,2] = (sqrt(dN_avg[:,2]/self.tot_nev - dN_avg[:,1]**2)
                       /sqrt(self.tot_nev))
        
        #interpolate results to desired pT range
        dNdyinterp = exp(interp(pT_range, dN_avg[:,0], log(dN_avg[:,1]+eps)))
        dNdyinterp_err = exp(
            interp(pT_range, dN_avg[:,0], log(dN_avg[:,2]+eps)))
        results = array([pT_range, dNdyinterp, dNdyinterp_err])
        return transpose(results)

    def get_particle_yield_vs_rap(self, particle_name, rap_type = 'rapidity', 
                                  rap_range = linspace(-2.0, 2.0, 41)):
        """
            It returns event averaged particle yield vs rapidity
            or pseudorapidity. The range is specified by the user.
            It returns (rap, dN/(drap), dN/(drap)_err)
        """
        pid = self.pid_lookup[particle_name]
        if rap_type == 'rapidity':
            analyzed_table_name = 'particle_yield_vs_rap'
        elif rap_type == 'pseudorapidity':
            analyzed_table_name = 'particle_yield_vs_psedurap'
        else:
            raise ValueError("unrecognized rap_type: %s" % rap_type)
        
        nrap = len(self.db.executeSQLquery(
            "select rap from %s where hydro_event_id = %d and "
            "urqmd_event_id = %d and pid = %d" 
            % (analyzed_table_name, 1, 1, pid)).fetchall())
        
        dN_avg = zeros([nrap, 3])
        #fetch data
        for ibin in range(1, self.nev_bin):
            print("processing events %d to %d ..." 
                % ((ibin-1)*self.process_nev, ibin*self.process_nev))
            hydro_ev_bound_low = self.event_bound_hydro[ibin-1]
            hydro_ev_bound_high = self.event_bound_hydro[ibin]
            urqmd_ev_bound_low = self.event_bound_urqmd[ibin-1]
            urqmd_ev_bound_high = self.event_bound_urqmd[ibin]
            if hydro_ev_bound_low == hydro_ev_bound_high:
                temp_data = array(self.db.executeSQLquery(
                    "select rap, dN_drap from %s where pid = %d and "
                    "hydro_event_id = %d and (%d <= urqmd_event_id "
                    "and urqmd_event_id < %d)"
                    % (analyzed_table_name, pid, hydro_ev_bound_low, 
                       urqmd_ev_bound_low, urqmd_ev_bound_high)).fetchall())
            else:
                temp_data = array(self.db.executeSQLquery(
                    "select rap, dN_drap from %s where pid = %d and "
                    "((hydro_event_id = %d and urqmd_event_id >= %d) "
                    " or (%d < hydro_event_id and hydro_event_id < %d) "
                    " or (hydro_event_id = %d and urqmd_event_id < %d))"
                    % (analyzed_table_name, pid, hydro_ev_bound_low, 
                       urqmd_ev_bound_low, hydro_ev_bound_low, 
                       hydro_ev_bound_high, hydro_ev_bound_high, 
                       urqmd_ev_bound_high)).fetchall())
            for irap in range(nrap):
                dN_avg[irap,0] += sum(
                    temp_data[irap::nrap,0]*temp_data[irap::nrap,1])
                dN_avg[irap,1] += sum(temp_data[irap::nrap,1])
                dN_avg[irap,2] += sum(temp_data[irap::nrap,1]**2)

        # calculate mean rap, <dN/drap>, and <dN/drap>_err 
        dN_avg[:,0] = dN_avg[:,0]/dN_avg[:,1]
        dN_avg[:,1] = dN_avg[:,1]/self.tot_nev
        dN_avg[:,2] = (sqrt(dN_avg[:,2]/self.tot_nev - dN_avg[:,1]**2)
                       /sqrt(self.tot_nev))
        
        #interpolate results to desired rap range
        dNdyinterp = interp(rap_range, dN_avg[:,0], dN_avg[:,1])
        dNdyinterp_err = interp(rap_range, dN_avg[:,0], dN_avg[:,2])
        results = array([rap_range, dNdyinterp, dNdyinterp_err])
        return transpose(results)

    def get_particle_yield(self, particle_name, rap_type = 'rapidity', 
                           rap_range = (-0.5, 0.5)):
        """
            return the integrated particle yield N of particle species 
            "particle_name" within the given rapidity or pseudorapidity 
            range by users
        """
        pid = self.pid_lookup[particle_name]
        if rap_type == 'rapidity':
            analyzed_table_name = 'particle_yield_vs_rap'
        elif rap_type == 'pseudorapidity':
            analyzed_table_name = 'particle_yield_vs_psedurap'
        else:
            raise ValueError("unrecognized rap_type: %s" % rap_type)
        
        nrap = len(self.db.executeSQLquery(
            "select rap from %s where hydro_event_id = %d and "
            "urqmd_event_id = %d and pid = %d" 
            % (analyzed_table_name, 1, 1, pid)).fetchall())
        drap = 0.1
        
        mean_rap = 0.0
        particle_yield = 0.0
        particle_yield_err = 0.0
        
        #fetch data
        for ibin in range(1, self.nev_bin):
            print("processing events %d to %d ..." 
                % ((ibin-1)*self.process_nev, ibin*self.process_nev))
            hydro_ev_bound_low = self.event_bound_hydro[ibin-1]
            hydro_ev_bound_high = self.event_bound_hydro[ibin]
            urqmd_ev_bound_low = self.event_bound_urqmd[ibin-1]
            urqmd_ev_bound_high = self.event_bound_urqmd[ibin]
            if hydro_ev_bound_low == hydro_ev_bound_high:
                temp_data = array(self.db.executeSQLquery(
                    "select rap, dN_drap from %s where pid = %d and "
                    "hydro_event_id = %d and (%d <= urqmd_event_id "
                    "and urqmd_event_id < %d)"
                    % (analyzed_table_name, pid, hydro_ev_bound_low, 
                       urqmd_ev_bound_low, urqmd_ev_bound_high)).fetchall())
            else:
                temp_data = array(self.db.executeSQLquery(
                    "select rap, dN_drap from %s where pid = %d and "
                    "((hydro_event_id = %d and urqmd_event_id >= %d) "
                    " or (%d < hydro_event_id and hydro_event_id < %d) "
                    " or (hydro_event_id = %d and urqmd_event_id < %d))"
                    % (analyzed_table_name, pid, hydro_ev_bound_low, 
                       urqmd_ev_bound_low, hydro_ev_bound_low, 
                       hydro_ev_bound_high, hydro_ev_bound_high, 
                       urqmd_ev_bound_high)).fetchall())
            cut_data = temp_data[temp_data[:,0] > rap_range[0],:]
            cut_data2 = cut_data[cut_data[:,0] < rap_range[1],:]
            mean_rap += sum(cut_data2[:,0]*cut_data2[:,1])
            particle_yield += sum(cut_data2[:,1])
            particle_yield_err += sum(cut_data2[:,1]**2)

        # calculate mean rap, <dN>, and <dN>_err 
        mean_rap = mean_rap/particle_yield
        particle_yield = particle_yield/self.tot_nev*drap
        particle_yield_err = (
            sqrt(particle_yield_err/self.tot_nev - particle_yield**2)
            /sqrt(self.tot_nev))*drap
        return(mean_rap, particle_yield, particle_yield_err)

    ###########################################################################
    # functions to collect particle emission function
    ########################################################################### 
    def get_particle_yield_vs_spatial_variable(
        self, particle_name, sv_type, sv_range, rap_type):
        """
            This function performs event average for particle yield per
            spatial variable from the database and interpolates to the 
            desired sv values specified by the users.
            It returns (sv, dN/(dsv), dN/(dsv)_err)
        """
        eps = 1e-15
        pid = self.pid_lookup[particle_name]
        if rap_type == 'rapidity':
            analyzed_table_name = 'particle_emission_d%s' % sv_type
        elif rap_type == 'pseudorapidity':
            analyzed_table_name = 'particle_emission_d%s_eta' % sv_type
        else:
            raise ValueError("AnalyzedDataReader.get_particle_yield_vs_"
                            "spatial_variable: invalid input rap_type : %s" 
                            % rap_type)

        nsv = len(self.db.executeSQLquery(
            "select %s from %s where hydro_event_id = %d and "
            "urqmd_event_id = %d and pid = %d" 
            % (sv_type, analyzed_table_name, 1, 1, pid)).fetchall())

        dN_avg = zeros([nsv, 3])

        #fetch data
        for ibin in range(1, self.nev_bin):
            print("processing events %d to %d ..." 
                % ((ibin-1)*self.process_nev, ibin*self.process_nev))
            hydro_ev_bound_low = self.event_bound_hydro[ibin-1]
            hydro_ev_bound_high = self.event_bound_hydro[ibin]
            urqmd_ev_bound_low = self.event_bound_urqmd[ibin-1]
            urqmd_ev_bound_high = self.event_bound_urqmd[ibin]
            if hydro_ev_bound_low == hydro_ev_bound_high:
                temp_data = array(self.db.executeSQLquery(
                    "select %s, dN_d%s from %s where pid = %d and "
                    "hydro_event_id = %d and (%d <= urqmd_event_id "
                    "and urqmd_event_id < %d)"
                    % (sv_type, sv_type, analyzed_table_name, pid, 
                       hydro_ev_bound_low, urqmd_ev_bound_low, 
                       urqmd_ev_bound_high)
                ).fetchall())
            else:
                temp_data = array(self.db.executeSQLquery(
                    "select %s, dN_d%s from %s where pid = %d and "
                    "((hydro_event_id = %d and urqmd_event_id >= %d) "
                    " or (%d < hydro_event_id and hydro_event_id < %d) "
                    " or (hydro_event_id = %d and urqmd_event_id < %d))"
                    % (sv_type, sv_type, analyzed_table_name, pid, 
                       hydro_ev_bound_low, urqmd_ev_bound_low, 
                       hydro_ev_bound_low, hydro_ev_bound_high, 
                       hydro_ev_bound_high, urqmd_ev_bound_high)
                ).fetchall())
            for isv in range(nsv):
                dN_avg[isv,0] += sum(
                    temp_data[isv::nsv,0]*temp_data[isv::nsv,1])
                dN_avg[isv,1] += sum(temp_data[isv::nsv,1])
                dN_avg[isv,2] += sum(temp_data[isv::nsv,1]**2)
        
        # calculate <sv>, <dN/dsv>, and <dN/dsv>_err 
        dN_avg[:,0] = dN_avg[:,0]/(dN_avg[:,1] + eps)
        dN_avg[:,1] = dN_avg[:,1]/self.tot_nev
        dN_avg[:,2] = (sqrt(dN_avg[:,2]/self.tot_nev - dN_avg[:,1]**2)
                       /sqrt(self.tot_nev))
        
        #interpolate results to desired sv range
        dNdyinterp = interp(sv_range, dN_avg[:,0], dN_avg[:,1])
        dNdyinterp_err = interp(sv_range, dN_avg[:,0], dN_avg[:,2])
        results = array([sv_range, dNdyinterp, dNdyinterp_err])
        return transpose(results)

    ###########################################################################
    # functions to collect particle anisotropic flows
    ########################################################################### 
    def get_avg_diffvn_flow(
        self, particle_name, order, psi_r=0., pT_range=linspace(0.0, 3.0, 31)):
        """
            compute <cos(n*(phi_i - psiR))>_ev over all events and interpolate
            the results to desired pT points given by the user
        """
        print("collect averaged diff vn flow of %s ..." % particle_name)
        pid = self.pid_lookup[particle_name]
        analyzed_table_name = 'flow_Qn_vectors_pTdiff'
        npT = len(self.db.executeSQLquery(
            "select pT from %s where hydro_event_id = %d and "
            "urqmd_event_id = %d and pid = %d and weight_type = '1' and n = %d" 
            % (analyzed_table_name, 1, 1, pid, order)).fetchall())

        vn_avg = zeros([npT, 3])
        vn_real = zeros(npT)
        vn_imag = zeros(npT)
        vn_real_err = zeros(npT)
        vn_imag_err = zeros(npT)
        totalN = zeros(npT)
        nev_pT = zeros(npT)
        #fetch data
        for hydroId in range(1, self.hydro_nev+1):
            urqmd_nev = self.db.executeSQLquery(
                "select Number_of_UrQMDevents from UrQMD_NevList where "
                "hydroEventId = %d " % hydroId).fetchall()[0][0]
            for urqmdId in range(1, urqmd_nev+1):
                temp_data = array(self.db.executeSQLquery(
                    "select pT, Nparticle, Qn, Qn_psi from %s "
                    "where hydro_event_id = %d and urqmd_event_id = %d and "
                    "pid = %d and weight_type = '1' and n = %d" 
                    % (analyzed_table_name, hydroId, urqmdId, pid, order)
                ).fetchall())
                
                vn_avg[:,0] += temp_data[:,0]*temp_data[:,1] #<pT>
                totalN += temp_data[:,1]
                for ipT in range(npT):
                    if(temp_data[ipT,1] > 0): 
                        nev_pT[ipT] += 1
                        vn_real[ipT] += (temp_data[ipT,2]
                                        *cos(order*(temp_data[ipT,3] - psi_r)))
                        vn_imag[ipT] += (temp_data[ipT,2]
                                        *sin(order*(temp_data[ipT,3] - psi_r)))
                        vn_real_err[ipT] += (
                            temp_data[ipT,2]*cos(
                                order*(temp_data[ipT,3] - psi_r)))**2.
                        vn_imag_err[ipT] += (
                            temp_data[ipT,2]*sin(
                                order*(temp_data[ipT,3] - psi_r)))**2.
        vn_avg[:,0] = vn_avg[:,0]/totalN
        vn_real = vn_real/nev_pT
        vn_imag = vn_imag/nev_pT
        vn_real_err = sqrt(vn_real_err/nev_pT - vn_real**2)/sqrt(nev_pT)
        vn_imag_err = sqrt(vn_imag_err/nev_pT - vn_imag**2)/sqrt(nev_pT)
        vn_avg[:,1] = sqrt(vn_real**2. + vn_imag**2.)
        vn_avg[:,2] = sqrt(vn_real_err**2. + vn_imag_err**2.)/vn_avg[:,1]
        
        #interpolate results to desired pT range
        vn_avg_interp = interp(pT_range, vn_avg[:,0], vn_avg[:,1])
        vn_avg_interp_err = interp(pT_range, vn_avg[:,0], vn_avg[:,2])
        results = array([pT_range, vn_avg_interp, vn_avg_interp_err])
        return transpose(results)
        
    def get_avg_intevn_flow(
        self, particle_name, order, psi_r=0., pT_range=(0.0, 3.0)):
        """
            compute pT integrated <cos(n*(phi_i - psiR))>_ev averaged over 
            all events the pT integrated range is given by the user
        """
        print("collect averaged inte vn flow of %s ..." % particle_name)
        pid = self.pid_lookup[particle_name]
        analyzed_table_name = 'flow_Qn_vectors_pTdiff'

        vn_avg = zeros(3)
        vn_real = 0.0
        vn_imag = 0.0
        vn_real_err = 0.0
        vn_imag_err = 0.0
        totalN = 0
        nev = 0
        #fetch data
        for hydroId in range(1, self.hydro_nev+1):
            urqmd_nev = self.db.executeSQLquery(
                "select Number_of_UrQMDevents from UrQMD_NevList where "
                "hydroEventId = %d " % hydroId).fetchall()[0][0]
            for urqmdId in range(1, urqmd_nev+1):
                temp_data = array(self.db.executeSQLquery(
                    "select pT, Nparticle, Qn, Qn_psi from %s "
                    "where hydro_event_id = %d and urqmd_event_id = %d and "
                    "pid = %d and weight_type = '1' and n = %d and "
                    "(%g <= pT and pT <= %g)"
                    % (analyzed_table_name, hydroId, urqmdId, pid, order,
                       pT_range[0], pT_range[1])
                ).fetchall())
               
                vn_avg[0] += sum(temp_data[:,0]*temp_data[:,1]) #<pT>
                nparticle = sum(temp_data[:,1])
                totalN += nparticle
                if nparticle > 0:
                    nev += 1
                    vn_real += (sum(temp_data[:,1]*temp_data[:,2]
                                    *cos(order*(temp_data[:,3] - psi_r)))
                                    /nparticle)
                    vn_imag += (sum(temp_data[:,1]*temp_data[:,2]
                                    *sin(order*(temp_data[:,3] - psi_r)))
                                    /nparticle)
                    vn_real_err += ((sum(temp_data[:,1]*temp_data[:,2]
                                     *cos(order*(temp_data[:,3] - psi_r))))**2.
                                     /nparticle)
                    vn_imag_err += ((sum(temp_data[:,1]*temp_data[:,2]
                                     *sin(order*(temp_data[:,3] - psi_r))))**2.
                                     /nparticle)
        vn_avg[0] = vn_avg[0]/totalN
        vn_real = vn_real/nev
        vn_imag = vn_imag/nev
        vn_real_err = sqrt(vn_real_err/totalN - vn_real**2)/sqrt(self.tot_nev)
        vn_imag_err = sqrt(vn_imag_err/totalN - vn_imag**2)/sqrt(self.tot_nev)
        vn_avg[1] = sqrt(vn_real**2. + vn_imag**2.)
        vn_avg[2] = sqrt(vn_real_err**2. + vn_imag_err**2.)/vn_avg[1]
        
        return vn_avg

#    def getFullplaneResolutionFactor(self, resolutionFactor_sub, Nfactor):
#        """
#            use binary search for numerical solution for R(chi_s) = resolutionFactor_sub
#            where R(chi) is defined in Eq. (7) in arXiv:0904.2315v3 and calculate 
#            resolutionFactor for the full event
#        """
#        # check
#        if(resolutionFactor_sub > 1.0):
#            print("error: resolutionFactor_sub = % g,  is larger than 1!" % resolutionFactor_sub)
#            exit(-1)
#        
#        tol = 1e-8 # accuracy
#        
#        #search boundary
#        left = 0.0; right = 2.0 # R(2.0) > 1.0
#        mid = (right + left)*0.5
#        dis = right - left
#        while dis > tol:
#            midval = self.resolution_Rfunction(mid)
#            diff = resolutionFactor_sub - midval
#            if abs(diff) < tol:
#                chi_s = mid
#                break
#            elif diff > tol :
#                left = mid
#            else:
#                right = mid
#            dis = right - left
#            mid = (right + left)*0.5
#        chi_s = mid
#        return(self.resolution_Rfunction(chi_s*sqrt(Nfactor)))
#        
#    def resolution_Rfunction(self, chi):
#        """
#            R(chi) for calculating resolution factor for the full event
#        """
#        chisq = chi*chi
#        result = sqrt(pi)/2*exp(-chisq/2)*chi*(scipy.special.i0(chisq/2) + scipy.special.i1(chisq/2))
#        return result
#
#    def calculcateResolutionfactor(self, particleName, order_range = [1, 6], oversampling = 1, weightType = 'pT', pT_range = [0.0, 3.0], rapType = "rapidity", rap_range = [-4.0, 4.0]):
#        """
#            calculate nth order full resolution factor using given species of
#            particles within pT and rapidity cuts from all events with option 
#            of oversampling.
#        """
#        Nev = self.totNev
#        pidString = self.getPidString(particleName)
#        num_of_order  = order_range[1] - order_range[0] + 1
#        resolutionFactor_sub = zeros(num_of_order)
#        for hydroId in range(1, self.hydroNev+1):
#            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
#            for UrQMDId in range(1, UrQMDNev+1):
#                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
#                particleList = array(self.db.executeSQLquery("select pT, phi_p from particle_list where hydroEvent_id = %d and UrQMDEvent_id = %d and (%s) and (%g <= pT and pT < %g) and (%g <= %s and %s <= %g)" % (hydroId, UrQMDId, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
#
#                #calculate resolution factor
#                pT = particleList[:,0]
#                phi = particleList[:,1]
#                if weightType == '1':
#                    weight = ones(len(pT))
#                elif weightType == 'pT':
#                    weight = pT
#                Nparticle = len(pT)
#                for iorder in range(num_of_order):
#                    order = order_range[0] + iorder
#                    for isample in range(oversampling):
#                        idx = range(Nparticle); shuffle(idx)
#                        idx_A = idx[0:int(Nparticle/2)]
#                        shuffle(idx)
#                        idx_B = idx[0:int(Nparticle/2)]
#
#                        subQnA_X = 0.0; subQnA_Y = 0.0
#                        for i in idx_A:
#                            subQnA_X += weight[i]*cos(order*phi[i])
#                            subQnA_Y += weight[i]*sin(order*phi[i])
#                        subPsi_nA = arctan2(subQnA_Y, subQnA_X)/order
#                        subQnB_X = 0.0; subQnB_Y = 0.0
#                        for i in idx_B:
#                            subQnB_X += weight[i]*cos(order*phi[i])
#                            subQnB_Y += weight[i]*sin(order*phi[i])
#                        subPsi_nB = arctan2(subQnB_Y, subQnB_X)/order
#                        resolutionFactor_sub[iorder] += cos(order*(subPsi_nA - subPsi_nB))
#        resolutionFactor_sub = sqrt(resolutionFactor_sub/Nev/oversampling)
#        resolutionFactor_full = []
#        for iorder in range(num_of_order):
#            resolutionFactor_full.append(self.getFullplaneResolutionFactor(resolutionFactor_sub[iorder], 2.0))
#        return(resolutionFactor_full)
#                    
#    def collectGlobalResolutionFactor(self):
#        """
#            collect and store the full resolution factors calculated using all charged particles
#            from all the events for order n = 1-6
#            R_sub = sqrt(<QnA*conj(QnB)/(|QnA||QnB|)>_ev), R_SP_sub = sqrt(<QnA*conj(QnB)>_ev)
#        """
#        weightTypes = ['1', 'pT']
#        norder = 6
#        if self.db.createTableIfNotExists("resolutionFactorR", (("weightType", "text"), ("n", "integer"), ("R","real"), ("R_sub", "real"), ("R_SP_sub", "real"))):
#            for iorder in range(1, norder+1):
#                for weightType in weightTypes:
#                    subQn_data = array(self.db.executeSQLquery("select subQnA_psi, subQnB_psi, subQnA, subQnB from globalQnvector where weightType = '%s' and n = %d" % (weightType, iorder)).fetchall())
#                    resolutionFactor_sub = sqrt(mean(cos(iorder*(subQn_data[:,0] - subQn_data[:,1]))))
#                    resolutionFactor_SP_sub = sqrt(mean(subQn_data[:,2]*subQn_data[:,3]*cos(iorder*(subQn_data[:,0] - subQn_data[:,1]))))
#                    resolutionFactor_full = self.getFullplaneResolutionFactor(resolutionFactor_sub, 2.0)
#                    self.db.insertIntoTable("resolutionFactorR", (weightType, iorder, resolutionFactor_full, resolutionFactor_sub, resolutionFactor_SP_sub))
#            self.db._dbCon.commit()  # commit changes
#        else:
#            print("Resolution factors from all charged particles are already collected!")
#            inputval = raw_input("Do you want to delete the existing one and collect again?")
#            if inputval.lower() == 'y' or inputval.lower() == 'yes':
#                self.db.executeSQLquery("drop table resolutionFactorR")
#                self.collectGlobalResolutionFactor()
#
#    def collectinteEventplaneflow(self, particleName = 'pion_p', pT_range = [0.0, 4.0], rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            collect nth order event plane harmonic flow, vn{EP}
#            event plane angle is determined by all charged particles in the whole event
#        """
#        Nev = self.totNev
#        weightTypes = ['1', 'pT']
#        norder = 6
#        pidString = self.getPidString(particleName)
#        vn_obs = zeros(norder); vn_obs_pTweight = zeros(norder)
#        vn_obs_sq = zeros(norder); vn_obs_pTweight_sq = zeros(norder)
#        for hydroId in range(1, self.hydroNev+1):
#            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
#            for UrQMDId in range(1, UrQMDNev+1):
#                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
#                particleList = array(self.db.executeSQLquery("select pT, phi_p from particle_list where hydroEvent_id = %d and UrQMDEvent_id = %d and (%s) and (%g <= pT and pT < %g) and (%g <= %s and %s <= %g)" % (hydroId, UrQMDId, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
#
#                pT = particleList[:,0]
#                phi = particleList[:,1]
#                Nparticle = len(pT)
#                for weightType in weightTypes:
#                    if weightType == '1':
#                        weight = ones(len(pT))
#                    elif weightType == 'pT':
#                        weight = pT
#                    for iorder in range(1, norder + 1):
#                        QnVector = array(self.db.executeSQLquery("select Nparticle, Qn, Qn_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = '%s' and n = %d" % (hydroId, UrQMDId, weightType, iorder)).fetchall())
#                        Qn_X = QnVector[0,0]*QnVector[0,1]*cos(iorder*QnVector[0,2])
#                        Qn_Y = QnVector[0,0]*QnVector[0,1]*sin(iorder*QnVector[0,2])
#
#                        #calculate event-plane vn
#                        vntemp = 0.0
#                        for ipart in range(Nparticle):
#                            # subtract self correlation
#                            Xn = Qn_X - weight[ipart]*cos(iorder*phi[ipart])
#                            Yn = Qn_Y - weight[ipart]*sin(iorder*phi[ipart])
#                            psi_n = arctan2(Yn, Xn)/iorder
#                            vntemp += cos(iorder*(phi[ipart] - psi_n))
#                        if weightType == '1':
#                            vn_obs[iorder-1] += vntemp/Nparticle
#                            vn_obs_sq[iorder-1] += (vntemp/Nparticle)**2
#                        elif weightType == 'pT':
#                            vn_obs_pTweight[iorder-1] += vntemp/Nparticle
#                            vn_obs_pTweight_sq[iorder-1] += (vntemp/Nparticle)**2
#        vn_obs /= Nev
#        vn_obs_err = sqrt(vn_obs_sq/Nev - (vn_obs)**2.)/sqrt(Nev)
#        vn_obs_pTweight /= Nev
#        vn_obs_pTweight_err = sqrt(vn_obs_pTweight_sq/Nev - (vn_obs_pTweight)**2.)/sqrt(Nev)
#        vnEP = zeros(norder); vnEP_pTweight = zeros(norder)
#        vnEP_err = zeros(norder); vnEP_pTweight_err = zeros(norder)
#        for iorder in range(1, norder + 1):
#            resolutionFactor = self.db.executeSQLquery("select R from resolutionFactorR where weightType = '1' and n = %d" % (iorder)).fetchall()[0][0]
#            resolutionFactor_pTweight = self.db.executeSQLquery("select R from resolutionFactorR where weightType = 'pT' and n = %d" % (iorder)).fetchall()[0][0]
#            vnEP[iorder-1] = vn_obs[iorder-1]/resolutionFactor
#            vnEP_err[iorder-1] = vn_obs_err[iorder-1]/resolutionFactor
#            vnEP_pTweight[iorder-1] = vn_obs_pTweight[iorder-1]/resolutionFactor_pTweight
#            vnEP_pTweight_err[iorder-1] = vn_obs_pTweight_err[iorder-1]/resolutionFactor_pTweight
#        
#        return(vnEP, vnEP_err, vnEP_pTweight, vnEP_pTweight_err)
#
#    def collectEventplaneflow(self, particleName = 'pion_p', pT_range = [0.0, 3.0], rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            collect nth order pT differential event plane and subevent plane harmonic 
#            flow, vn{EP}(pT) and vn{subEP}(pT). 
#            event plane angle is determined by all charged particles in the whole event
#            subevent plane angle is determined by all charged particles in one sub event
#            with eta gap
#        """
#        print("Collecting differential event plane flow vn{EP} and vn{subEP} for %s ..." % particleName)
#        Nev = self.totNev
#        norder = 6; npT = 30
#        weightTypes = ['1', 'pT']
#        pidString = self.getPidString(particleName)
#
#        Nev_vninte = 0.0
#        vn_inte_obs = zeros([norder]); vn_inte_obs_pTweight = zeros([norder])
#        vn_inte_obs_sq = zeros([norder]); vn_inte_obs_pTweight_sq = zeros([norder])
#        vn_inte_obs_sub = zeros([norder]); vn_inte_obs_pTweight_sub = zeros([norder])
#        vn_inte_obs_sq_sub = zeros([norder]); vn_inte_obs_pTweight_sq_sub = zeros([norder])
#        vn_inte_obs_err = zeros([norder]); vn_inte_obs_pTweight_err = zeros([norder])
#        vn_inte_obs_sub_err = zeros([norder]); vn_inte_obs_pTweight_sub_err = zeros([norder])
#
#        pTBoundary = linspace(pT_range[0], pT_range[1], npT+1)
#        pTmean = zeros([npT]); NevpT = zeros([npT])
#        vn_obs = zeros([norder, npT]); vn_obs_pTweight = zeros([norder, npT])
#        vn_obs_sq = zeros([norder, npT]); vn_obs_pTweight_sq = zeros([norder, npT])
#        vn_obs_sub = zeros([norder, npT]); vn_obs_pTweight_sub = zeros([norder, npT])
#        vn_obs_sq_sub = zeros([norder, npT]); vn_obs_pTweight_sq_sub = zeros([norder, npT])
#        vn_obs_err = zeros([norder, npT]); vn_obs_pTweight_err = zeros([norder, npT])
#        vn_obs_sub_err = zeros([norder, npT]); vn_obs_pTweight_sub_err = zeros([norder, npT])
#        
#        for hydroId in range(1, self.hydroNev+1):
#            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
#            for UrQMDId in range(1, UrQMDNev+1):
#                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
#                particleList = array(self.db.executeSQLquery("select pT, phi_p from particle_list where hydroEvent_id = %d and UrQMDEvent_id = %d and (%s) and (%g <= pT and pT < %g) and (%g <= %s and %s <= %g)" % (hydroId, UrQMDId, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
#
#                Qn_X = zeros([2, norder]); Qn_Y = zeros([2, norder])
#                subQnA_psi = zeros([2, norder])
#                for iorder in range(1, norder + 1):
#                    QnVector = array(self.db.executeSQLquery("select Nparticle, Qn, Qn_psi, subQnA_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = '1' and n = %d" % (hydroId, UrQMDId, iorder)).fetchall())
#                    Qn_X[0, iorder-1] = QnVector[0,0]*QnVector[0,1]*cos(iorder*QnVector[0,2])
#                    Qn_Y[0, iorder-1] = QnVector[0,0]*QnVector[0,1]*sin(iorder*QnVector[0,2])
#                    subQnA_psi[0, iorder-1] = QnVector[0,3]
#                    QnVector = array(self.db.executeSQLquery("select Nparticle, Qn, Qn_psi, subQnA_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = 'pT' and n = %d" % (hydroId, UrQMDId, iorder)).fetchall())
#                    Qn_X[1, iorder-1] = QnVector[0,0]*QnVector[0,1]*cos(iorder*QnVector[0,2])
#                    Qn_Y[1, iorder-1] = QnVector[0,0]*QnVector[0,1]*sin(iorder*QnVector[0,2])
#                    subQnA_psi[1, iorder-1] = QnVector[0,3]
#                
#                # calculate pT integrated vn{EP} and vn{subEP}
#                if(particleList.size != 0):
#                    Nev_vninte += 1
#                    pT = particleList[:,0]
#                    phi = particleList[:,1]
#                    Nparticle = len(pT)
#                    for weightType in weightTypes:
#                        if weightType == '1':
#                            weight = ones(Nparticle); iweight = 0
#                        elif weightType == 'pT':
#                            weight = pT; iweight = 1
#                        for iorder in range(1, norder + 1):
#                            Qn_X_local = Qn_X[iweight, iorder-1]
#                            Qn_Y_local = Qn_Y[iweight, iorder-1]
#                            subPsi_n = subQnA_psi[iweight, iorder-1]
#
#                            #calculate event-plane vn
#                            vntemp = 0.0
#                            for ipart in range(Nparticle):
#                                # subtract self correlation
#                                Xn = Qn_X_local - weight[ipart]*cos(iorder*phi[ipart])
#                                Yn = Qn_Y_local - weight[ipart]*sin(iorder*phi[ipart])
#                                psi_n = arctan2(Yn, Xn)/iorder
#                                vntemp += weight[ipart]*cos(iorder*(phi[ipart] - psi_n))
#                            vntemp = vntemp/Nparticle
#                            subvntemp = sum(weight*cos(iorder*(phi - subPsi_n)))/Nparticle
#                            if weightType == '1':
#                                vn_inte_obs[iorder-1] += vntemp
#                                vn_inte_obs_sq[iorder-1] += vntemp**2
#                                vn_inte_obs_sub[iorder-1] += subvntemp
#                                vn_inte_obs_sq_sub[iorder-1] += subvntemp**2
#                            elif weightType == 'pT':
#                                vn_inte_obs_pTweight[iorder-1] += vntemp
#                                vn_inte_obs_pTweight_sq[iorder-1] += vntemp**2
#                                vn_inte_obs_pTweight_sub[iorder-1] += subvntemp
#                                vn_inte_obs_pTweight_sq_sub[iorder-1] += subvntemp**2
#
#                #calculate pT differential vn{EP} and vn{subEP}
#                for ipT in range(npT):
#                    pTlow = pTBoundary[ipT]; pThigh = pTBoundary[ipT+1]
#                    tempList = particleList[particleList[:,0]<pThigh,:]
#                    particleList_pTcut = tempList[tempList[:,0]>=pTlow,:]
#                    
#                    if(particleList_pTcut.size == 0):  # no particle in the event
#                        pTmean[ipT] = (pTlow + pThigh)/2.
#                    else:
#                        NevpT[ipT] += 1
#                        pT = particleList_pTcut[:,0]
#                        phi = particleList_pTcut[:,1]
#                        pTmean[ipT] = mean(pT)
#                        Nparticle = len(pT)
#                        for weightType in weightTypes:
#                            if weightType == '1':
#                                weight = ones(len(pT))
#                                iweight = 0
#                            elif weightType == 'pT':
#                                weight = pT
#                                iweight = 1
#                            for iorder in range(1, norder + 1):
#                                Qn_X_local = Qn_X[iweight, iorder-1]
#                                Qn_Y_local = Qn_Y[iweight, iorder-1]
#                                subPsi_n = subQnA_psi[iweight, iorder-1]
#
#                                #calculate event-plane vn
#                                vntemp = 0.0
#                                for ipart in range(Nparticle):
#                                    # subtract self correlation
#                                    Xn = Qn_X_local - weight[ipart]*cos(iorder*phi[ipart])
#                                    Yn = Qn_Y_local - weight[ipart]*sin(iorder*phi[ipart])
#                                    psi_n = arctan2(Yn, Xn)/iorder
#                                    vntemp += weight[ipart]*cos(iorder*(phi[ipart] - psi_n))
#                                vntemp = vntemp/Nparticle
#                                subvntemp = sum(weight*cos(iorder*(phi - subPsi_n)))/Nparticle
#                                if weightType == '1':
#                                    vn_obs[iorder-1][ipT] += vntemp
#                                    vn_obs_sq[iorder-1][ipT] += vntemp**2
#                                    vn_obs_sub[iorder-1][ipT] += subvntemp
#                                    vn_obs_sq_sub[iorder-1][ipT] += subvntemp**2
#                                elif weightType == 'pT':
#                                    vn_obs_pTweight[iorder-1][ipT] += vntemp
#                                    vn_obs_pTweight_sq[iorder-1][ipT] += vntemp**2
#                                    vn_obs_pTweight_sub[iorder-1][ipT] += subvntemp
#                                    vn_obs_pTweight_sq_sub[iorder-1][ipT] += subvntemp**2
#
#        eps = 1e-15
#        vn_inte_obs = vn_inte_obs/(Nev_vninte + eps)
#        vn_inte_obs_err = sqrt(vn_inte_obs_sq/(Nev_vninte + eps) - (vn_inte_obs)**2.)/sqrt(Nev_vninte + eps)
#        vn_inte_obs_pTweight = vn_inte_obs_pTweight/(Nev_vninte + eps)
#        vn_inte_obs_pTweight_err = sqrt(vn_inte_obs_pTweight_sq/(Nev_vninte + eps) - (vn_inte_obs_pTweight)**2.)/sqrt(Nev_vninte + eps)
#        vn_inte_obs_sub = vn_inte_obs_sub/(Nev_vninte + eps)
#        vn_inte_obs_sub_err = sqrt(vn_inte_obs_sq_sub/(Nev_vninte + eps) - (vn_inte_obs_sub)**2.)/sqrt(Nev_vninte + eps)
#        vn_inte_obs_pTweight_sub = vn_inte_obs_pTweight_sub/(Nev_vninte + eps)
#        vn_inte_obs_pTweight_sub_err = sqrt(vn_inte_obs_pTweight_sq_sub/(Nev_vninte + eps) - (vn_inte_obs_pTweight_sub)**2.)/sqrt(Nev_vninte + eps)
#
#        vn_obs = vn_obs/(NevpT + eps)
#        vn_obs_err = sqrt(vn_obs_sq/(NevpT + eps) - (vn_obs)**2.)/sqrt(NevpT + eps)
#        vn_obs_pTweight = vn_obs_pTweight/(NevpT + eps)
#        vn_obs_pTweight_err = sqrt(vn_obs_pTweight_sq/(NevpT + eps) - (vn_obs_pTweight)**2.)/sqrt(NevpT + eps)
#        vn_obs_sub = vn_obs_sub/(NevpT + eps)
#        vn_obs_sub_err = sqrt(vn_obs_sq_sub/(NevpT + eps) - (vn_obs_sub)**2.)/sqrt(NevpT + eps)
#        vn_obs_pTweight_sub = vn_obs_pTweight_sub/(NevpT + eps)
#        vn_obs_pTweight_sub_err = sqrt(vn_obs_pTweight_sq_sub/(NevpT + eps) - (vn_obs_pTweight_sub)**2.)/sqrt(NevpT + eps)
#
#        vninteEP = zeros([norder]); vninteEP_pTweight = zeros([norder])
#        vninteEP_err = zeros([norder]); vninteEP_pTweight_err = zeros([norder])
#        vnintesubEP = zeros([norder]); vnintesubEP_pTweight = zeros([norder])
#        vnintesubEP_err = zeros([norder]); vnintesubEP_pTweight_err = zeros([norder])
#
#        vnEP = zeros([norder, npT]); vnEP_pTweight = zeros([norder, npT])
#        vnEP_err = zeros([norder, npT]); vnEP_pTweight_err = zeros([norder, npT])
#        vnsubEP = zeros([norder, npT]); vnsubEP_pTweight = zeros([norder, npT])
#        vnsubEP_err = zeros([norder, npT]); vnsubEP_pTweight_err = zeros([norder, npT])
#
#        for iorder in range(1, norder + 1):
#            resolutionFactor = self.db.executeSQLquery("select R from resolutionFactorR where weightType = '1' and n = %d" % (iorder)).fetchall()[0][0]
#            resolutionFactor_pTweight = self.db.executeSQLquery("select R from resolutionFactorR where weightType = 'pT' and n = %d" % (iorder)).fetchall()[0][0]
#            resolutionFactor_sub = self.db.executeSQLquery("select R_sub from resolutionFactorR where weightType = '1' and n = %d" % (iorder)).fetchall()[0][0]
#            resolutionFactor_pTweight_sub = self.db.executeSQLquery("select R_sub from resolutionFactorR where weightType = 'pT' and n = %d" % (iorder)).fetchall()[0][0]
#
#            vninteEP[iorder-1] = vn_inte_obs[iorder-1]/resolutionFactor
#            vninteEP_err[iorder-1] = vn_inte_obs_err[iorder-1]/resolutionFactor
#            vninteEP_pTweight[iorder-1] = vn_inte_obs_pTweight[iorder-1]/resolutionFactor_pTweight
#            vninteEP_pTweight_err[iorder-1] = vn_inte_obs_pTweight_err[iorder-1]/resolutionFactor_pTweight
#            vnintesubEP[iorder-1] = vn_inte_obs_sub[iorder-1]/resolutionFactor_sub
#            vnintesubEP_err[iorder-1] = vn_inte_obs_sub_err[iorder-1]/resolutionFactor_sub
#            vnintesubEP_pTweight[iorder-1] = vn_inte_obs_pTweight_sub[iorder-1]/resolutionFactor_pTweight_sub
#            vnintesubEP_pTweight_err[iorder-1] = vn_inte_obs_pTweight_sub_err[iorder-1]/resolutionFactor_pTweight_sub
#
#            for ipT in range(npT):
#                vnEP[iorder-1, ipT] = vn_obs[iorder-1, ipT]/resolutionFactor
#                vnEP_err[iorder-1, ipT] = vn_obs_err[iorder-1, ipT]/resolutionFactor
#                vnEP_pTweight[iorder-1, ipT] = vn_obs_pTweight[iorder-1, ipT]/resolutionFactor_pTweight
#                vnEP_pTweight_err[iorder-1, ipT] = vn_obs_pTweight_err[iorder-1, ipT]/resolutionFactor_pTweight
#                vnsubEP[iorder-1, ipT] = vn_obs_sub[iorder-1, ipT]/resolutionFactor_sub
#                vnsubEP_err[iorder-1, ipT] = vn_obs_sub_err[iorder-1, ipT]/resolutionFactor_sub
#                vnsubEP_pTweight[iorder-1, ipT] = vn_obs_pTweight_sub[iorder-1, ipT]/resolutionFactor_pTweight_sub
#                vnsubEP_pTweight_err[iorder-1, ipT] = vn_obs_pTweight_sub_err[iorder-1, ipT]/resolutionFactor_pTweight_sub
#        
#        return(vninteEP, vninteEP_err, vninteEP_pTweight, vninteEP_pTweight_err, vnintesubEP, vnintesubEP_err, vnintesubEP_pTweight, vnintesubEP_pTweight_err, pTmean, vnEP, vnEP_err, vnEP_pTweight, vnEP_pTweight_err, vnsubEP, vnsubEP_err, vnsubEP_pTweight, vnsubEP_pTweight_err)
#
#    def getEventplaneflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            retrieve nth order event plane and subevent plane flow data from database 
#            for the given species of particles with given pT_range.
#            if no results is found, it will collect vn{EP}(pT) and vn{subEP}(pT) and 
#            store them into the database
#        """
#        pid = self.pid_lookup[particleName]
#        collectFlag = False
#        if rapType == 'rapidity':
#            tableName = "diffvnEP"; tableName_vninte = "intevnEP"
#        elif rapType == 'pseudorapidity':
#            tableName = "diffvnEPeta"; tableName_vninte = "intevnEPeta"
#        if self.db.createTableIfNotExists(tableName, (("pid", "integer"), ("weightType", "text"), ("n", "integer"), ("pT", "real"), ("vnEP", "real"), ("vnEP_err", "real"), ("vnsubEP", "real"), ("vnsubEP_err", "real") )):
#            collectFlag = True
#        else:
#            vnEPdata = array(self.db.executeSQLquery("select pT, vnEP, vnEP_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName, pid, weightType, order)).fetchall())
#            vnsubEPdata = array(self.db.executeSQLquery("select pT, vnsubEP, vnsubEP_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName, pid, weightType, order)).fetchall())
#            if vnEPdata.size == 0: collectFlag = True
#            if vnsubEPdata.size == 0: collectFlag = True
#        if self.db.createTableIfNotExists(tableName_vninte, (("pid", "integer"), ("weightType", "text"), ("n", "integer"), ("vnEP", "real"), ("vnEP_err", "real"), ("vnsubEP", "real"), ("vnsubEP_err", "real") )):
#            collectFlag = True
#        else:
#            vninteEPdata = array(self.db.executeSQLquery("select vnEP, vnEP_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName_vninte, pid, weightType, order)).fetchall())
#            vnintesubEPdata = array(self.db.executeSQLquery("select vnsubEP, vnsubEP_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName_vninte, pid, weightType, order)).fetchall())
#            if vninteEPdata.size == 0: collectFlag = True
#            if vnintesubEPdata.size == 0: collectFlag = True
#        if collectFlag:
#            vninteEP, vninteEP_err, vninteEP_pTweight, vninteEP_pTweight_err, vnintesubEP, vnintesubEP_err, vnintesubEP_pTweight, vnintesubEP_pTweight_err, pT, vnEP, vnEP_err, vnEP_pTweight, vnEP_pTweight_err, vnsubEP, vnsubEP_err, vnsubEP_pTweight, vnsubEP_pTweight_err = self.collectEventplaneflow(particleName = particleName, rapType = rapType) 
#            for iorder in range(len(vnEP[:,0])):
#                self.db.insertIntoTable(tableName_vninte, (pid, '1', iorder+1, vninteEP[iorder], vninteEP_err[iorder], vnintesubEP[iorder], vnintesubEP_err[iorder]))
#                self.db.insertIntoTable(tableName_vninte, (pid, 'pT', iorder+1, vninteEP_pTweight[iorder], vninteEP_pTweight_err[iorder], vnintesubEP_pTweight[iorder], vnintesubEP_pTweight_err[iorder]) )
#                for ipT in range(len(pT)):
#                    self.db.insertIntoTable(tableName, (pid, '1', iorder+1, pT[ipT], vnEP[iorder, ipT], vnEP_err[iorder, ipT], vnsubEP[iorder, ipT], vnsubEP_err[iorder, ipT]))
#                    self.db.insertIntoTable(tableName, (pid, 'pT', iorder+1, pT[ipT], vnEP_pTweight[iorder, ipT], vnEP_pTweight_err[iorder, ipT], vnsubEP_pTweight[iorder, ipT], vnsubEP_pTweight_err[iorder, ipT]) )
#            self.db._dbCon.commit()
#            if weightType == '1':
#                vninteEPdata = array([vninteEP[order-1], vninteEP_err[order-1]]).transpose()
#                vnintesubEPdata = array([vnsubEP[order-1], vnsubEP_err[order-1]]).transpose()
#                vnEPdata = array([pT, vnEP[order-1,:], vnEP_err[order-1,:]]).transpose()
#                vnsubEPdata = array([pT, vnsubEP[order-1,:], vnsubEP_err[order-1,:]]).transpose()
#            elif weightType == 'pT':
#                vninteEPdata = array([vninteEP_pTweight[order-1], vninteEP_pTweight_err[order-1]]).transpose()
#                vnintesubEPdata = array([vnsubEP_pTweight[order-1], vnsubEP_pTweight_err[order-1]]).transpose()
#                vnEPdata = array([pT, vnEP_pTweight[order-1], vnEP_pTweight_err[order-1,:]]).transpose()
#                vnsubEPdata = array([pT, vnsubEP_pTweight[order-1], vnsubEP_pTweight_err[order-1,:]]).transpose()
#            if vnEPdata.size == 0 or vnsubEPdata.size == 0:
#                print("There is no record for different event plane flow vn for %s" % particleName)
#                return None
#      
#        #interpolate results to desired pT range
#        vnEPpTinterp = interp(pT_range, vnEPdata[:,0], vnEPdata[:,1])
#        vnEPpTinterp_err = interp(pT_range, vnEPdata[:,0], vnEPdata[:,2])
#        vnsubEPpTinterp = interp(pT_range, vnsubEPdata[:,0], vnsubEPdata[:,1])
#        vnsubEPpTinterp_err = interp(pT_range, vnsubEPdata[:,0], vnsubEPdata[:,2])
#        results = array([pT_range, vnEPpTinterp, vnEPpTinterp_err, vnsubEPpTinterp, vnsubEPpTinterp_err])
#        return(array([vninteEPdata, vnintesubEPdata]).transpose(), transpose(results))
#                
#    def getdiffEventplaneflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            return vn{EP} and vn{subEP}
#        """
#        vninteEP, vndiffEP = self.getEventplaneflow(particleName = particleName, order = order, weightType = weightType, pT_range = pT_range, rap_range = rap_range, rapType = rapType)
#        return(vndiffEP)
#    
#    def getinteEventplaneflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            return vn{EP} and vn{subEP}
#        """
#        vninteEP, vndiffEP = self.getEventplaneflow(particleName = particleName, order = order, weightType = weightType, pT_range = pT_range, rap_range = rap_range, rapType = rapType)
#        return(vninteEP)
#
#    def collectScalarProductflow(self, particleName = 'pion_p', pT_range = [0.0, 3.0], rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            collect nth order pT differential scalar product harmonic flow, vn{SP}(pT)
#            event plane angle is determined by all charged particles in the whole event
#        """
#        print("Collecting differential scalar product flow vn{SP} for %s ..." % particleName)
#        Nev = self.totNev
#        norder = 6; npT = 30
#        weightTypes = ['1', 'pT']
#        pidString = self.getPidString(particleName)
#
#        Nev_vninte = 0.0
#        vn_inte_obs = zeros([norder]); vn_inte_obs_pTweight = zeros([norder])
#        vn_inte_obs_sq = zeros([norder]); vn_inte_obs_pTweight_sq = zeros([norder])
#
#        pTBoundary = linspace(pT_range[0], pT_range[1], npT+1)
#        pTmean = zeros(npT); NevpT = zeros(npT)
#        vn_obs = zeros([norder, npT]); vn_obs_pTweight = zeros([norder, npT])
#        vn_obs_sq = zeros([norder, npT]); vn_obs_pTweight_sq = zeros([norder, npT])
#
#        for hydroId in range(1, self.hydroNev+1):
#            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
#            for UrQMDId in range(1, UrQMDNev+1):
#                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
#                particleList = array(self.db.executeSQLquery("select pT, phi_p from particle_list where hydroEvent_id = %d and UrQMDEvent_id = %d and (%s) and (%g <= pT and pT < %g) and (%g <= %s and %s <= %g)" % (hydroId, UrQMDId, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
#
#                QnA_X = zeros([2,norder]); QnA_Y = zeros([2,norder])
#                for iorder in range(1, norder + 1):
#                    subQnVectors = array(self.db.executeSQLquery("select subQnA, subQnA_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = '1' and n = %d" % (hydroId, UrQMDId, iorder)).fetchall())
#                    QnA_X[0,iorder-1] = subQnVectors[0,0]*cos(iorder*subQnVectors[0,1])
#                    QnA_Y[0,iorder-1] = subQnVectors[0,0]*sin(iorder*subQnVectors[0,1])
#                    subQnVectors = array(self.db.executeSQLquery("select subQnA, subQnA_psi from globalQnvector where hydroEvent_id = %d and UrQMDEvent_id = %d and weightType = 'pT' and n = %d" % (hydroId, UrQMDId, iorder)).fetchall())
#                    QnA_X[1,iorder-1] = subQnVectors[0,0]*cos(iorder*subQnVectors[0,1])
#                    QnA_Y[1,iorder-1] = subQnVectors[0,0]*sin(iorder*subQnVectors[0,1])
#              
#                # pT integrate vn{SP}
#                if(particleList.size != 0):
#                    Nev_vninte += 1
#                    pT = particleList[:,0]
#                    phi = particleList[:,1]
#                    Nparticle = len(pT)
#                    for weightType in weightTypes:
#                        if weightType == '1':
#                            weight = ones(Nparticle); iweight = 0
#                        elif weightType == 'pT':
#                            weight = pT; iweight = 1
#                        for iorder in range(1, norder + 1):
#                            QnA_X_local = QnA_X[iweight, iorder-1]
#                            QnA_Y_local = QnA_Y[iweight, iorder-1]
#
#                            # Qn vector for the interested particles
#                            Qn_X = sum(weight*cos(iorder*phi))/Nparticle
#                            Qn_Y = sum(weight*sin(iorder*phi))/Nparticle
#
#                            temp = Qn_X*QnA_X_local + Qn_Y*QnA_Y_local
#                            if weightType == '1':
#                                vn_inte_obs[iorder-1] += temp
#                                vn_inte_obs_sq[iorder-1] += temp*temp
#                            elif weightType == 'pT':
#                                vn_inte_obs_pTweight[iorder-1] += temp
#                                vn_inte_obs_pTweight_sq[iorder-1] += temp*temp
#
#                # pT differential vn{SP}
#                for ipT in range(npT):
#                    pTlow = pTBoundary[ipT]; pThigh = pTBoundary[ipT+1]
#                    tempList = particleList[particleList[:,0]<pThigh,:]
#                    particleList_pTcut = tempList[tempList[:,0]>=pTlow,:]
#                    
#                    if(particleList_pTcut.size == 0):  # no particle in the bin
#                        pTmean[ipT] = (pTlow + pThigh)/2.
#                    else:
#                        NevpT[ipT] += 1
#                        pT = particleList_pTcut[:,0]
#                        phi = particleList_pTcut[:,1]
#                        pTmean[ipT] = mean(pT)
#                        Nparticle = len(pT)
#
#                        for weightType in weightTypes:
#                            if weightType == '1':
#                                weight = ones(len(pT))
#                                iweight = 0
#                            elif weightType == 'pT':
#                                weight = pT
#                                iweight = 1
#                            for iorder in range(1, norder + 1):
#                                QnA_X_local = QnA_X[iweight, iorder-1]
#                                QnA_Y_local = QnA_Y[iweight, iorder-1]
#
#                                # Qn vector for the interested particles
#                                Qn_X = sum(weight*cos(iorder*phi))/Nparticle
#                                Qn_Y = sum(weight*sin(iorder*phi))/Nparticle
#
#                                temp = Qn_X*QnA_X_local + Qn_Y*QnA_Y_local
#                                if weightType == '1':
#                                    vn_obs[iorder-1, ipT] += temp
#                                    vn_obs_sq[iorder-1, ipT] += temp*temp
#                                elif weightType == 'pT':
#                                    vn_obs_pTweight[iorder-1, ipT] += temp
#                                    vn_obs_pTweight_sq[iorder-1, ipT] += temp*temp
#        eps = 1e-15
#        vn_inte_obs = vn_inte_obs/(Nev_vninte + eps)
#        vn_inte_obs_err = sqrt(vn_inte_obs_sq/(Nev_vninte + eps) - (vn_inte_obs)**2.)/sqrt(Nev_vninte + eps)
#        vn_inte_obs_pTweight = vn_inte_obs_pTweight/(Nev_vninte + eps)
#        vn_inte_obs_pTweight_err = sqrt(vn_inte_obs_pTweight_sq/(Nev_vninte + eps) - (vn_inte_obs_pTweight)**2.)/sqrt(Nev_vninte + eps)
#        
#        vn_obs = vn_obs/(NevpT + eps)
#        vn_obs_err = sqrt(vn_obs_sq/(NevpT + eps) - (vn_obs)**2.)/sqrt(NevpT + eps)
#        vn_obs_pTweight = vn_obs_pTweight/(NevpT + eps)
#        vn_obs_pTweight_err = sqrt(vn_obs_pTweight_sq/(NevpT + eps) - (vn_obs_pTweight)**2.)/sqrt(NevpT + eps)
#
#        vninteSP = zeros([norder]); vninteSP_pTweight = zeros([norder])
#        vninteSP_err = zeros([norder]); vninteSP_pTweight_err = zeros([norder])
#        
#        vnSP = zeros([norder, npT]); vnSP_pTweight = zeros([norder, npT])
#        vnSP_err = zeros([norder, npT]); vnSP_pTweight_err = zeros([norder, npT])
#
#        for iorder in range(1, norder + 1):
#            resolutionFactor = self.db.executeSQLquery("select R_SP_sub from resolutionFactorR where weightType = '1' and n = %d" % (iorder)).fetchall()[0][0]
#            resolutionFactor_pTweight = self.db.executeSQLquery("select R_SP_sub from resolutionFactorR where weightType = 'pT' and n = %d" % (iorder)).fetchall()[0][0]
#
#            vninteSP[iorder-1] = vn_inte_obs[iorder-1]/resolutionFactor
#            vninteSP_err[iorder-1] = vn_inte_obs_err[iorder-1]/resolutionFactor
#            vninteSP_pTweight[iorder-1] = vn_inte_obs_pTweight[iorder-1]/resolutionFactor_pTweight
#            vninteSP_pTweight_err[iorder-1] = vn_inte_obs_pTweight_err[iorder-1]/resolutionFactor_pTweight
#
#            vnSP[iorder-1, :] = vn_obs[iorder-1, :]/resolutionFactor
#            vnSP_err[iorder-1, :] = vn_obs_err[iorder-1, :]/resolutionFactor
#            vnSP_pTweight[iorder-1, :] = vn_obs_pTweight[iorder-1, :]/resolutionFactor_pTweight
#            vnSP_pTweight_err[iorder-1, :] = vn_obs_pTweight_err[iorder-1, :]/resolutionFactor_pTweight
#        
#        return(vninteSP, vninteSP_err, vninteSP_pTweight, vninteSP_pTweight_err, pTmean, vnSP, vnSP_err, vnSP_pTweight, vnSP_pTweight_err)
#
#    def getScalarProductflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            retrieve nth order scalar product flow data from database for the given species 
#            of particles with given pT_range.
#            if no results is found, it will collect vn{SP}(pT) and store it into database
#        """
#        pid = self.pid_lookup[particleName]
#        collectFlag = False
#        if rapType == 'rapidity':
#            tableName = "diffvnSP"; tableName_vninte = "intevnSP"
#        elif rapType == 'pseudorapidity':
#            tableName = "diffvnSPeta"; tableName_vninte = "intevnSPeta"
#        if self.db.createTableIfNotExists(tableName, (("pid", "integer"), ("weightType", "text"), ("n", "integer"), ("pT", "real"), ("vn", "real"), ("vn_err", "real") )):
#            collectFlag = True
#        else:
#            vnSPdata = array(self.db.executeSQLquery("select pT, vn, vn_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName, pid, weightType, order)).fetchall())
#            if vnSPdata.size == 0:
#                collectFlag = True
#        if self.db.createTableIfNotExists(tableName_vninte, (("pid", "integer"), ("weightType", "text"), ("n", "integer"), ("vn", "real"), ("vn_err", "real") )):
#            collectFlag = True
#        else:
#            vninteSPdata = array(self.db.executeSQLquery("select vn, vn_err from %s where pid = %d and weightType = '%s' and n = %d" % (tableName_vninte, pid, weightType, order)).fetchall())
#            if vninteSPdata.size == 0:
#                collectFlag = True
#        if collectFlag:
#            vninteSP, vninteSP_err, vninteSP_pTweight, vninteSP_pTweight_err, pT, vnSP, vnSP_err, vnSP_pTweight, vnSP_pTweight_err = self.collectScalarProductflow(particleName = particleName, rapType = rapType) 
#            for iorder in range(len(vnSP[:,0])):
#                self.db.insertIntoTable(tableName_vninte, (pid, '1', iorder+1, vninteSP[iorder], vninteSP_err[iorder]))
#                self.db.insertIntoTable(tableName_vninte, (pid, 'pT', iorder+1, vninteSP_pTweight[iorder], vninteSP_pTweight_err[iorder]))
#                for ipT in range(len(pT)):
#                    self.db.insertIntoTable(tableName, (pid, '1', iorder+1, pT[ipT], vnSP[iorder, ipT], vnSP_err[iorder, ipT]))
#                    self.db.insertIntoTable(tableName, (pid, 'pT', iorder+1, pT[ipT], vnSP_pTweight[iorder, ipT], vnSP_pTweight_err[iorder, ipT]))
#            self.db._dbCon.commit()
#            if weightType == '1':
#                vninteSPdata = array([vninteSP[order-1], vninteSP_err[order-1]]).transpose()
#                vnSPdata = array([pT, vnSP[order-1,:], vnSP_err[order-1,:]]).transpose()
#            elif weightType == 'pT':
#                vninteSPdata = array([vninteSP_pTweight[order-1], vninteSP_pTweight_err[order-1]]).transpose()
#                vnSPdata = array([pT, vnSP_pTweight[order-1,:], vnSP_pTweight_err[order-1,:]]).transpose()
#            if vnSPdata.size == 0:
#                print("There is no record for different event plane flow vn for %s" % particleName)
#                return None
#      
#        #interpolate results to desired pT range
#        vnpTinterp = interp(pT_range, vnSPdata[:,0], vnSPdata[:,1])
#        vnpTinterp_err = interp(pT_range, vnSPdata[:,0], vnSPdata[:,2])
#        results = array([pT_range, vnpTinterp, vnpTinterp_err])
#        return(vninteSPdata, transpose(results))
#    
#    def getdiffScalarProductflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            return vn{SP}(pT)
#        """
#        vninteSP, vndiffSP = self.getScalarProductflow(particleName = particleName, order = order, weightType = weightType, pT_range = pT_range, rap_range = rap_range, rapType = rapType)
#        return(vndiffSP)
#    
#    def getinteScalarProductflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            return vn{SP}
#        """
#        vninteSP, vndiffSP = self.getScalarProductflow(particleName = particleName, order = order, weightType = weightType, pT_range = pT_range, rap_range = rap_range, rapType = rapType)
#        return(vninteSP)
#    
#    def collectTwoparticleCumulantflow(self, particleName = 'pion_p', pT_range = [0.0, 3.0], rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            collect nth order pT differential two particle cumulant harmonic flow, vn{2}(pT)
#            the two particles are taken from the same bin
#        """
#        print("Collecting differential two particle cumulant flow vn{2} for %s ..." % particleName)
#        norder = 6; npT = 30
#        weightTypes = ['1', 'pT']
#        pidString = self.getPidString(particleName)
#
#        pTBoundary = linspace(pT_range[0], pT_range[1], npT+1)
#        Nev_vn = 0.0
#        vn_inte_obs = zeros([norder]); vn_inte_obs_sq = zeros([norder])
#        vn_inte_obs_err = zeros([norder])
#
#        pTmean = zeros([npT]); NevpT = zeros(npT)
#        vn_obs = zeros([norder, npT]); vn_obs_sq = zeros([norder, npT])
#        vn_obs_err = zeros([norder, npT])
#
#        for hydroId in range(1, self.hydroNev+1):
#            UrQMDNev = self.db.executeSQLquery("select Number_of_UrQMDevents from UrQMD_NevList where hydroEventId = %d " % hydroId).fetchall()[0][0]
#            for UrQMDId in range(1, UrQMDNev+1):
#                print("processing event: (%d, %d) " % (hydroId, UrQMDId))
#                particleList = array(self.db.executeSQLquery("select pT, phi_p from particle_list where hydroEvent_id = %d and UrQMDEvent_id = %d and (%s) and (%g <= pT and pT < %g) and (%g <= %s and %s <= %g)" % (hydroId, UrQMDId, pidString, pT_range[0], pT_range[1], rap_range[0], rapType, rapType, rap_range[1])).fetchall())
#                
#                #calculate pT integrated v_n{2}
#                particleList_phi = particleList[:,1]
#                if(len(particleList_phi) > 1):
#                    Nev_vn += 1
#                    Nparticle = len(particleList_phi)
#                    Npairs = Nparticle*(Nparticle - 1)/2.
#                    for iorder in range(1, norder + 1):
#                        temp = 0.0
#                        for ipart in range(Nparticle-1):
#                            for iassopart in range(ipart+1, Nparticle):
#                                temp += cos(iorder*(particleList_phi[ipart] - particleList_phi[iassopart]))
#                        vn_inte_obs[iorder-1] += temp/Npairs
#                        vn_inte_obs_sq[iorder-1] += temp*temp/Npairs
#
#                #calculate pT differential v_n{2}(pT)
#                for ipT in range(npT):
#                    pTlow = pTBoundary[ipT]; pThigh = pTBoundary[ipT+1]
#                    tempList = particleList[particleList[:,0]<pThigh,:]
#                    particleList_pTcut = tempList[tempList[:,0]>=pTlow,:]
#                    
#                    if(particleList_pTcut.size == 0):
#                        pTmean[ipT] = (pTlow + pThigh)/2.
#                    elif(len(particleList[:,0]) == 1):
#                        pTmean[ipT] = particleList_pTcut[0,0]
#                    else:
#                        NevpT[ipT] += 1
#                        pT = particleList_pTcut[:,0]
#                        pTmean[ipT] = mean(pT)
#                        phi = particleList_pTcut[:,1]
#                        Nparticle = len(pT)
#                        Npairs = Nparticle*(Nparticle - 1)/2.
#
#                        for iorder in range(1, norder + 1):
#                            temp = 0.0
#                            for ipart in range(Nparticle-1):
#                                for iassopart in range(ipart+1, Nparticle):
#                                    temp += cos(iorder*(phi[ipart] - phi[iassopart]))
#
#                            vn_obs[iorder-1, ipT] += temp/Npairs
#                            vn_obs_sq[iorder-1, ipT] += temp*temp/Npairs
#        eps = 1e-15
#        vn_inte_obs = vn_inte_obs/(Nev_vn + eps)
#        vn_inte_obs_err = sqrt(vn_inte_obs_sq/(Nev_vn + eps) - (vn_inte_obs)**2.)/sqrt(Nev_vn + eps)
#
#        vn_obs = vn_obs/(NevpT + eps)
#        vn_obs_err = sqrt(vn_obs_sq/(NevpT + eps) - (vn_obs)**2.)/sqrt(NevpT + eps)
#        
#        vn_2_inte = zeros([norder]); vn_2_inte_err = zeros([norder])
#        vn_2 = zeros([norder, npT]); vn_2_err = zeros([norder, npT])
#        for iorder in range(norder):
#            if vn_inte_obs[iorder] > 0.0:
#                vn_2_inte[iorder] = sqrt(vn_inte_obs[iorder])
#                vn_2_inte_err[iorder] = vn_inte_obs_err[iorder]/(2.*sqrt(vn_inte_obs[iorder]))
#            for ipT in range(npT):
#                if vn_obs[iorder, ipT] > 0.0:
#                    vn_2[iorder, ipT] = sqrt(vn_obs[iorder, ipT])
#                    vn_2_err[iorder, ipT] = vn_obs_err[iorder, ipT]/(2.*sqrt(vn_obs[iorder, ipT]))
#        return(vn_2_inte, vn_2_inte_err, pTmean, vn_2, vn_2_err)
#
#    def getTwoparticlecumulantflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            retrieve nth order two particle cumulant flow data from database for the 
#            given species of particles with given pT_range.
#            if no results is found, it will collect vn{2}(pT) and store it into database
#        """
#        pid = self.pid_lookup[particleName]
#        collectFlag = False
#        if rapType == 'rapidity':
#            tableName = "diffvn2"; tableName_vninte = "intevn2"
#        elif rapType == 'pseudorapidity':
#            tableName = "diffvn2eta"; tableName_vninte = "intevn2eta"
#        if self.db.createTableIfNotExists(tableName, (("pid", "integer"), ("n", "integer"), ("pT", "real"), ("vn", "real"), ("vn_err", "real") )):
#            collectFlag = True
#        else:
#            vn2data = array(self.db.executeSQLquery("select pT, vn, vn_err from %s where pid = %d and n = %d" % (tableName, pid, order)).fetchall())
#            if vn2data.size == 0:
#                collectFlag = True
#        if self.db.createTableIfNotExists(tableName_vninte, (("pid", "integer"), ("n", "integer"), ("vn", "real"), ("vn_err", "real") )):
#            collectFlag = True
#        else:
#            vninte2data = array(self.db.executeSQLquery("select vn, vn_err from %s where pid = %d and n = %d" % (tableName_vninte, pid, order)).fetchall())
#            if vninte2data.size == 0:
#                collectFlag = True
#
#        if collectFlag:
#            vn2_inte, vn2_inte_err, pT, vn2, vn2_err = self.collectTwoparticleCumulantflow(particleName = particleName, rapType = rapType) 
#            for iorder in range(len(vn2[:,0])):
#                self.db.insertIntoTable(tableName_vninte, (pid, iorder+1, vn2_inte[iorder], vn2_inte_err[iorder]))
#                for ipT in range(len(pT)):
#                    self.db.insertIntoTable(tableName, (pid, iorder+1, pT[ipT], vn2[iorder, ipT], vn2_err[iorder, ipT]))
#            self.db._dbCon.commit()
#            vninte2data = array([vn2_inte[order-1], vn2_inte_err[order-1]]).transpose()
#            vn2data = array([pT, vn2[order-1,:], vn2_err[order-1,:]]).transpose()
#            if vn2data.size == 0:
#                print("There is no record for different event plane flow vn for %s" % particleName)
#                return None
#      
#        #interpolate results to desired pT range
#        vnpTinterp = interp(pT_range, vn2data[:,0], vn2data[:,1])
#        vnpTinterp_err = interp(pT_range, vn2data[:,0], vn2data[:,2])
#        results = array([pT_range, vnpTinterp, vnpTinterp_err])
#        return(vninte2data, transpose(results))
#    
#    def getdiffTwoparticlecumulantflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            return vn{2}(pT)
#        """
#        vn2inte, vn2diff = self.getTwoparticlecumulantflow(particleName = particleName, order = order, weightType = weightType, pT_range = pT_range, rap_range = rap_range, rapType = rapType)
#        return(vn2diff)
#    
#    def getinteTwoparticlecumulantflow(self, particleName = 'pion_p', order = 2, weightType = '1', pT_range = linspace(0.05, 2.5, 20), rap_range = [-0.5, 0.5], rapType = "rapidity"):
#        """
#            return vn{2}
#        """
#        vn2inte, vn2diff = self.getTwoparticlecumulantflow(particleName = particleName, order = order, weightType = weightType, pT_range = pT_range, rap_range = rap_range, rapType = rapType)
#        return(vn2inte)
#
#    def collectTwoparticleCorrelation(self, particleName = 'pion_p', pT_range = [1.0, 2.0]):
#        """
#            collect two particle correlation function C(\delta phi, \delta eta) from 
#            all the events within given pT range for a given particle species, "particleName"
#        """
#        pass

def printHelpMessageandQuit():
    print "Usage : "
    print "AnalyzedDataReader.py databaseName"
    exit(0)

if __name__ == "__main__":
    if len(argv) < 2:
        printHelpMessageandQuit()
    test = AnalyzedDataReader(str(argv[1]))
    print(test.get_particle_spectra('pion_p', pT_range=linspace(0.1, 2.5, 20), rap_type = 'pseudorapidity'))
    print(test.get_particle_yield_vs_rap('pion_p', rap_type = 'rapidity', rap_range=linspace(-1.0, 1.0, 15)))
    print(test.get_particle_yield('pion_p', rap_type = 'rapidity', rap_range=(-0.5, 0.5)))
    print(test.get_particle_yield_vs_spatial_variable('pion_p', 'tau', 
          linspace(0.6, 10, 50), rap_type = 'rapidity'))
   # print(test.get_avg_diffvn_flow('pion_p', 2, psi_r = 0., 
   #       pT_range = linspace(0.0, 2.0, 21)))
   # print(test.get_avg_intevn_flow('pion_p', 2, psi_r = 0., 
   #       pT_range = (0.3, 3.0)))
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
    #print(test.getdiffTwoparticlecumulantflow(particleName = 'pion_p'))
    #print(test.getdiffTwoparticlecumulantflow(particleName = 'kaon_p'))
    #print(test.getdiffTwoparticlecumulantflow(particleName = 'proton'))
    #test.collectEventplaneflow('charged', 2)
    #test.collectTwoparticleCorrelation()

