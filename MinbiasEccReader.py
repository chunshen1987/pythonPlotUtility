#! /usr/bin/env python

import sys
from os import path
from DBR import SqliteDB
from numpy import *
from CSplottools import getBinnedAveragedDatawithErrorbars

# define colors
purple = "\033[95m"
green = "\033[92m"
red = "\033[91m"
normal = "\033[0m"


class MinbiasEccReader(object):
    """
        This class contains functions to perform statistical analysis on
        initial eccentricities from minimum bias events generated from superMC
    """

    def __init__(self, database_name):
        """
            Register a sqlite database
        """
        # setup database
        database = None
        if isinstance(database_name, str):
            if path.exists(database_name):
                database = SqliteDB(database_name)
                self.db_name = database_name
            else:
                raise ValueError(
                    "EbeDBReader.__init__: the input argument must be an "
                    "existing database file. %s can not be found"
                    % (red + database_name + normal))
        if isinstance(database, SqliteDB):
            self.db = database
        else:
            raise TypeError(
                "EbeDBReader.__init__: the input argument must be a string "
                "or a SqliteDB database.")

        # define centrality boundaries
        #self.centrality_boundaries = [
        #    (0, 0.2), (0, 1), (0, 5), (5, 10), (10, 20), (20, 30), (30, 40),
        #    (40, 50), (50, 60), (60, 70), (70, 80)]

        self.centrality_boundaries = [
            (0, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, 30), 
            (30, 35), (35, 40), (40, 45), (45, 50), (50, 55), (55, 60), 
            (60, 65), (65, 70), (70, 75), (75, 80)]

        # get total number of events
        self.nev = self.get_number_of_events()

    def get_number_of_events(self):
        nev = self.db.executeSQLquery(
            "select count(*) from collisionParameters").fetchall()[0][0]
        return nev

    def cut_centralities_with_ecc_statistics(
            self, cut_type, multiplicity_factor=1.0):
        """
            this function cut the centralities and also output event averaged
            ecc_n with statistical error
        """
        print ('cutting centralities from database %s according to %s ....'
               % (purple + self.db_name + normal, green + cut_type + normal))
        centrality_output = open('centralityCut_%s.dat' % cut_type, 'w')
        eccn_stat_ed_output = open('eccnStatistics_ed_%s.dat' % cut_type, 'w')
        eccn_stat_sd_output = open('eccnStatistics_sd_%s.dat' % cut_type, 'w')
        nevent = self.nev

        for icen in range(len(self.centrality_boundaries)):
            lowerbound = self.centrality_boundaries[icen][0]
            upperbound = self.centrality_boundaries[icen][1]
            nsample = int(nevent * (upperbound - lowerbound) / 100) - 1
            noffset = int(nevent * lowerbound / 100)
            if cut_type == 'b':
                fetched_data = array(self.db.executeSQLquery(
                    "select Npart, b, total_entropy from collisionParameters "
                    "order by %s limit %d offset %d"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            else:
                fetched_data = array(self.db.executeSQLquery(
                    "select Npart, b, total_entropy from collisionParameters "
                    "order by -%s limit %d offset %d"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            cen_central = (upperbound + lowerbound) / 2.
            npart_mean = mean(fetched_data[:, 0])
            npart_min = min(fetched_data[:, 0])
            npart_max = max(fetched_data[:, 0])
            bmin = min(fetched_data[:, 1])
            bmax = max(fetched_data[:, 1])
            dsdymin = min(fetched_data[:, 2]) / multiplicity_factor
            dsdymax = max(fetched_data[:, 2]) / multiplicity_factor
            centrality_output.write(
                "%6.4f  %d  %d  %d  %18.8e  %18.8e  %18.8e  %18.8e \n"
                % (cen_central, npart_mean, npart_min, npart_max, dsdymin,
                   dsdymax, bmin, bmax)
            )
            if cut_type == 'b':
                fetched_data = array(self.db.executeSQLquery(
                    "select ecc_id, n, ecc_real, ecc_imag from eccentricities "
                    "where event_id in (select event_id from "
                    "collisionParameters order by collisionParameters.%s "
                    "limit %d offset %d)"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            else:
                fetched_data = array(self.db.executeSQLquery(
                    "select ecc_id, n, ecc_real, ecc_imag from eccentricities "
                    "where event_id in (select event_id from "
                    "collisionParameters order by -collisionParameters.%s "
                    "limit %d offset %d)"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            for ecc_type in range(1, 3):
                tempidx = (fetched_data[:, 0] == ecc_type)
                temp_data = fetched_data[tempidx, :]
                ecc_output = []
                for iorder in range(1, 10):
                    idx = (temp_data[:, 1] == iorder)
                    eccn2 = sqrt(mean(temp_data[idx, 2] ** 2
                                      + temp_data[idx, 3] ** 2))
                    eccn2err = (std(temp_data[idx, 2] ** 2
                                    + temp_data[idx, 3] ** 2)
                                / (2. * eccn2) / sqrt(nsample))
                    ecc_output += [eccn2, eccn2err]
                if ecc_type == 1:
                    eccn_stat_sd_output.write("%6.4f  " % cen_central)
                    for tempecc in ecc_output:
                        eccn_stat_sd_output.write("%18.8e  " % tempecc)
                    eccn_stat_sd_output.write("\n")
                elif ecc_type == 2:
                    eccn_stat_ed_output.write("%6.4f  " % cen_central)
                    for tempecc in ecc_output:
                        eccn_stat_ed_output.write("%18.8e  " % tempecc)
                    eccn_stat_ed_output.write("\n")

        centrality_output.close()
        eccn_stat_sd_output.close()
        eccn_stat_ed_output.close()

    def centrality_cut_with_avg_collisional_parameters(self, cut_type):
        """
            this function cut the centralities and output event averaged
            collisional parameters
        """
        print ('cutting centralities from database %s according to %s ....'
               % (purple + self.db_name + normal, green + cut_type + normal))
        centrality_output = open(
            'avgcollisionParameters_centralityCut_%s.dat' % cut_type, 'w')
        centrality_output.write(
            "#centrality <Npart> Npart_err <b> b_err <dS/dy> dS/dy_err "
            "<Ncoll> Ncoll_err \n")
        nevent = self.nev

        for icen in range(len(self.centrality_boundaries)):
            lowerbound = self.centrality_boundaries[icen][0]
            upperbound = self.centrality_boundaries[icen][1]
            nsample = int(nevent * (upperbound - lowerbound) / 100) - 1
            noffset = int(nevent * lowerbound / 100)
            if cut_type == 'b':
                fetched_data = array(self.db.executeSQLquery(
                    "select Npart, b, total_entropy, Ncoll from "
                    "collisionParameters order by %s limit %d offset %d"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            else:
                fetched_data = array(self.db.executeSQLquery(
                    "select Npart, b, total_entropy, Ncoll from "
                    "collisionParameters order by -%s limit %d offset %d"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            cen_central = (upperbound + lowerbound) / 2.
            nev = len(fetched_data[:, 0])
            npart_mean = mean(fetched_data[:, 0])
            b_mean = mean(fetched_data[:, 1])
            dsdy_mean = mean(fetched_data[:, 2])
            ncoll_mean = mean(fetched_data[:, 3])
            npart_err = std(fetched_data[:, 0]) / sqrt(nev)
            b_err = std(fetched_data[:, 1]) / sqrt(nev)
            dsdy_err = std(fetched_data[:, 2]) / sqrt(nev)
            ncoll_err = std(fetched_data[:, 3]) / sqrt(nev)
            centrality_output.write(
                "%6.4e  %18.8e  %18.8e  %18.8e  %18.8e  %18.8e  %18.8e  "
                "%18.8e  %18.8e \n"
                % (cen_central, npart_mean, npart_err, b_mean, b_err,
                   dsdy_mean, dsdy_err, ncoll_mean, ncoll_err)
            )
        centrality_output.close()

    def centrality_cut_with_avg_collisional_parameters_latex(self, cut_type):
        """
            this function cut the centralities and output event averaged
            collisional parameters
        """
        print ('cutting centralities from database %s according to %s ....'
               % (purple + self.db_name + normal, green + cut_type + normal))
        centrality_output = open(
            'avgcollisionParameters_centralityCut_latex_%s.dat' % cut_type,
            'w')
        centralitycut_output = open(
            'centralityCut_latex_%s.dat' % cut_type, 'w')
        centrality_output.write(
            "centrality (\%) & $\langle N_\mathrm{part} \\rangle$ & "
            "$\langle b \\rangle$ (fm) & $\langle dS/dy \\rangle$  & "
            "$\langle N_\mathrm{coll} \\rangle$ \\\\ \hline \n")
        centralitycut_output.write(
            "centrality (\%) & $N_\mathrm{min}$ & $N_\mathrm{max}$ & "
            "$dS/dy_\mathrm{min}$ & $dS/dy_\mathrm{max}$ & "
            "$b_\mathrm{min}$ (fm) & $b_\mathrm{max}$ (fm) \\\\ \hline \n")
        nevent = self.nev

        for icen in range(len(self.centrality_boundaries)):
            lowerbound = self.centrality_boundaries[icen][0]
            upperbound = self.centrality_boundaries[icen][1]
            nsample = int(nevent * (upperbound - lowerbound) / 100) - 1
            noffset = int(nevent * lowerbound / 100)
            if cut_type == 'b':
                fetched_data = array(self.db.executeSQLquery(
                    "select Npart, b, total_entropy, Ncoll from "
                    "collisionParameters order by %s limit %d offset %d"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            else:
                fetched_data = array(self.db.executeSQLquery(
                    "select Npart, b, total_entropy, Ncoll from "
                    "collisionParameters order by -%s limit %d offset %d"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            nev = len(fetched_data[:, 0])
            npart_mean = mean(fetched_data[:, 0])
            npart_min = min(fetched_data[:, 0])
            npart_max = max(fetched_data[:, 0])
            b_mean = mean(fetched_data[:, 1])
            bmin = min(fetched_data[:, 1])
            bmax = max(fetched_data[:, 1])
            dsdy_mean = mean(fetched_data[:, 2])
            dsdymin = min(fetched_data[:, 2])
            dsdymax = max(fetched_data[:, 2])
            ncoll_mean = mean(fetched_data[:, 3])
            npart_err = std(fetched_data[:, 0]) / sqrt(nev)
            b_err = std(fetched_data[:, 1]) / sqrt(nev)
            dsdy_err = std(fetched_data[:, 2]) / sqrt(nev)
            ncoll_err = std(fetched_data[:, 3]) / sqrt(nev)
            centrality_output.write(
                "%g-%g & %10.2f $\pm$ %10.2f & %10.2f $\pm$ %10.2f & "
                "%10.2f $\pm$ %10.2f & %10.2f $\pm$ %10.2f \\\\ \hline \n"
                % (lowerbound, upperbound, npart_mean, npart_err,
                   b_mean, b_err, dsdy_mean, dsdy_err, ncoll_mean, ncoll_err)
            )
            centralitycut_output.write(
                "%g-%g &  %d & %d & %10.2f & %10.2f & %10.2f & %10.2f "
                "\\\\ \hline \n"
                % (lowerbound, upperbound, npart_min, npart_max, dsdymin,
                   dsdymax, bmin, bmax)
            )
        centrality_output.close()
        centralitycut_output.close()

    def get_distribution(self, dis_type='total_entropy', nbin=30,
                         cut_type='total_entropy', centrality_bound=(0, 100)):
        """
            this function bin and output distribution of a given quantity
            in a given centrality bin output format:
            X, P(X), P(X)_err, dP(X)/dX, dP(X)/dX_err
        """
        nevent = self.nev
        lowerbound = centrality_bound[0]
        upperbound = centrality_bound[1]
        nsample = int(nevent * (upperbound - lowerbound) / 100) - 1
        noffset = int(nevent * lowerbound / 100)
        if dis_type in ['Npart', 'Ncoll', 'b', 'total_entropy']:
            if cut_type in ['Npart', 'Ncoll', 'total_entropy']:
                fetched_data = array(self.db.executeSQLquery(
                    "select %s from collisionParameters order by -%s limit "
                    "%d offset %d"
                    % (dis_type, cut_type, nsample, noffset)
                ).fetchall())
            elif cut_type in ['b']:
                fetched_data = array(self.db.executeSQLquery(
                    "select %s from collisionParameters order by %s limit "
                    "%d offset %d"
                    % (dis_type, cut_type, nsample, noffset)
                ).fetchall())
        elif 'ecc' in dis_type:
            temp = dis_type.split('_')
            eccorder = int(temp[1])
            if cut_type in ['Npart', 'Ncoll', 'total_entropy']:
                temp_data = array(self.db.executeSQLquery(
                    "select ecc_real, ecc_imag from eccentricities where "
                    "ecc_id = 2 and n = %d and event_id in (select "
                    "event_id from collisionParameters order by "
                    "-collisionParameters.%s limit %d offset %d)"
                    % (eccorder, cut_type, nsample, noffset)
                ).fetchall())
            elif cut_type in ['b']:
                temp_data = array(self.db.executeSQLquery(
                    "select ecc_real, ecc_imag from eccentricities where "
                    "ecc_id = 2 and n = %d and event_id in (select "
                    "event_id from collisionParameters order by "
                    "collisionParameters.%s limit %d offset %d)"
                    % (eccorder, cut_type, nsample, noffset)).fetchall())
            fetched_data = sqrt(temp_data[:, 0] ** 2 + temp_data[:, 1] ** 2)
        elif 'deformed' in dis_type:
            temp = dis_type.split('_')
            dis_quantity = temp[1]
            if cut_type in ['Npart', 'Ncoll', 'total_entropy']:
                fetched_data = array(self.db.executeSQLquery(
                    "select %s from deformationParameters where event_id "
                    "in (select event_id from collisionParameters order "
                    "by -collisionParameters.%s limit %d offset %d)"
                    % (dis_quantity, cut_type, nsample, noffset)
                ).fetchall())
            elif cut_type in ['b']:
                fetched_data = array(self.db.executeSQLquery(
                    "select %s from deformationParameters where event_id "
                    "in (select event_id from collisionParameters order "
                    "by collisionParameters.%s limit %d offset %d)"
                    % (dis_quantity, cut_type, nsample, noffset)
                ).fetchall())
        else:
            raise TypeError("unrecognized distriubtion type.")

        binned_data, binned_data_err = getBinnedAveragedDatawithErrorbars(
            fetched_data, nbin)
        dis_output = open('%s_distribution_C%g-%g_%s.dat'
                          % (dis_type, centrality_bound[0],
                             centrality_bound[1], cut_type), 'w')
        for i in range(nbin):
            dis_output.write(
                "%18.8e  %18.8e  %18.8e  %18.8e   %18.8e\n"
                % (binned_data[i, 0], binned_data[i, 1], binned_data_err[i, 1],
                   binned_data[i, 2], binned_data_err[i, 2])
            )
        dis_output.close()
        return array([binned_data[:, 0], binned_data[:, 1],
                      binned_data_err[:, 1], binned_data[:, 2],
                      binned_data_err[:, 2]]).transpose()

    def output_centrality_cut_table_iebe_package(self, cut_type):
        """
        This function output centrality cut table for iEBE package
        :param cut_type: centrality cut type
        """
        print ('cutting centralities from database %s according to %s ....'
               % (purple + self.db_name + normal, green + cut_type + normal))
        centrality_output = (
            open('iebe_centralityCut_%s_%s.dat'
                 % (cut_type, self.db_name.split('/')[-1].split('.db')[0]),
                 'w'))
        if cut_type == 'total_entropy':
            centrality_output.write(
                "#centrality dS/dy Npart_min Npart_max b_min(fm) b_max(fm)\n")
            centrality_output.write(
                "%6.4e  %18.8e  %18.8e  %18.8e  %18.8e  %18.8e \n"
                % (0.0, 1000000, 500, 500, 0.0, 0.0))
        elif cut_type == 'Npart':
            centrality_output.write("#centrality Npart b_min(fm) b_max(fm)\n")
            centrality_output.write("%6.4e  %18.8e  %18.8e  %18.8e \n"
                % (0.0, 500, 0.0, 0.0))
        elif cut_type == 'b':
            centrality_output.write("#centrality b(fm)\n")
            centrality_output.write("%6.4e %18.8e \n" % (0.0, 0.0))
        else:
            print "invalid cutType", red + cut_type + normal

        nevent = self.nev

        centrality_bound = [x*0.1 for x in range(0, 10)]+range(1, 101)

        for icen in range(1,len(centrality_bound)):
            lowerbound = centrality_bound[icen-1]
            upperbound = centrality_bound[icen]
            nsample = int(nevent * (upperbound - lowerbound) / 100) - 1
            noffset = int(nevent * lowerbound / 100)
            if cut_type == 'b':
                fetched_data = array(self.db.executeSQLquery(
                    "select Npart, b, total_entropy from "
                    "collisionParameters order by %s limit %d offset %d"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            else:
                fetched_data = array(self.db.executeSQLquery(
                    "select Npart, b, total_entropy from "
                    "collisionParameters order by -%s limit %d offset %d"
                    % (cut_type, nsample, noffset)
                ).fetchall())
            npart_min = min(fetched_data[:, 0])
            npart_max = max(fetched_data[:, 0])
            b_min = min(fetched_data[:, 1])
            b_max = max(fetched_data[:, 1])
            dsdy_min = min(fetched_data[:, 2])
            if cut_type == 'total_entropy':
                centrality_output.write(
                    "%6.4e  %18.8e  %18.8e  %18.8e  %18.8e  %18.8e \n"
                    % (upperbound, dsdy_min, npart_min, npart_max,
                       b_min, b_max)
                )
            elif cut_type == 'Npart':
                centrality_output.write("%6.4e  %18.8e  %18.8e  %18.8e \n"
                    % (upperbound, npart_min, b_min, b_max)
                )
            elif cut_type == 'b':
                centrality_output.write("%6.4e  %18.8e \n"
                    % (upperbound, b_max)
                )
        centrality_output.close()


if __name__ == "__main__":
    reader = MinbiasEccReader(str(sys.argv[1]))
    reader.cut_centralities_with_ecc_statistics('total_entropy', 1.0)
