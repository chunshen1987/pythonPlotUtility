#! /usr/bin/env python

from sys import argv, exit
from MinbiasEccReader import MinbiasEccReader

#define colors
purple = "\033[95m"
green = "\033[92m"
red = "\033[91m"
normal = "\033[0m"


def print_help_message_and_quit():
    print "Usage : "
    print "minbias_CentralityCut databaseName -type cutType"
    print "Usage of minbias_CentralityCut command line arguments: "
    print("-type   the type of quantity used to cut centrality: " + green
          + " Npart, Ncoll, b, total_entropy" + normal)
    print "-h | -help   This message"
    exit(0)


if __name__ == "__main__":
    multiplicity_factor = 1.0
    if len(argv) <= 2:
        print_help_message_and_quit()
    db_name = str(argv[1])
    while len(argv) > 2:
        option = argv[2]
        del argv[2]
        if option == '-type':
            cut_type = str(argv[2])
            del argv[2]
            if not cut_type in ['Npart', 'Ncoll', 'b', 'total_entropy']:
                print argv[0], ": invalid cutType", red + cut_type + normal
                print_help_message_and_quit()
        elif option == '-mult':
            multiplicity_factor = float(argv[2])
            del argv[2]
        elif option == '-h' or option == '-help':
            print_help_message_and_quit()
        else:
            print argv[0], ': invalid option', option
            print_help_message_and_quit()
    reader = MinbiasEccReader(db_name)
    #reader.centrality_cut_with_avg_collisional_parameters_latex(cut_type)
    #reader.cut_centralities_with_ecc_statistics(cut_type, multiplicity_factor)
    reader.centrality_cut_with_avg_collisional_parameters(cut_type)
    #reader.output_centrality_cut_table_iebe_package(cut_type)
