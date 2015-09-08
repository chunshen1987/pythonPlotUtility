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
    print "minbias_ALICE_standardcandle.py databaseName -type cutType"
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
        elif option == '-h' or option == '-help':
            print_help_message_and_quit()
        else:
            print argv[0], ': invalid option', option
            print_help_message_and_quit()
    reader = MinbiasEccReader(db_name)
    reader.ALICE_standard_candles_vs_centrality(3, 2, cut_type)
    reader.ALICE_standard_candles_vs_centrality(4, 2, cut_type)
