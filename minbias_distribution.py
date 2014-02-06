#! /usr/bin/env python

from sys import argv, exit
from MinbiasEccReader import MinbiasEccReader

#define colors
purple = "\033[95m"
green = "\033[92m"
red = "\033[91m"
normal = "\033[0m"


def print_help_message_and_quit():
    print("minbias_distribution outputs" + purple +
          " X, P(X), P(X)_err, dP(X)/dX, dP(X)/dX_err" + normal)
    print "Usage : "
    print("minbias_distribution %s -disType disType -nbin nbin "
          "-cutType cutType -cen centralityBoundary"
          % (purple + 'databaseName' + normal))
    print "Usage of minbias_CentralityCut command line arguments: "
    print("-disType   the type of quantity of the output distribution: "
          + green
          + " Npart, Ncoll, b, total_entropy, ecc_n, deformed_cosTheta1"
          + normal)
    print "-nbin      set number of bins "
    print("-cutType   the type of quantity used to cut centrality: "
          + green + " Npart, Ncoll, b, total_entropy" + normal)
    print("-cen       the boundary of centrality cut: " + purple +
          " e.g. 0-5" + normal)
    print "-h | -help   This message"
    exit(0)


if __name__ == "__main__":
    # set default values
    cut_type = 'total_entropy'
    cen = [0, 100]
    nbin = 30
    if len(argv) <= 2:
        print_help_message_and_quit()
    db_name = str(argv[1])
    while len(argv) > 2:
        option = argv[2]
        del argv[2]
        if option == '-disType':
            dis_type = str(argv[2])
            del argv[2]
            if(not dis_type in ['Npart', 'Ncoll', 'b', 'total_entropy']
               and not 'ecc' in dis_type and not 'deformed' in dis_type):
                print argv[0], ': invalid disType', red + dis_type + normal
                print_help_message_and_quit()
        elif option == '-cutType':
            cut_type = str(argv[2])
            del argv[2]
            if not cut_type in ['Npart', 'Ncoll', 'b', 'total_entropy']:
                print argv[0], ': invalid cutType', red + cut_type + normal
                print_help_message_and_quit()
        elif option == '-cen':
            cen = str(argv[2]).split('-')
            del argv[2]
            cen = [float(x) for x in cen]
        elif option == '-nbin':
            nbin = int(argv[2])
            del argv[2]
        elif option == '-h' or option == '-help':
            print_help_message_and_quit()
        else:
            print argv[0], ': invalid option', option
            print_help_message_and_quit()
    reader = MinbiasEccReader(db_name)
    print(reader.get_distribution(dis_type, nbin, cut_type, centrality_bound=cen))
