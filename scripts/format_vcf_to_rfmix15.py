#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
import vcfpy
import re, csv
import os, errno
import sys, getopt
import os.path

def main(argv):
    inputfile = ''
    outfile = ''

    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print(os.path.basename(__file__) + ' -i <phased vcf to format> -o <outputfile>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(os.path.basename(__file__) + ' -i <phased vcf to format> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outfile = arg

    if inputfile=='' or outfile=='':
        print("ERROR: You must provide input and output filenames")
        sys.exit(2)

    j = 0
    fo = open(outfile,'w')
    print('Processing file "' + inputfile + '"')
    with open(inputfile,'r') as f:
        rows = csv.reader(f, delimiter='\t')
        for row in rows:
            #pos = 0
            if not row[0].startswith('#'):
                chr = str(row[0])
                pos = str(row[1])
                ref = str(row[3])
                alt = str(row[4])
                snp = chr+':'+pos+':'+ref+':'+alt
                hap = ''
                for i in range(9,len(row)):
                    gtype = str(row[i]).replace('/','|')
                    gtypeg = re.match('(^[0-1]\|[0-1])',gtype)
                    x = gtypeg.groups()[0].split('|')
                    hap += str(x[0]) + str(x[1])

                fo.write(hap + '\n')

                if j > 0 and j % 1000 == 0:
                    print('Reading ' + str(j) + ' records and writing haplotypes')
                j += 1

        print(str(j) + ' positions were read')
    f.close()
    fo.close()

if __name__ == "__main__":
    main(sys.argv[1:])
