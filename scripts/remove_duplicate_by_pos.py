#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
import vcfpy
import re, csv
import os, errno
import sys, getopt
from shutil import copyfile

def main(argv):
    inputfile = ''
    outputfile = ''

    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print(os.path.basename(__file__) + ' -i <file.vcf> -o <outputfile>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(os.path.basename(__file__) + ' -i <file.vcf>  -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    if inputfile=='' or outputfile=='':
        print("ERROR: You must provide input and output filenames")
        sys.exit(2)

    #VCF = vcfpy.Reader.from_path(inputfile) # to check VCF integrity

    RECORD = {}
    s = 0
    print('Processing file "' + inputfile + '"')
    with open(inputfile,'r') as f:
        rows = csv.reader(f, delimiter='\t')
        i = 0
        for row in rows:
            pos = 0
            if not row[0].startswith('#'):
                chr = row[0]
                pos = row[1]
                if pos in RECORD:
                    RECORD[pos] += 1
                else:
                    RECORD[pos] = 1
                s += 1
                if i > 0 and i % 10000 == 0:
                    print('Reading ' + str(i) + ' records')
                i += 1
        print('Reading a total of ' + str(i) + ' records')
    f.close()

    if s == len(RECORD):
        copyfile(inputfile, outputfile)
        print('Writing a total of ' + str(i) + ' records')
        exit(0)


    fo = open(outputfile,'w+')
    fcsv = csv.writer(fo, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='', escapechar='\\')
    with open(inputfile,'r') as f:
        rows = csv.reader(f, delimiter='\t')
        i = 0
        for row in rows:
            pos = 0
            if not row[0].startswith('#'):
                chr = row[0]
                pos = row[1]
                if RECORD[pos] == 1:
                    fcsv.writerow(row)
                    if i % 10000 == 0 and i > 0:
                        print('Writing ' + str(i) + ' records')
                    i += 1
                else:
                    print('WARNING: both positions ' + str(pos) + ' of duplicate SNP were removed')
            else:
                fcsv.writerow(row)

        print('Writing a total of ' + str(i) + ' records')
    fo.close()
    f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
