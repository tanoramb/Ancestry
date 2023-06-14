#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Tomas Norambuena A.
# email: tanoramb@gmail.com

import pysam
import vcfpy
import re, csv
import os, errno
import sys, getopt
import os.path

def compl(b):
	if str(b).upper() == 'A':
		return 'T'
	elif str(b).upper() == 'C':
		return 'G'
	elif str(b).upper() == 'G':
		return 'C'
	elif str(b).upper() == 'T':
		return 'A'
	else:
		return '.'

def main(argv):
	inputfiles = ''
	outputfile = ''
	refsnps = ''
	printintfiles = False

	try:
		opts, args = getopt.getopt(argv,"hi:o:ps:",["ifiles=","ofile=","--printintfiles","--refsnps"])
	except getopt.GetoptError:
		print(os.path.basename(__file__) + '-s <chr_pos_ref_alt.txt> -i <file1.vcf,file2.vcf,file3.vcf,...,fileN.vcf> -o <outputfile>')
		sys.exit(2)
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print(os.path.basename(__file__) + '-s <chr_pos_ref_alt.txt> -i <file1.vcf,file2.vcf,file3.vcf,...,fileN.vcf>  -o <outputfile>')
			sys.exit()
		elif opt in ("-i", "--ifiles"):
			inputfiles = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
		elif opt in ("-s", "--refsnps"):
			refsnps = arg
		elif opt in ("-p", "--printintfiles"):
			printintfiles = True
	
	if inputfiles=='' or outputfile=='' or refsnps=='':
		print("ERROR: You must provide input and output, as well as refsnps filenames")
		sys.exit(2)
	
	#print('Input file is "', inputfile)
	#print('Output file is "', outputfile)

	POS = {}
	POS['REF'] = []
	print('Processing REFERENCE SNPs file "' + refsnps + '"')
	
	with open(refsnps,'r') as f:
		rows = csv.reader(f, delimiter='\t')
		for row in rows:
			#pos = 0
			if not row[0].startswith('#'):
				chr = str(row[0])
				pos = str(row[1])
				ref = str(row[2])
				alt = str(row[3])
				
				snp = chr+'_'+pos+ref+alt
				if alt == ".":
					print('WARNING: Reference SNP '+snp+' is monomorphic. Monomorphic SNPs are allowed in files to be normalized not in reference SNPs.')
				else:
					POS['REF'].append(snp)
		print(str(len(POS['REF'])) + ' reference positions were read')
	f.close()

	FILES=inputfiles.split(",")
	print(FILES)
	i = 0
	for inputfile in FILES:
		#VCF = vcfpy.Reader.from_path(inputfile) # verify VCF integrity
		if not os.path.isfile(inputfile):
			print('ERROR: File "' + inputfile + '" does not exist')
			exit(1)
	
	for inputfile in FILES:
		j = 0
		print('Processing file "' + inputfile + '"')
		POS[inputfile] = [0] * len(POS['REF'])
		
		with open(inputfile,'r') as f:
			rows = csv.reader(f, delimiter='\t')
			for row in rows:
				#pos = 0
				if not row[0].startswith('#'):
					chr = str(row[0])
					pos = str(row[1])
					ref = str(row[3])
					alt = str(row[4])
					if len(ref) == 1 and len(alt) == 1:
						# POLYMORPHIC
						if alt != ".":
							snp = chr+'_'+pos+ref+alt
							_snp = snp
							if snp in POS['REF']:
								index = POS['REF'].index(snp)
								POS[inputfile][index] = 1
							else:
								snp = chr+'_'+pos+alt+ref
								if snp in POS['REF']:
									index = POS['REF'].index(snp)
									POS[inputfile][index] = 1
									print('WARNING: ALT<->REF '+_snp+' (read) '+snp+' (written)') 
						# MONOMORPHIC
						else:
							snp = chr+'_'+pos+ref+alt
							_snp = snp
							snpl = [chr+'_'+pos+ref+'A', chr+'_'+pos+ref+'C', chr+'_'+pos+ref+'T', chr+'_'+pos+ref+'G']
							hits = [x in POS['REF'] for x in snpl]
							if any(hits):
								index = hits.index(True)
								snp = snpl[index]
								index = POS['REF'].index(snp)
								POS[inputfile][index] = 1
								print('WARNING: MONOMORPH '+_snp+' (read) '+snp+' (written)') 
							else:
								snpl = [chr+'_'+pos+'A'+ref, chr+'_'+pos+'C'+ref, chr+'_'+pos+'T'+ref, chr+'_'+pos+'G'+ref]
								hits = [x in POS['REF'] for x in snpl]
								if any(hits):
									index = hits.index(True)
									snp = snpl[index]
									index = POS['REF'].index(snp)
									POS[inputfile][index] = 1
									print('WARNING: MONOMORPH & ALT<->REF '+_snp+' (read) '+snp+' (written)') 
						
					if j > 0 and j % 10000 == 0:
						print('Reading ' + str(j) + ' records')
					j += 1
					
			print(str(j) + ' positions were read')
		f.close()
	
	HIT = {}
	i = 0
	for pos in POS['REF']:
		HIT[i] = [0] * (len(FILES) - 0)
		i += 1
	
	for filename, hits in POS.items():
		if filename != "REF":
			index = FILES.index(filename) - 0
			i = 0
			for hit in hits:
				HIT[i][index] = hit
				i += 1

	f = open(outputfile,"w+")
	n = len(FILES) - 0
	HITPOS = []
	j = 0
	for i, hits in HIT.items():
		if sum(hits) == n:
			HITPOS.append(POS['REF'][i])
			f.write(str(POS['REF'][i])[:-2] + '\n')
			j += 1
	
	f.close()
	print(str(j) + ' positions are common')

	if printintfiles:
		print('Printing Intersected files')
		for inputfile in FILES:
			outputfile = 'intersected.norm.' + inputfile.split("/")[-1]
			print('Printing file "' + outputfile + '"')
			fo = open(outputfile,'w+', newline='\n', encoding='utf-8')
			fcsv = csv.writer(fo, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='', escapechar='\\')
			with open(inputfile,'r') as f:
				rows = csv.reader(f, delimiter='\t')
				i = 0
				for row in rows:
					#pos = 0
					if not row[0].startswith('#'):
						printi = False
						chr = str(row[0])
						pos = str(row[1])
						ref = str(row[3])
						alt = str(row[4])
						if len(ref) == 1 and len(alt) == 1:
							if alt != ".":
								snp = chr+'_'+pos+ref+alt
								if snp in HITPOS:
									for k in range(9,len(row)):
										#row[k] = str(row[k]).replace('|','/')
										if str(row[k]) == '1/0':
											row[k] = '0/1'
									fcsv.writerow(row)
									i += 1
									printi = True
								else:
									snp = chr+'_'+pos+alt+ref
									if snp in HITPOS:
										row[3] = alt
										row[4] = ref
										for k in range(9,len(row)):
											#row[k] = str(row[k]).replace('|','/')
											if str(row[k]) == '0/0':
												row[k] = '1/1'
											elif str(row[k]) == '1/1':
												row[k] = '0/0'
											elif str(row[k]) == '1/0':
												row[k] = '0/1'
											elif str(row[k]) == '0|0':
												row[k] = '1|1'
											elif str(row[k]) == '1|1':
												row[k] = '0|0'
											elif str(row[k]) == '1|0':
												row[k] = '0|1'
											elif str(row[k]) == '0|1':
												row[k] = '1|0'
										fcsv.writerow(row)
										#print('WARNING: ALT<->REF '+_snp+' '+snp)
										i += 1
										printi = True
								
							else:
								snp = chr+'_'+pos+ref+alt
								_snp = snp
								snpl = [chr+'_'+pos+ref+'A', chr+'_'+pos+ref+'C', chr+'_'+pos+ref+'T', chr+'_'+pos+ref+'G']
								hits = [x in HITPOS for x in snpl]
								if any(hits):
									index = hits.index(True)
									snp = snpl[index]
									row[3] = ref
									row[4] = snp[-1]
									for k in range(9,len(row)):
										#row[k] = str(row[k]).replace('|','/')
										if "/" in  str(row[k]) and str(row[k]) != '0/0':
											row[k] = '0/0'
											print('WARNING: polymorphic found in monomorphic SNP '+_snp+' (read) '+snp+' (written) <- Set to monomorphic')
										if "|" in str(row[k]) and str(row[k]) != '0|0':
											row[k] = '0|0'
											print('WARNING: polymorphic found in monomorphic SNP '+_snp+' (read) '+snp+' (written) <- Set to monomorphic')
									fcsv.writerow(row)
									#print('WARNING: MONOMORPH '+_snp+' '+snp)
									i += 1
									printi = True
								
								else:
									snpl = [chr+'_'+pos+'A'+ref, chr+'_'+pos+'C'+ref, chr+'_'+pos+'T'+ref, chr+'_'+pos+'G'+ref]
									hits = [x in HITPOS for x in snpl]
									if any(hits):
										index = hits.index(True)
										snp = snpl[index]
										row[3] = snp[-2]
										row[4] = ref
										for k in range(9,len(row)):
											#row[k] = str(row[k]).replace('|','/')
											if "/" in  str(row[k]) and str(row[k]) != '0/0':
												print('WARNING: polymorphic found in monomorphic SNP '+_snp+' (read) '+snp+' (written) <- Set to monomorphic')
												row[k] = '1/1'
											if "|" in str(row[k]) and str(row[k]) != '0|0':
												print('WARNING: polymorphic found in monomorphic SNP '+_snp+' (read) '+snp+' (written) <- Set to monomorphic')
												row[k] = '1|1'

										fcsv.writerow(row)
										#print('WARNING: MONOMORPH & ALT<->REF '+_snp+' '+snp)
										i += 1
										printi = True			
																					
						if i % 1000 == 0 and i > 0 and printi:
							print('Writing ' + str(i) + ' records')
					else:
						fcsv.writerow(row)
						
				print('Writing a total of ' + str(i) + ' records')
			fo.close()
			f.close()	
	
	#print(HIT)

if __name__ == "__main__":
	main(sys.argv[1:])
