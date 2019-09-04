#!/usr/bin/env python

"""

Name: output_sample_fastas.py

Author: Michael G. Harvey
Date: 4 September 2019


"""

import os
import sys
import re
import random
import argparse
import copy


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"seq_file",
			type=str,
			help="""The Velvet contigs file for this individual"""
		)
	parser.add_argument(
			"map_file",
			type=str,
			help="""The map of contigs to target loci for this sample from SqCL"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory of fasta files for this individual"""
		)
	return parser.parse_args()

def main():
	args = get_args()
	sample_name = map_file.split[-1]
	sample_name = sample_name.split("_matches.csv")[0]
	contigsfile = open("{0}".format(args.seq_file))
	mapfile = open("{0}".format(args.map_file))
	contigs = contigsfile.readlines()
	maps = mapfile.readlines()
	for contig in contigs: # For each line in contig file
		if contig.startswith('>'): # If header			
			contig_name = line.split('>').rstrip()	
			if seq: # If there was a prior locus
				
				# Repeat the following at end
				for map in maps: # For line in mapfile
					if not map.startswith("contig"): # Bypass header
						parts = map.split(',')
						map_contig = parts[0]
						if map_contig == contig_name:
							map_locus = parts[1]
							outfile = os.path.join(args.out_dir, map_locus, ".fasta")
							# Check if output file exists
							if outfile:
								outfile = open(outfile), 'a')
								outfile.write(">{0}".format(sample_name))
								outfile.write(''.join(seq))
							else:
								outfile = open(outfile), 'wb')
								outfile.write(">{0}".format(sample_name))
								outfile.write(''.join(seq))
							outfile.close()
					
			seq = list()
		else:
			seq.append(line.rstrip())
			
	# Repeat sequence output stuff here	
	
	outfile.close()
	contigsfile.close()
	mapfile.close()

if __name__ == '__main__':
    main()