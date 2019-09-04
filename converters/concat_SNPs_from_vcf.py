#!/usr/bin/env python

"""
Name: concat_SNPs_from_vcf.py
Author: Michael G. Harvey
Date: 21 July 2015
Description: This script concatenated the SNPs in a vcf file into an alignment.
Usage: python adegenet_from_vcf.py in_file out_file

"""

import os
import sys
import argparse
import random


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_file",
			type=str,
			help="""The input vcf file from freebayes"""
		)
	parser.add_argument(
			"out_file",
			type=str,
			help="""The output fasta file"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	infile = open("{0}".format(args.in_file), 'r')
	outfile = open("{0}".format(args.out_file), 'wb')
	lines = infile.readlines()
	columns_a = list()
	columns_b = list()
	for line in lines:
		if line.startswith("#CHROM"):
			parts = line.split()
			names = parts[9:len(parts)]
		elif not line.startswith("#"):
			alleles_a = list()
			alleles_b = list()
			parts = line.split()
			alleles = list(parts[3])
			alt_alleles = parts[4].split(",")
			for alt_allele in alt_alleles:
				alleles.append(alt_allele)
			print alleles
			for i, name in enumerate(names):
				ind = parts[i+9].split(":")[0]
				allele_a = ind.split("/")[0]
				allele_b = ind.split("/")[1]
				if allele_a == ".":
					alleles_a.append("N")
					alleles_b.append("N")
				else:
					print allele_a
					# Randomly assign alleles to either genome copy
					rand = random.randint(0,1)
					if rand == 0:
						alleles_a.append(alleles[int(allele_a)])
						alleles_b.append(alleles[int(allele_b)])
					elif rand == 1:
						alleles_a.append(alleles[int(allele_b)])
						alleles_b.append(alleles[int(allele_a)])
			if not "N" in alleles_a:
				if not "N" in alleles_b:			
					columns_a.append(alleles_a)
					columns_b.append(alleles_b)
			print alleles_a
			print alleles_b
			
	outfile.write("{0}	{1}\n".format((len(columns_a[1])*2), len(columns_a)))
	for i in range(len(columns_a[1])):
		outfile.write("{0}a	".format(names[i]))
		for column_a in columns_a:
			outfile.write(column_a[i])
		outfile.write("\n")		
		outfile.write("{0}b	".format(names[i]))
		for column_b in columns_b:
			outfile.write(column_b[i])
		outfile.write("\n")
	infile.close()
	outfile.close()
									
if __name__ == '__main__':
    main()	
	