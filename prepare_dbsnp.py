#!/usr/bin/python
"""
This scripts prepares the dbsnpfiles used by the imprinting pipeline.
It requires files to be downloaded from ftp://ftp.ensembl.org/pub/release-94/variation/vcf/homo_sapiens/homo_sapiens-chr*.vcf.gz
"""

"""
Imports
"""
import sys
import os
import traceback
import re
import itertools
import glob
from multiprocessing import Pool
import argparse

"""
Parameters
"""

parser = argparse.ArgumentParser(description='Create dbSNP files compatible with imprinting pipeline.')
parser.add_argument('--ensembl', dest="ensembl", help="ensembl variation vcf files.")
args = parser.parse_args()
"""
Locations
"""

dbsnp_in = "../data/"
dbsnp_out = "../dbsnp/"


"""
Placeholder variables
"""
vcf_files = []

"""
Functions
"""
def parseindels(original, new, info):
    if info == "deletion":
        newsnp = original[1:] + "/-"
    elif info == "insertion":
        newsplit= new.split(",")
        newnew = ""
        for split in newsplit:
            newnew = newnew + "/" + split[1:]
        newsnp = "-" + newnew
    else:
        newsnp = original + "/" + new.replace(',', '/')
    return newsnp

def parseEvidence(info):
    '''
    :param info: info string of the vcf
    :return: evidence fields as a string

    https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html
    '''
    evidences = ["E_Cited", "E_Multiple_observations", "E_Freq", "E_TOPMed", "E_Hapmap", "E_Phenotype_or_Disease", "E_ESP", "E_gnomAD", "E_1000G", "E_ExAC"]
    elements = info.split(';')
    evidencestring = ""
    snptype = ""
    for el in elements:
        if el in evidences:
            evidencestring = evidencestring + el + ";"
        elif "TSA" in el:
            snptype = el.split("=")[1]

    return evidencestring, snptype

def parseVCF(vcf,chr_annotations):
    with open(vcf, "r") as vcf_file:
        for line in vcf_file.readlines():
            if line.__contains__('#'):
                continue
            else:
                array = line.split("\t")
                chromosome = array[0]
                filename = "chr_%s.txt" % chromosome
                outfile = open(filename, "a+")
                annotations = chr_anotations.get(chromosome)
                geneannotation = annotations.get(array[2], "")
                info = array[7].split(';')[1].split("=")[1]
                snp = parseindels(array[3], array[4], info)
                evidence, snptype = parseEvidence(array[7].rstrip())
                newline = array[0] + "\t" + array[1] + "\t" + snp + "\t" + snptype + "\t" + array[2] + "\t" + geneannotation + "\t" + evidence + "\n"
                outfile.write(newline)

"""
Get file list
"""
"""
for filename in glob.iglob(args.ensembl + "/*.vcf"):
    vcf_files.append(filename)
for file in vcf_files:
    parseVCF(file)
"""
chromosomes = ["1", "2", "3", "4", "5", "6","7","8","9","10","11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
chr_anotations = {}
for chromosome in chromosomes:
    dbsnpfile = "./dbsnp/chr_%s.txt" % chromosome
    dict = {}
    with open(dbsnpfile, "r") as genefile:
        for geneline in genefile.readlines():
            genearray = geneline.split("\t")
            if genearray.__len__() > 12:
                dict['rs' + genearray[0]] = genearray[12]
        chr_anotations[chromosome] = dict

parseVCF("/data/dbsnp/00-common_all.vcf", chr_anotations)



