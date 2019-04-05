import configparser
import sys
import string
import os
import glob
import subprocess
from multiprocessing import Pool

def LoadConfig(file):
    """
    returns a dictionary with keys of the form
    <section>.<option> and the corresponding values
    """
    config = configparser.ConfigParser()
    config.read(file)
    inputs = config.items('INPUTS')
    inputDict = {}
    for list in inputs:
        print(list)
        inputDict[list[0]] = list[1]
    outputs = config.items("OUTPUTS")
    outputDict = {}
    for list in outputs:
        print(list)
        outputDict[list[0]] = list[1]
    return inputDict, outputDict

def getfiles(inputs):
    files={}
    globpattern = inputs.get('bam_path').strip("\'") + "/*.bam"
    print(globpattern)
    for file in glob.glob(globpattern):
        Sample = file.split("/")[-1].split(".")[0]
        files[Sample] = [file,]
    return files

def mkdirs(outputs):
    for key, value in outputs.items():
        cmd = "mkdir -p %s" %value
        os.system(cmd)

def flagstatLauncher(files, inputs, outputs):
    threads = int(inputs.get("threads"))
    bampath = inputs.get("bam_path")
    pool = Pool(processes=1)
    [pool.apply_async(flagstat, args=(files, key, value[0])) for key, value in files.items()]
    pool.close()
    pool.join()
def flagstat(files, file, path):
    cmd = "samtools flagstat %s" %path
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    values = files.get(file)
    for line in p.stdout:
        lala = str(line)
        lala2 = lala[2:-3]
        array = lala2.split(" ")
        values.append(array[0])
    files[file] = values
    print(files[file])


def indexfiles(files):
    for key, value in files.items():
        #check for index file
        index = value[0] + '.bai'
        if not os.path.isfile(index):
            cmd = "samtools index %s " % value[0]
            print(cmd)
            os.system(cmd)
        else:
            print(key + " already indexed.")

def snpCallingLauncher(files, inputs, outputs):
    filestring = ""
    for key, value in files.items():
        filestring = filestring + " " + value[0]
    chromosomes = inputs.get("chromosomes").split(",")
    bcfpath = outputs.get('bcf')
    vcfpath = outputs.get('vcf')
    pool = Pool(processes=int(inputs.get("threads")))
    [pool.apply_async(snpCalling, args=(x, filestring, inputs.get("reference_genome"), bcfpath, vcfpath )) for x in chromosomes]
    pool.close()
    pool.join()

def snpCalling(ch, bam_files, ref_genome, bcfloc, vcfloc):
    print("Launching mpileup for chromosome " + ch)
    cmd = "bcftools mpileup -f %s -r chr%s %s | bcftools call -mv -Oz -o %s/var_ch%s.bcf" %(ref_genome, ch, bam_files,bcfloc,ch)
    #print(cmd)
    #os.system(cmd)
    cmd2 = "bcftools view %s/var_ch%s.bcf | vcfutils.pl varFilter -d10 > %s/var_ch%s.vcf" % (bcfloc,ch, vcfloc, ch) #-d10 parameter
    #os.system(cmd2)

def bam2samLauncher(files, inputs, outputs):
    #for each bam launch bam2sam
    chromosomes = inputs.get("chromosomes").split(",")
    threads = int(inputs.get("threads"))
    sampath = outputs.get("sam")
    for key, value in files.items():
        pool = Pool(processes=threads)
        [pool.apply_async(bam2sam, args=(x, value[0], key, sampath)) for x in chromosomes]
        pool.close()
        pool.join()

def bam2sam(ch, bam,key, sampath):
    file = "%s/%s_chr%s_U.sam" %(sampath,key, ch)
    cmd = "samtools view -h %s chr%s | grep -w 'NH:i:1' > %s" %(bam,ch,file)
    print(cmd)
    #os.system(cmd)

def snpExtractLauncer(inputs, outputs):
    chromosomes = inputs.get("chromosomes").split(",")
    threads = int(inputs.get("threads"))
    vcf_path = outputs.get("vcf")
    dbsnp_path = inputs.get("dbsnp_path")
    variant_path = outputs.get("variant")
    pool = Pool(processes=threads)
    [pool.apply_async(snpExtract, args=(x, vcf_path, dbsnp_path, variant_path)) for x in chromosomes]
    pool.close()
    pool.join()


def snpExtract(chr, vcf_path, dbsnp_path, variant_path):
    vcf_file = "%s/var_ch%s.vcf" %(vcf_path, chr)
    dbsn_file = "%s/chr_%s.txt" %(dbsnp_path, chr)
    variant_file = "%s/variants_chr%s.txt" %(variant_path, chr)
    txt = "got %s %s %s %s " % (vcf_file, dbsn_file, variant_file, chr)
    print(txt)
    outputfile = open(variant_file, 'w')
    dbSNP = open(dbsn_file, "r")
    SNP = {}
    for line in dbSNP:
        line = line.rstrip("\r\n")
        a = line.split("\t")
        SNP[a[1]] = [a[2], a[3], a[4], a[5]]
    dbSNP.close()
    inputfile = open(vcf_file, "r")
    for line in inputfile:
        line = line.rstrip("\r\n")
        a = line.split("\t")
        if "chr" in a[0]:
            if "chr" not in chr:
                chr = "chr" + chr

        print(a[0] + "vs" + chr)
        if a[0] == chr:
            print(a)
            if a[1] in SNP:
                print(a[1] + "in SNP")
                info = a[7].split(";")
                if len(a[4]) == 1 and info[0] != "INDEL":
                    outputfile.write(
                        chr + "\t" + a[1] + "\t" + SNP[a[1]][0] + "\t" + SNP[a[1]][1] + "\t" + SNP[a[1]][
                            2] + "\t" + SNP[a[1]][3] + "\n")
    inputfile.close()
    outputfile.close()


def getSeqLauncher(inputs, outputs):
    chromosomes = inputs.get("chromosomes").split(",")
    threads = int(inputs.get("threads"))
    variant_path = outputs.get("variant")
    output_path = outputs.get("sequences")
    sam_path = outputs.get("sam")


    return



def main():
    config = sys.argv[1]
    inputs, outputs = LoadConfig(config)
    files = getfiles(inputs)
    mkdirs(outputs)
    flagstatLauncher(files, inputs, outputs)
    exit()
    indexfiles(files)
    snpCallingLauncher(files, inputs, outputs)
    #bam2samLauncher(files,inputs, outputs)
    snpExtractLauncer(inputs,outputs)
    #getSeqLauncher("1", files, outputs)


if __name__ == '__main__':
    main()
