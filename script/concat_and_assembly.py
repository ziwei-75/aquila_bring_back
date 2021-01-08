import subprocess
import portion as p
import numpy as np
import os
import argparse
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool 

parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', type=str, default="less_covered_regions", help='output directory for assembly and merging results for less covered regions')
parser.add_argument('--num_threads', type=int, default="6", help='number of threads for the process')

args = parser.parse_args()
output_dir = args.output_dir
num_threads = args.num_threads

def loadFasta(filename):
    """ Parses a classically formatted and possibly 
        compressed FASTA file into two lists. One of 
        headers and a second list of sequences.
        The ith index of each list correspond."""
    fp = open(filename, 'r')
    # split at headers
    data = fp.read().split('>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)     
    headers = []
    sequences = []
    for sequence in data:
        lines = sequence.split('\n')
        headers.append(lines.pop(0))
        # add an extra "+" to make string "1-referenced"
        sequences.append( ''.join(lines))
    return (headers, sequences)

def concat_reads(chromo,hp1_less_covered_regions,hp2_less_covered_regions):
    for region in hp1_less_covered_regions:
        start,end = region.split(",")
        start,end = int(start),int(end)
        read_dir = "./%s/chromo%s_phased_reads/"%(output_dir,chromo)
        hp1_reads = read_dir + "%s_%s_reads/%s_%s_hp1.fastq"%(start,end,start,end)
        discard_reads = read_dir + "%s_%s_reads/%s_%s_not_found.fastq"%(start,end,start,end)
        concat_reads = open(hp1_reads).read() + open(discard_reads).read()

        concat_dir = "./%s/chromo%s_concat_reads/%s_%s_reads/"%(output_dir,chromo,start,end)
        process = subprocess.Popen("mkdir %s"%concat_dir,shell=True)
        process.wait()
        concat_file = concat_dir + "%s_%s_hp1.fastq"%(start,end)
        open(concat_file,"w").write(concat_reads)
        
    for region in hp2_less_covered_regions:
        start,end = region.split(",")
        start,end = int(start),int(end)
        read_dir = "./%s/chromo%s_phased_reads/"%(output_dir,chromo)
        hp2_reads = read_dir + "%s_%s_reads/%s_%s_hp2.fastq"%(start,end,start,end)
        discard_reads = read_dir + "%s_%s_reads/%s_%s_not_found.fastq"%(start,end,start,end)
        concat_reads = open(hp2_reads).read() + open(discard_reads).read()

        concat_dir = "./%s/chromo%s_concat_reads/%s_%s_reads/"%(output_dir,chromo,start,end)
        if os.path.exists(concat_dir) == False:
            process = subprocess.Popen("mkdir %s"%concat_dir,shell=True)
            process.wait()
        concat_file = concat_dir + "%s_%s_hp2.fastq"%(start,end)
        open(concat_file,"w").write(concat_reads)

# assembly
def perform_assembly(chromo,chromo_concat_dir):
    readDir = os.listdir(chromo_concat_dir)
    hp1_count,hp2_count = 0,0
    for r_dir in readDir:
        r_dir = chromo_concat_dir + r_dir
        fastqFile = os.listdir(r_dir)
        for f_file in fastqFile:
            start,end,hp = f_file.split("_")
            f_file = r_dir + "/" + f_file
            hp = hp[:3]
            if hp == "hp1":
                hp1_count +=1 
            elif hp == "hp2":
                hp2_count += 1
            process = subprocess.Popen("mkdir ./%s/chromo%s_assembly/%s_%s_assembly/"%(output_dir,chromo,start,end),shell=True)
            process.wait()
            assemblyPath = "./%s/chromo%s_assembly/%s_%s_assembly/%s"%(output_dir,chromo,start,end,hp)
            contig_file = assemblyPath + "/contigs.fasta"
            spades_cmd = "python ./script/SPAdes-3.13.0-Linux/bin/spades.py -m 50 --phred-offset 33 --12 %s -o %s" %(f_file,assemblyPath)

            if os.path.exists(contig_file) == False:
                process = subprocess.Popen(spades_cmd,shell=True)
                process.wait()
        
chromoList = range(1,24)
def merge_reads_generate_bams(offset,stride):
    for i in range(offset, len(chromoList), stride):
        chromo = chromoList[i]
        hp1_less_covered_regions = open("./%s/chromo%s_hp1_merged_less_covered_regions.txt"%(output_dir,chromo)).read().split("\n")
        hp2_less_covered_regions = open("./%s/chromo%s_hp2_merged_less_covered_regions.txt"%(output_dir,chromo)).read().split("\n")
        chromo_concat_dir = "./%s/chromo%s_concat_reads/"%(output_dir,chromo)
        process = subprocess.Popen("mkdir %s"%chromo_concat_dir,shell=True)
        process.wait()
        concat_reads(chromo,hp1_less_covered_regions,hp2_less_covered_regions)
        assembly_path = "./%s/chromo%s_assembly/"%(output_dir,chromo)
        process = subprocess.Popen("mkdir %s"%assembly_path,shell=True)
        process.wait()
        perform_assembly(chromo,chromo_concat_dir)
    
def helper(i):
    return merge_reads_generate_bams(i, num_threads)

pool = ThreadPool(num_threads) 
results = pool.map(helper, range(num_threads))


    