import os
import subprocess
import pysam
import portion as p
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool 

parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', type=str, default="less_covered_regions", help='output directory for assembly and merging results for less covered regions')
parser.add_argument('--correct_bam_output_dir', type=str, help='output directory for correct bam file')
parser.add_argument('--num_threads', type=int, default="6", help='number of threads for the process')

args = parser.parse_args()
correct_bam_output_dir = args.correct_bam_output_dir
output_dir = args.output_dir
num_threads = args.num_threads

def map_contigs(chromo):
    assembly_dirs = os.listdir("./%s/chromo%s_assembly/"%(output_dir,chromo))
    failed_assemblies = []
    for a_dir in assembly_dirs:
        complete_a_dir = "./%s/chromo%s_assembly/"%(output_dir,chromo) + a_dir
        for hp in os.listdir(complete_a_dir):
            contig_file = "./%s/chromo%s_assembly/%s/%s/contigs.fasta"%(output_dir,chromo,a_dir,hp)
            if os.path.exists(contig_file) == False:
                failed_assemblies += [contig_file]
            output_sam = "./%s/chromo%s_assembly/%s/%s/output.sam"%(output_dir,chromo,a_dir,hp)
            if chromo == 23:
                indexFile = "./bam_index/chrX.mmi"
            else:
                indexFile = "./bam_index/chr%s.mmi"%chromo
#             print(indexFile)
            map_cmd = "minimap2 -a %s %s > %s " %(indexFile,contig_file,output_sam)
            process = subprocess.Popen(map_cmd,shell=True)
            process.wait()
            sorted_output_sam = "./%s/chromo%s_assembly/%s/%s/sorted_output.sam"%(output_dir,chromo,a_dir,hp)
            sort_cmd = "samtools sort %s -o %s"%(output_sam,sorted_output_sam)
            process = subprocess.Popen(sort_cmd,shell=True)
            process.wait()

def cal_length(interval):
    return interval.upper - interval.lower

def calculate_covered_bases(reads,region_interval):
    intervals = []
    
    count = 0
    for r in reads:
        if r.is_secondary:
            continue
        read_start = r.reference_start
        read_end = r.reference_end
        # unmapped reads
        if read_start == -1:
            continue
        else:
            length = read_end - read_start
            if length < 500:
                continue
            count += 1
            if p.closed(read_start,read_end).overlaps(region_interval):
                interval = p.closed(read_start,read_end).intersection(region_interval)
                intervals += [interval]
            else:
                continue
    
    if len(intervals) > 0:
        merged_interval = intervals[0]
        for i in range(1,len(intervals)):
            interval = intervals[i]
            merged_interval = merged_interval.union(interval)
        try: 
            covered_length = sum([cal_length(i) for i in merged_interval])
        except:
            print(interval.lower==+inf,interval.upper==-inf)
    else:
        covered_length = 0
    return count, covered_length

def cal_concat_bases(chromo,hp1_less_covered_regions,hp2_less_covered_regions):
    hp1_concat_bases_convered = {}
    hp1_concat_number_reads = {}
    for start,end in hp1_less_covered_regions:
    #     start,end = region.split(",")
        start,end = int(start),int(end)
        b_file = "./%s/chromo%s_assembly/%s_%s_assembly/%s/sorted_output.sam"%(output_dir,chromo,start,end,"hp1")
        bamFile = pysam.AlignmentFile(b_file)
        reads = bamFile.fetch()
        region_interval = p.closed(start,end)
        count, bases = calculate_covered_bases(reads,region_interval)
        hp1_concat_bases_convered[(start,end)] = bases/(end-start)
        hp1_concat_number_reads[(start,end)] = count

    hp2_concat_bases_convered = {}
    hp2_concat_number_reads = {}
    for start,end in hp2_less_covered_regions:
    #     start,end = region.split(",")
        start,end = int(start),int(end)
        bamFile = "./%s/chromo%s_assembly/%s_%s_assembly/%s/sorted_output.sam"%(output_dir,chromo,start,end,"hp2")
        bamFile = pysam.AlignmentFile(bamFile)
        reads = bamFile.fetch()
        region_interval = p.closed(start,end)
        count, bases = calculate_covered_bases(reads,region_interval)
        hp2_concat_bases_convered[(start,end)] = bases/(end-start)
        hp2_concat_number_reads[(start,end)] = count
    return hp1_concat_bases_convered,hp2_concat_bases_convered,hp1_concat_number_reads,hp2_concat_number_reads
    
    
def cal_aquila_reads(chromo,hp1_less_covered_regions,hp2_less_covered_regions,aquila_bamFile1,aquila_bamFile2):
    hp1_bases_covered = {}
    hp1_number_reads = {}
    for start,end in hp1_less_covered_regions:
        start,end = int(start),int(end)
        if str(chromo) == '23':
            reads = aquila_bamFile2.fetch('chrX',start,end)
        else:
            reads = aquila_bamFile2.fetch('chr' + str(chromo),start,end)
        region_interval = p.closed(start,end)
        count, bases = calculate_covered_bases(reads,region_interval)
        hp1_bases_covered[(start,end)] = bases/(end-start)
        hp1_number_reads[(start,end)] = count

    hp2_bases_covered = {}
    hp2_number_reads = {}
    for start,end in hp2_less_covered_regions:
        start,end = int(start),int(end)
        if str(chromo) == '23':
            reads = aquila_bamFile2.fetch('chrX',start,end)
        else:
            reads = aquila_bamFile2.fetch('chr' + str(chromo),start,end)
        region_interval = p.closed(start,end)
        count, bases = calculate_covered_bases(reads,region_interval)
        hp2_bases_covered[(start,end)] = bases/(end-start)
        hp2_number_reads[(start,end)] = count
    return hp1_bases_covered,hp2_bases_covered,hp1_number_reads,hp2_number_reads

def compute_improved_assemblies(chromo):
    hp1_less_covered_regions = open("./%s/chromo%s_hp1_merged_less_covered_regions.txt"%(output_dir,chromo)).read().split("\n")
    hp2_less_covered_regions = open("./%s/chromo%s_hp2_merged_less_covered_regions.txt"%(output_dir,chromo)).read().split("\n")
    hp1_less_covered_regions = [region.split(",") for region in hp1_less_covered_regions]
    hp2_less_covered_regions = [region.split(",") for region in hp2_less_covered_regions]

    hp1_concat_bases_convered,hp2_concat_bases_convered,hp1_concat_number_reads,hp2_concat_number_reads = cal_concat_bases(chromo,hp1_less_covered_regions,hp2_less_covered_regions)

    aquila_bamFile1 = pysam.AlignmentFile("./%s/Aquila_Contig_chr%s_hp1_sorted.bam"%(correct_bam_output_dir,chromo), "rb")
    aquila_bamFile2 = pysam.AlignmentFile("./%s/Aquila_Contig_chr%s_hp2_sorted.bam"%(correct_bam_output_dir,chromo), "rb")

    hp1_bases_covered,hp2_bases_covered,hp1_number_reads,hp2_number_reads = cal_aquila_reads(chromo,hp1_less_covered_regions,hp2_less_covered_regions,aquila_bamFile1,aquila_bamFile2)

    # hp1 assembly
    hp1_percentage, hp1_concat_percentage = [],[]
    hp1_reads, hp1_concat_reads = [],[]
    for start,end in hp1_less_covered_regions:
        start,end = int(start),int(end)
        hp1_percentage += [hp1_bases_covered[(start,end)]]
        hp1_concat_percentage += [hp1_concat_bases_convered[(start,end)]]
        hp1_reads += [hp1_number_reads[(start,end)]]
        hp1_concat_reads += [hp1_concat_number_reads[(start,end)]]

    # hp2 assembly
    hp2_percentage, hp2_concat_percentage = [],[]
    hp2_reads, hp2_concat_reads = [],[]
    for start,end in hp2_less_covered_regions:
        start,end = int(start),int(end)
        hp2_percentage += [hp2_bases_covered[(start,end)]]
        hp2_concat_percentage += [hp2_concat_bases_convered[(start,end)]]
        hp2_reads += [hp2_number_reads[(start,end)]]
        hp2_concat_reads += [hp2_concat_number_reads[(start,end)]]


    hp1_summaries = []
    for i in range(len(hp1_less_covered_regions)):
        start,end = hp1_less_covered_regions[i]
        start,end = int(start),int(end)
        percentage = hp1_bases_covered[(start,end)]
        concat_percentage = hp1_concat_bases_convered[(start,end)]
        reads = hp1_number_reads[(start,end)]
        concat_reads = hp1_concat_number_reads[(start,end)]
        hp1_summaries += [[start,end,percentage,concat_percentage,reads,concat_reads]]

    hp1_improved_assemblies = []
    for summary in hp1_summaries:
        start,end,percentage,concat_percentage,read,concat_read = summary
        # a lot of improvement in coverage 
        if concat_percentage > percentage + 0.5:
            hp1_improved_assemblies += [summary]
        elif concat_percentage - percentage > 0.2:
            if concat_read - read < 5:
                hp1_improved_assemblies += [summary]
        elif concat_percentage >= percentage:
            if concat_read < read:
                hp1_improved_assemblies += [summary]
                
    hp2_summaries = []
    for i in range(len(hp2_less_covered_regions)):
        start,end = hp2_less_covered_regions[i]
        start,end = int(start),int(end)
        percentage = hp2_bases_covered[(start,end)]
        concat_percentage = hp2_concat_bases_convered[(start,end)]
        reads = hp2_number_reads[(start,end)]
        concat_reads = hp2_concat_number_reads[(start,end)]
        hp2_summaries += [[start,end,percentage,concat_percentage,reads,concat_reads]]

    hp2_improved_assemblies = []
    for summary in hp2_summaries:
        start,end,percentage,concat_percentage,read,concat_read = summary
        # a lot of improvement in coverage 
        if concat_percentage > percentage + 0.5:
            hp2_improved_assemblies += [summary]
        elif concat_percentage - percentage > 0.2:
            if concat_read - read < 5:
                hp2_improved_assemblies += [summary]
        elif concat_percentage >= percentage:
            if concat_read < read:
                hp2_improved_assemblies += [summary]
    return hp1_less_covered_regions,hp2_less_covered_regions,hp1_improved_assemblies,hp2_improved_assemblies
   
chromoList = range(1,24)
def map_and_decide(offset,stride):
    for i in range(offset, len(chromoList), stride):
        chromo = chromoList[i]
        map_contigs(chromo)
        hp1_less_covered_regions,hp2_less_covered_regions,hp1_improved_assemblies,hp2_improved_assemblies = compute_improved_assemblies(chromo)

        hp1_improved_file = "./%s/improved_assemblies/chr%s_hp1_improved_assemblies.txt"%(output_dir,chromo)
        with open(hp1_improved_file,"w") as f:
            for start,end,_,_,_,_ in hp1_improved_assemblies:
                f.write(",".join([str(start),str(end)]) + "\n")

        hp2_improved_file = "./%s/improved_assemblies/chr%s_hp2_improved_assemblies.txt"%(output_dir,chromo)
        with open(hp2_improved_file,"w") as f:
            for start,end,_,_,_,_ in hp2_improved_assemblies:
                f.write(",".join([str(start),str(end)]) + "\n")
            
process = subprocess.Popen("mkdir ./%s/improved_assemblies/"%output_dir,shell=True)
process.wait()

def helper(i):
    return map_and_decide(i, num_threads)

pool = ThreadPool(num_threads) 
results = pool.map(helper, range(num_threads))
