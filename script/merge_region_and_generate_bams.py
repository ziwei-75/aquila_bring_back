import pysam
import os
import subprocess
import portion as p
import argparse
from collections import defaultdict
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool 


parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', type=str, help='input directory that contains aquila assembly')
parser.add_argument('--output_dir', type=str, default="less_covered_regions", help='output directory for assembly and merging results for less covered regions')
parser.add_argument('--num_threads', type=int, default="6", help='number of threads for the process')

args = parser.parse_args()
output_dir = args.output_dir
input_dir = args.input_dir
num_threads = args.num_threads

def parse_phase_blocks(read_fasta_dir):
    phase_blocks = os.listdir(read_fasta_dir)
    hp1_phase_blocks = []
    hp2_phase_blocks = []
    for p_f in phase_blocks:
        start,end = p_f.split("_")[2:4]
        if p_f[-9:-6] == "hp1":
            hp1_phase_blocks += [(int(start),int(end))]
        else:
            assert p_f[-9:-6] == "hp2"
            hp2_phase_blocks += [(int(start),int(end))]
        
    # filter overlap_regions
    overlapped_regions = set()
    for i in range(len(hp1_phase_blocks)):
        region_i = p.open(hp1_phase_blocks[i][0],hp1_phase_blocks[i][1])
        for j in range(i+1, len(hp1_phase_blocks)):
            region_j = p.open(hp1_phase_blocks[j][0],hp1_phase_blocks[j][1])
            if region_i.overlaps(region_j):
    #             print(region_i,region_j)
                overlapped_regions.add(hp1_phase_blocks[i])
                overlapped_regions.add(hp1_phase_blocks[j])
            
        # filter overlap_regions
    for i in range(len(hp2_phase_blocks)):
        region_i = p.open(hp2_phase_blocks[i][0],hp2_phase_blocks[i][1])
        for j in range(i+1, len(hp2_phase_blocks)):
            region_j = p.open(hp2_phase_blocks[j][0],hp2_phase_blocks[j][1])
            if region_i.overlaps(region_j):
    #             print(region_i,region_j)
                overlapped_regions.add(hp2_phase_blocks[i])
                overlapped_regions.add(hp2_phase_blocks[j])
            
    # remove overlap regions
    for overlap in overlapped_regions:
        try:
            hp1_phase_blocks.remove(overlap)
        except:
            continue

    for overlap in overlapped_regions:
        try:
            hp2_phase_blocks.remove(overlap)
        except:
            continue
    
    bothHap = []
    for block1 in hp1_phase_blocks:
        for block2 in hp2_phase_blocks:
            if block1 == block2:
                bothHap += [block1]
    totalHap = set(hp1_phase_blocks) | set(hp2_phase_blocks)
    return hp1_phase_blocks,hp2_phase_blocks,bothHap,totalHap


def merge_intervals(hp1_less_covered_regions,hp2_less_covered_regions,bothHap):
    for i in range(len(hp1_less_covered_regions[::])):
        region = hp1_less_covered_regions[i]
        start,end = region.split(",")
        start,end = int(start),int(end)
        hp1_less_covered_regions[i] = [start,end]

    hp1_less_covered_regions_dict = {}
    for i in range(len(hp1_less_covered_regions)):
        current_start,current_end = hp1_less_covered_regions[i]
        less_cover_interval = p.closed(current_start,current_end)
        for start,end in bothHap:
            hap_interval = p.closed(start,end)
            if hap_interval.contains(less_cover_interval):
                hp1_less_covered_regions_dict[hap_interval] = hp1_less_covered_regions_dict.get(hap_interval,[]) + [less_cover_interval]

    hp1_merged_intervals = []
    for hap, less_interval in hp1_less_covered_regions_dict.items():
        less_interval.sort(key=lambda x: x[0])
        merged_intervals = []
        my_start,my_end = less_interval[0].lower, less_interval[0].upper

        for i in range(len(less_interval)-1):
            current_start,current_end = less_interval[i].lower,less_interval[i].upper
            next_start,next_end = less_interval[i+1].lower,less_interval[i+1].upper

            # merge 
            if current_end == next_start:
                my_end = next_end
            else:
                merged_intervals += [ [my_start,my_end]]
                my_start = next_start
                my_end = next_end
            if i == len(less_interval) - 2:
                merged_intervals += [[my_start,my_end]]
        if len(less_interval) == 1:
            merged_intervals = [[less_interval[0].lower, less_interval[0].upper]]
        hp1_merged_intervals += merged_intervals

    for i in range(len(hp2_less_covered_regions[::])):
        region = hp2_less_covered_regions[i]
        start,end = region.split(",")
        start,end = int(start),int(end)
        hp2_less_covered_regions[i] = [start,end]

    hp2_less_covered_regions_dict = {}
    for i in range(len(hp2_less_covered_regions)):
        current_start,current_end = hp2_less_covered_regions[i]
        less_cover_interval = p.closed(current_start,current_end)
        for start,end in bothHap:
            hap_interval = p.closed(start,end)
            if hap_interval.contains(less_cover_interval):
                hp2_less_covered_regions_dict[hap_interval] = hp2_less_covered_regions_dict.get(hap_interval,[]) + [less_cover_interval]

    hp2_merged_intervals = []
    for hap, less_interval in hp2_less_covered_regions_dict.items():
        less_interval.sort(key=lambda x: x[0])
        merged_intervals = []
        my_start,my_end = less_interval[0].lower, less_interval[0].upper

        for i in range(len(less_interval)-1):
            current_start,current_end = less_interval[i].lower,less_interval[i].upper
            next_start,next_end = less_interval[i+1].lower,less_interval[i+1].upper

            # merge 
            if current_end == next_start:
                my_end = next_end
            else:
                merged_intervals += [ [my_start,my_end]]
                my_start = next_start
                my_end = next_end
            if i == len(less_interval) - 2:
                merged_intervals += [[my_start,my_end]]
        if len(less_interval) == 1:
            merged_intervals = [[less_interval[0].lower, less_interval[0].upper]]
        hp2_merged_intervals += merged_intervals

    hp1_filtered_intervals = []
    for start,end in hp1_merged_intervals:
        if end - start >= 10000:
            start,end = str(start),str(end)
            hp1_filtered_intervals += [",".join([start,end])]

    hp2_filtered_intervals = []
    for start,end in hp2_merged_intervals:
        if end - start >= 10000:
            start,end = str(start),str(end)
            hp2_filtered_intervals += [",".join([start,end])]

    total_filtered_intervals = set(hp1_filtered_intervals).union(set(hp2_filtered_intervals))
    
    return hp1_filtered_intervals,hp2_filtered_intervals,total_filtered_intervals


def write_merged_intervals(chromo,hp1_filtered_intervals,hp2_filtered_intervals,total_filtered_intervals):
    hp1_filtered_intervals = "\n".join(hp1_filtered_intervals)
    hp2_filtered_intervals = "\n".join(hp2_filtered_intervals)
    total_filtered_intervals = "\n".join(total_filtered_intervals)
    open("./%s/chromo%s_hp1_merged_less_covered_regions.txt"%(output_dir,chromo),"w").write(hp1_filtered_intervals)
    open("./%s/chromo%s_hp2_merged_less_covered_regions.txt"%(output_dir,chromo),"w").write(hp2_filtered_intervals)
    open("./%s/chromo%s_total_merged_less_covered_regions.txt"%(output_dir,chromo),"w").write(total_filtered_intervals)
    
    
def generate_region_reads_bams(chromo):
    assembly_breaks = open("./%s/chromo%s_total_merged_less_covered_regions.txt"%(output_dir,chromo)).read().split("\n")
    bamFile = pysam.AlignmentFile("./%s/possorted_bam.bam"%input_dir, "rb")
    process = subprocess.Popen("mkdir ./%s/chromo%s_read_bam"%(output_dir,chromo),shell=True)
    process.wait()
    for region in assembly_breaks:
        region = region.split(",")
        start,end = region
        start,end = int(start), int(end)
        outputFile = "./%s/chromo%s_read_bam/%s_%s_reads.bam"%(output_dir,chromo,start,end)
        if chromo == 23:
            reads = bamFile.fetch('chrX', start, end)
        else:
            reads = bamFile.fetch('chr' + str(chromo), start, end)
        
        newFile = pysam.AlignmentFile(outputFile, "wb", template=bamFile)

        for read in reads:
            newFile.write(read)

        newFile.close()
    bamFile.close()
    
    
chromoList = range(1,24)
if os.path.exists("./%s/possorted_bam.bam.bai"%input_dir) == False:
    index_p = subprocess.Popen("samtools index %s/possorted_bam.bam"%input_dir, shell=True)
    index_p.wait()
    print("finish indexing")
def merge_reads_generate_bams(offset,stride):
    for i in range(offset, len(chromoList), stride):
        chromo = chromoList[i]
        read_fasta_dir = "./%s/Local_Assembly_by_chunks/chr%s_files_cutPBHC"%(input_dir,chromo)
        hp1_phase_blocks,hp2_phase_blocks,bothHap,totalHap = parse_phase_blocks(read_fasta_dir)
        hp1_less_covered_regions = open("./%s/chromo%s_hp1_less_covered_regions.txt"%(output_dir,chromo)).read().split("\n")
        hp2_less_covered_regions = open("./%s/chromo%s_hp2_less_covered_regions.txt"%(output_dir,chromo)).read().split("\n")
        hp1_filtered_intervals,hp2_filtered_intervals,total_filtered_intervals = merge_intervals(hp1_less_covered_regions,hp2_less_covered_regions,bothHap)
        write_merged_intervals(chromo,hp1_filtered_intervals,hp2_filtered_intervals,total_filtered_intervals)
        generate_region_reads_bams(chromo)
    

def helper(i):
    return merge_reads_generate_bams(i, num_threads)

pool = ThreadPool(num_threads) 
results = pool.map(helper, range(num_threads))
    