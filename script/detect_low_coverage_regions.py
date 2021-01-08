import os
import portion as p
import pysam
import argparse
import subprocess
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool 

parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', type=str, help='input directory that contains aquila assembly')
parser.add_argument('--correct_bam_output_dir', type=str, help='output directory for correct bam file')
parser.add_argument('--output_dir', type=str, default="less_covered_regions", help='output directory for assembly and merging results for less covered regions')
parser.add_argument('--num_threads', type=int, default="6", help='number of threads for the process')

args = parser.parse_args()
correct_bam_output_dir = args.correct_bam_output_dir
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

def cal_length(interval):
    return interval.upper - interval.lower

def calculate_covered_bases(reads,region_interval):
    intervals = []
    
    count = 0
    for r in reads:
        count += 1
        read_start = r.reference_start
        read_end = r.reference_end
        # unmapped reads
        if read_start == -1:
            continue
        else:
            interval = p.closed(read_start,read_end).intersection(region_interval)

        intervals += [interval]
    
    if len(intervals) > 0:
        merged_interval = intervals[0]
        for i in range(1,len(intervals)):
            interval = intervals[i]
            merged_interval = merged_interval.union(interval)
        covered_length = sum([cal_length(i) for i in merged_interval])
    else:
        covered_length = 0
    return count, covered_length

def cal_covered_bases(chromo,bamFile1,bamFile2,bothHap):
    chromo = str(chromo)
    num_covered_bases1 = {}
    for block in bothHap:
        s,e = block
        if e - s < 10000:
            if chromo == '23':
                reads = bamFile1.fetch('chrX',s,e)
            else:
                reads = bamFile1.fetch('chr' + chromo,s,e)   
            region_interval = p.closed(s,e)
            count, bases = calculate_covered_bases(reads,region_interval)
            num_covered_bases1[(s,e)] = bases
        else:
            for i in range(s,e,10000):
                start = i 
                if start + 10000 < e:
                    end = start + 10000
                else:
                    end = e
                if chromo == '23':
                    reads = bamFile1.fetch('chrX',start,end)
                else:
                    reads = bamFile1.fetch('chr' + chromo,start,end)
                region_interval = p.closed(start,end)
                count, bases = calculate_covered_bases(reads,region_interval)
                num_covered_bases1[(start,end)] = bases
    num_covered_bases2 = {}
    for block in bothHap:
        s,e = block
        if e - s < 10000:
            if chromo == '23':
                reads = bamFile2.fetch('chrX',s,e)
            else:
                reads = bamFile2.fetch('chr' + chromo,s,e)   
            region_interval = p.closed(s,e)
            count, bases = calculate_covered_bases(reads,region_interval)
            num_covered_bases2[(s,e)] = bases
        else:
            for i in range(s,e,10000):
                start = i 
                if start + 10000 < e:
                    end = start + 10000
                else:
                    end = e
                if chromo == '23':
                    reads = bamFile2.fetch('chrX',start,end)
                else:
                    reads = bamFile2.fetch('chr' + chromo,start,end)
                region_interval = p.closed(start,end)
                count, bases = calculate_covered_bases(reads,region_interval)
                num_covered_bases2[(start,end)] = bases
    return num_covered_bases1,num_covered_bases2

def find_less_covered_regions(num_covered_bases1,num_covered_bases2):
    hp1_less_covered_regions = []
    hp2_less_covered_regions = []
    for region, bases1 in num_covered_bases1.items():
        bases2 = num_covered_bases2[region]
        try:
            percentage1 = bases1/(region[1] - region[0])
            percentage2 = bases2/(region[1] - region[0])
            if percentage1 < 0.8:
                hp1_less_covered_regions += [region]
            if percentage2 < 0.8: 
                hp2_less_covered_regions += [region]
        except:
            print(region)
    return hp1_less_covered_regions,hp2_less_covered_regions

def write_less_covered_regions(chromo,hp1_less_covered_regions,hp2_less_covered_regions):
    write_hp1_less_covered_regions = []
    for start,end in hp1_less_covered_regions:
        write_hp1_less_covered_regions += [",".join([str(start),str(end)])]
    write_hp1_less_covered_regions = "\n".join(write_hp1_less_covered_regions)

    write_hp2_less_covered_regions = []
    for start,end in hp2_less_covered_regions:
        write_hp2_less_covered_regions += [",".join([str(start),str(end)])]
    write_hp2_less_covered_regions = "\n".join(write_hp2_less_covered_regions)
    
    total_less_covered_regions = set(hp1_less_covered_regions).union(set(hp2_less_covered_regions))
    write_total_less_covered_regions = []
    for start,end in total_less_covered_regions:
        write_total_less_covered_regions += [",".join([str(start),str(end)])]
    write_total_less_covered_regions = "\n".join(write_total_less_covered_regions)
    
    open("%s/chromo%s_hp1_less_covered_regions.txt"%(output_dir,chromo),"w").write(write_hp1_less_covered_regions)
    open("%s/chromo%s_hp2_less_covered_regions.txt"%(output_dir,chromo),"w").write(write_hp2_less_covered_regions)
    open("%s/chromo%s_total_less_covered_regions.txt"%(output_dir,chromo),"w").write(write_total_less_covered_regions)

    
chromoList = range(1,24)
def main_less_covered_regions(offset,stride):
    for i in range(offset, len(chromoList), stride):
        chromo = chromoList[i]
        print("find less covered region for chromo %s"%chromo)
        read_fasta_dir = "./%s/Local_Assembly_by_chunks/chr%s_files_cutPBHC"%(input_dir,chromo)
        hp1_phase_blocks,hp2_phase_blocks,bothHap,totalHap = parse_phase_blocks(read_fasta_dir)
        aquila_bamFile1 = pysam.AlignmentFile("./%s/Aquila_Contig_chr%s_hp1_sorted.bam"%(correct_bam_output_dir,chromo), "rb")
        aquila_bamFile2 = pysam.AlignmentFile("./%s/Aquila_Contig_chr%s_hp2_sorted.bam"%(correct_bam_output_dir,chromo), "rb")
        num_covered_bases1,num_covered_bases2 = cal_covered_bases(chromo,aquila_bamFile1,aquila_bamFile2,bothHap)
        hp1_less_covered_regions,hp2_less_covered_regions = find_less_covered_regions(num_covered_bases1,num_covered_bases2)
        write_less_covered_regions(chromo,hp1_less_covered_regions,hp2_less_covered_regions)
  
process = subprocess.Popen("mkdir %s"%output_dir,shell=True)
process.wait()

def helper(i):
    return main_less_covered_regions(i, num_threads)

pool = ThreadPool(num_threads) 
results = pool.map(helper, range(num_threads))

