import pysam
import os
import subprocess
import portion as P
from collections import defaultdict
import argparse
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool 
import psutil

parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', type=str, help='input directory that contains aquila assembly')
parser.add_argument('--output_dir', type=str, default="less_covered_regions", help='output directory for assembly and merging results for less covered regions')
parser.add_argument('--num_threads', type=int, default="6", help='number of threads for the process')

args = parser.parse_args()
output_dir = args.output_dir
input_dir = args.input_dir
num_threads = args.num_threads

def load_phase_blocks(read_fasta_dir):
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
            
    hp1_fastq_intervals = []
    for start, end in hp1_phase_blocks:
        hp1_fastq_intervals += [P.closed(start,end)]

    hp2_fastq_intervals = []
    for start, end in hp2_phase_blocks:
        hp2_fastq_intervals += [P.closed(start,end)]
    return hp1_fastq_intervals,hp2_fastq_intervals

def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch():
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
        
def merge_fastq_bam(chromo, hp1_overlap_fastq,hp2_overlap_fastq,bam_interval):
    start, end = bam_interval.lower,bam_interval.upper
    phased_reads_path = "./%s/chromo%s_phased_reads/%s_%s_reads"%(output_dir,chromo,start,end)
    readName1 = []
    if len(hp1_overlap_fastq) != 0:
        hp1_overlap_fastq = list(hp1_overlap_fastq)[0]
        start1, end1 = hp1_overlap_fastq.lower, hp1_overlap_fastq.upper
        fastq_file1 = "./%s/Local_Assembly_by_chunks/chr%s_files_cutPBHC/fastq_by_%s_%s_hp1.fastq"%(input_dir,chromo,start1,end1)
        fastq_file1 = open(fastq_file1).read().split("\n")[:-1]
        readName1 = []
        for i in range(len(fastq_file1)//4):
            name = fastq_file1[4*i]
            readName1 += [name[1:]]
    
    readName2 = []
    if len(hp2_overlap_fastq) != 0:
        hp2_overlap_fastq = list(hp2_overlap_fastq)[0]
        start2, end2 = hp2_overlap_fastq.lower, hp2_overlap_fastq.upper
        fastq_file2 = "./%s/Local_Assembly_by_chunks/chr%s_files_cutPBHC/fastq_by_%s_%s_hp2.fastq"%(input_dir,chromo,start1,end1)
        fastq_file2 = open(fastq_file2).read().split("\n")[:-1]
        readName2 = []
        for i in range(len(fastq_file2)//4):
            name = fastq_file2[4*i]
            readName2 += [name[1:]]


    samfile = pysam.AlignmentFile("./%s/chromo%s_read_bam/%s_%s_reads.bam"%(output_dir,chromo,start,end), "rb")
    hp1_reads = []
    hp2_reads = []
    notFound = []
    for read in samfile.fetch():
        readName = read.qname
        if readName in readName1:
            assert not readName in readName2
            hp1_reads += [readName]
        elif readName in readName2:
            assert not readName in readName1
            hp2_reads += [readName]
        else:
            notFound += [readName]

    # write reads
    mkdir_cmd = "mkdir ./%s/chromo%s_phased_reads/%s_%s_reads"%(output_dir,chromo, start,end)
    process = subprocess.Popen(mkdir_cmd,shell=True)
    process.wait()
    
    readPairs = read_pair_generator(samfile)
    hp1_reads_fastq = ""
    hp2_reads_fastq = ""
    notFound_fastq = ""
    for r1, r2 in readPairs:
        if r1.qname in hp1_reads: 
            hp1_reads_fastq += '@' + str(r1.qname)  + "\n" + str(r1.query_sequence) + "\n" + "+" + "\n" + str(r1.qual) + "\n"
            hp1_reads_fastq += '@' + str(r2.qname)  + "\n" + str(r2.query_sequence) + "\n" + "+" + "\n" + str(r2.qual) + "\n"

        elif r1.qname in hp2_reads: 
            hp2_reads_fastq += '@' + str(r1.qname) + "\n" + str(r1.query_sequence) + "\n" + "+" + "\n" + str(r1.qual) + "\n"
            hp2_reads_fastq += '@' + str(r2.qname) + "\n" + str(r2.query_sequence) + "\n" + "+" + "\n" + str(r2.qual) + "\n"

        else: 
            try:
                assert r1.qname in notFound
            except:
                print(r1.qname)
            notFound_fastq += '@' + str(r1.qname) + "\n" + str(r1.query_sequence) + "\n" + "+" + "\n" + str(r1.qual) + "\n"
            notFound_fastq += '@' + str(r2.qname) + "\n" + str(r2.query_sequence) + "\n" + "+" + "\n" + str(r2.qual) + "\n"

    with open("./%s/chromo%s_phased_reads/%s_%s_reads/%s_%s_hp1.fastq"%(output_dir,chromo,start,end,start,end),"w") as f:
        f.write(hp1_reads_fastq)

    with open("./%s/chromo%s_phased_reads/%s_%s_reads/%s_%s_hp2.fastq"%(output_dir,chromo,start,end,start,end),"w") as f:
        f.write(hp2_reads_fastq)

    with open("./%s/chromo%s_phased_reads/%s_%s_reads/%s_%s_not_found.fastq"%(output_dir,chromo,start,end,start,end),"w") as f:
        f.write(notFound_fastq)

    return 

def perform_phasing(chromo,hp1_fastq_intervals,hp2_fastq_intervals,bamFiles):
    for i in range(len(bamFiles)):
        f = bamFiles[i]
        if f[-3:] == "bam":
            start, end = f.split("_")[:2]
            bam_interval = P.open(int(start),int(end))
            hp1_overlap_fastq = set()
            hp2_overlap_fastq = set()
            for f_i in hp1_fastq_intervals:
                if f_i.contains(bam_interval):
                    hp1_overlap_fastq.add(f_i)

            for f_i in hp2_fastq_intervals:
                if f_i.contains(bam_interval):
                    hp2_overlap_fastq.add(f_i)

            if len(hp1_overlap_fastq) == 1 and len(hp2_overlap_fastq) == 1:
                merge_fastq_bam(chromo,hp1_overlap_fastq,hp2_overlap_fastq,bam_interval)

            elif len(hp1_overlap_fastq) == 1 or len(hp2_overlap_fastq) == 1:
                merge_fastq_bam(chromo,hp1_overlap_fastq,hp2_overlap_fastq,bam_interval)
            else:
                print(f_i)
                print("not right ")
                break

                
chromoList = range(1,24)
def merge_reads_generate_bams(offset,stride):
    for i in range(offset, len(chromoList), stride):
        chromo = chromoList[i]
        read_fasta_dir = "./%s/Local_Assembly_by_chunks/chr%s_files_cutPBHC"%(input_dir,chromo)
        hp1_fastq_intervals,hp2_fastq_intervals = load_phase_blocks(read_fasta_dir)
        reads = os.listdir("%s/chromo%s_read_bam/"%(output_dir,chromo))
        for r in reads:
            index_cmd = "samtools index ./%s/chromo%s_read_bam/%s"%(output_dir,chromo,r)
            process = subprocess.Popen(index_cmd,shell=True)
            process.wait()

        process = subprocess.Popen("mkdir ./%s/chromo%s_phased_reads"%(output_dir,chromo),shell=True)
        process.wait()
        bamFiles = os.listdir("./%s/chromo%s_read_bam/"%(output_dir,chromo))
        perform_phasing(chromo,hp1_fastq_intervals,hp2_fastq_intervals,bamFiles)
    

def helper(i):
    return merge_reads_generate_bams(i, num_threads)

pool = ThreadPool(num_threads) 
results = pool.map(helper, range(num_threads))
    
    