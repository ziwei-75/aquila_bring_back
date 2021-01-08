import subprocess
import portion as p
import numpy as np
import gzip
import pysam
import argparse
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser()
parser.add_argument('--chromosome', type=int, nargs='+',
                    help='chromosome')

parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', type=str, help='input directory that contains aquila assembly')
parser.add_argument('--correct_bam_output_dir', type=str, help='output directory for correct bam file')
parser.add_argument('--output_dir', type=str, default="less_covered_regions", help='output directory for assembly and merging results for less covered regions')
parser.add_argument('--merge_supplementary', type=str2bool, default="True", help='whether to merge on supplementary reads')
parser.add_argument('--num_threads', type=int, default="6", help='number of threads for the process')

args = parser.parse_args()
input_dir = args.input_dir
output_dir = args.output_dir
correct_bam_output_dir = args.correct_bam_output_dir
merge_supplementary = args.merge_supplementary
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

# fetch leftmost reads
def find_boundary_reads(bamfile,start,end):
    reads = []
    for read in bamfile.fetch():
        reads += [read]
        
    read_start = []
    read_end = []
    for read in reads:
        if read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start > start -5000 and read.reference_start < end + 5000:
            try:
                read_length = read.reference_end - read.reference_start
            except:
                read_length = 0
            if read_length > 500:
                read_start += [read.reference_start]
                read_end += [read.reference_end]

    leftmost_coords = np.min(read_start)
    rightmost_coords = np.max(read_end)
    for read in reads:
        if read.reference_start == leftmost_coords:
            leftmost_read = read
        if read.reference_end == rightmost_coords:
            rightmost_read = read
    return leftmost_read,rightmost_read

def get_nearby_contigs_aquila(chromo,aquila_bamFile,leftmost_read,rightmost_read):
    # get nearby contigs in aquila
    leftmost_coord = leftmost_read.reference_start
    rightmost_coord = rightmost_read.reference_end
    aquila_left_read = []
    aquila_right_read = []
    if chromo == 23:
        leftmost_reads = aquila_bamFile.fetch("chrX",leftmost_coord-1,leftmost_coord+1)
    else:
        leftmost_reads = aquila_bamFile.fetch('chr' + str(chromo),leftmost_coord-1,leftmost_coord+1)

    if chromo == 23:
        rightmost_reads = aquila_bamFile.fetch("chrX",rightmost_coord-1,rightmost_coord+1)
    else:
        rightmost_reads = aquila_bamFile.fetch('chr' + str(chromo),rightmost_coord-1,rightmost_coord+1)

    for read in leftmost_reads:
        if read.reference_end - read.reference_start > 1000:
            if merge_supplementary:  
                if read.is_secondary or read.is_supplementary:
                    continue
            else:
                if check_read_supplemntary(read,aquila_bamFile,chromo) or read.is_secondary:
                    continue
                
            aquila_left_read += [read]
    for read in rightmost_reads:
        if read.reference_end -  read.reference_start > 1000:
            if merge_supplementary:  
                if read.is_secondary or read.is_supplementary:
                    continue
            else:
                if check_read_supplemntary(read,aquila_bamFile,chromo) or read.is_secondary:
                    continue
            aquila_right_read += [read]
    return aquila_left_read,aquila_right_read

def adjust_clipping(cigar,concat_index):
    clipping = 0
    current_index = 0
    for operation,length in cigar:
        if operation == 4 or operation == 5:
            clipping += length
        current_index += length
        if current_index >= concat_index:
            return clipping

# merge left side
def merge_left_boundary(hp,leftmost_read,aquila_left_read):
    reverse=False
    if len(aquila_left_read) == 0:
        return ""
    elif len(aquila_left_read) == 1:
        concat_index = -1
        aquila_left_read = aquila_left_read[0]
        aquila_end = aquila_left_read.reference_end
        if aquila_left_read.is_supplementary or aquila_left_read.is_secondary:
            return ""
        aligned_pairs = leftmost_read.aligned_pairs
        for i,j in aligned_pairs:
            if j == aquila_end:
                if i != None:
                    concat_index = i
                    
        if concat_index == -1:
            # search in the other direction
            my_start = leftmost_read.reference_start
            aligned_pairs = aquila_left_read.aligned_pairs
            for i,j in aligned_pairs:
                if j == my_start:
                    if i != None:
                        concat_index = i
                        reverse=True
            if concat_index == -1:
                return ""
            
        # check if the sequence is null
        if aquila_left_read.reference_start > leftmost_read.reference_start:
            return ""
        aquila_left_read_seq = aquila_left_read.query
        leftmost_read_seq = leftmost_read.query
        if leftmost_read_seq == None or aquila_left_read_seq == None:
            return ""
        
        if reverse:
            cigar = aquila_left_read.cigar
            clipping = adjust_clipping(cigar,concat_index)
            concat_index = concat_index - clipping
        else:
            cigar = leftmost_read.cigar
            clipping = adjust_clipping(cigar,concat_index)
            concat_index = concat_index - clipping
        
        
#         if leftmost_read_seq == None:
#             leftmost_read_seq = recover_repeat_seq(hp,leftmost_read.query_name)
#         if aquila_left_read_seq == None:
#             aquila_left_read_seq = recover_repeat_seq(hp,aquila_left_read.query_name)
            
        return aquila_left_read_seq,leftmost_read_seq,concat_index, aquila_left_read.query_name, reverse
#         merged_left = aquila_left_read_seq + leftmost_read_seq[concat_index:]
    else:
        return ""

def merge_right_boundary(hp,rightmost_read,aquila_right_read):
    reverse=False
    if len(aquila_right_read) == 0:
        return ""
    elif len(aquila_right_read) == 1:
        concat_index = -1
        aquila_right_read = aquila_right_read[0]
        if aquila_right_read.is_supplementary or aquila_right_read.is_secondary:
            return ""
        aquila_start = aquila_right_read.reference_start
        aligned_pairs = rightmost_read.aligned_pairs
        for i,j in aligned_pairs:
            if j == aquila_start:
                if i != None:
                    concat_index = i
                    
        if concat_index == -1:
            # search in the other direction
            my_end = rightmost_read.reference_end
            aligned_pairs = aquila_right_read.aligned_pairs
            for i,j in aligned_pairs:
                if j == my_end:
                    if i != None:
                        concat_index = i
                        reverse=True
            if concat_index == -1:
                return ""
        
        # check if the sequence is null
        if aquila_right_read.reference_end < rightmost_read.reference_end:
            return ""
        aquila_right_read_seq = aquila_right_read.query
        rightmost_read_seq = rightmost_read.query
        if rightmost_read_seq == None or aquila_right_read_seq == None:
            return ""
        
        if reverse:
            cigar = aquila_right_read.cigar
            clipping = adjust_clipping(cigar,concat_index)
            concat_index = concat_index - clipping
        else:
            cigar = rightmost_read.cigar
            clipping = adjust_clipping(cigar,concat_index)
            concat_index = concat_index - clipping
        
        
#             rightmost_read_seq = recover_repeat_seq(hp,rightmost_read.query_name)
#         if aquila_right_read_seq == None:
#             aquila_right_read_seq = recover_repeat_seq(hp,aquila_right_read.query_name)
            
#         merged_right = rightmost_read_seq[:concat_index] + aquila_right_read_seq
        return aquila_right_read_seq, rightmost_read_seq, concat_index, aquila_right_read.query_name, reverse

    else:
        return ""

def check_read_supplemntary(read,aquila_bamFile,chromo):
    supplementary = False
    # check if it is supplemntary 
    if chromo == 23:
        total_reads = aquila_bamFile.fetch("chrX")
    else:
        total_reads = aquila_bamFile.fetch('chr' + str(chromo))

    for r in total_reads:
        if r.query_name == read.query_name:
            if r.is_supplementary:
                supplementary = True
    
    return supplementary

chromoList = range(1,24)
def concat_with_original_assembly(offset,stride):        
    hpList = ["hp1","hp2"]
    for i in range(offset, len(chromoList), stride):
        chromo = chromoList[i]
        print("concat for chromo %s"%chromo)
        for hp in hpList:
            improved_assemblies = "./%s/improved_assemblies/chr%s_%s_improved_assemblies.txt"%(output_dir,chromo,hp)
            improved_assemblies = open(improved_assemblies).read().split("\n")[:-1]

            for i in range(len(improved_assemblies)):
                region = improved_assemblies[i]
                print(region)

                if i == 0:
                    aquila_bamFile = pysam.AlignmentFile("./%s/Aquila_Contig_chr%s_%s_sorted.bam"%(correct_bam_output_dir,chromo,hp),"rb")
                    aquila_fasta = "./%s/Assembly_Contigs_files/Aquila_Contig_chr%s_%s.fasta"%(input_dir,chromo,hp)
                else:
                    aquila_bamFile = pysam.AlignmentFile("%s/merged_assembly/Aquila_Contig_chr%s_%s_sorted.bam"%(output_dir,chromo,hp))
                    aquila_fasta = "%s/merged_assembly/Aquila_Contig_chr%s_%s.fasta"%(output_dir,chromo,hp)
                aquila_header,aquila_sequence = loadFasta(aquila_fasta)

                aquila_header_dict = {}
                for i in range(len(aquila_header)):
                    h = aquila_header[i]
                    aquila_header_dict[i] = h

                start,end = region.split(",")
                start,end = int(start),int(end)

                bamfile = pysam.AlignmentFile("./%s/chromo%s_assembly/%s_%s_assembly/%s/sorted_output.sam"%(output_dir,chromo,start,end,hp),"rb")
                leftmost_read,rightmost_read = find_boundary_reads(bamfile,start,end)
                aquila_left_read, aquila_right_read = get_nearby_contigs_aquila(chromo,aquila_bamFile,leftmost_read,rightmost_read)

    #             left_supplementary,right_supplementary = check_left_right_supplemntary(aquila_left_read,aquila_right_read,aquila_bamFile,chromo)

    #             if left_supplementary:
    #                 merged_left = ""
    #             else:
    #                 merged_left = merge_left_boundary(hp,leftmost_read,aquila_left_read)
    #             if right_supplementary:
    #                 merged_right = ""
    #             else:
    #                 merged_right = merge_right_boundary(hp,rightmost_read,aquila_right_read)

                headers, sequences = loadFasta("./%s/chromo%s_assembly/%s_%s_assembly/%s/contigs.fasta"%(output_dir,chromo,start,end,hp))

                for i in range(len(headers)):
                    if headers[i] == leftmost_read.query_name:
                        leftmost_index = i

                for i in range(len(headers)):
                    if headers[i] == rightmost_read.query_name:
                        rightmost_index = i

                # perform merging
                left_no_merge,right_no_merge = False,False
                left_queryname,right_queryname = "",""

                merged_left = merge_left_boundary(hp,leftmost_read,aquila_left_read)
                merged_right = merge_right_boundary(hp,rightmost_read,aquila_right_read)

                if len(merged_left) == 0:
                    left_no_merge = True
                else:
                    aquila_left_read_seq,leftmost_read_seq,left_concat_index,left_queryname, left_reverse = merged_left
                    if left_reverse:
                        merged_left = aquila_left_read_seq[:left_concat_index] + leftmost_read_seq
                    else:
                        merged_left = aquila_left_read_seq + leftmost_read_seq[left_concat_index:]
                    sequences[leftmost_index] = merged_left

                if len(merged_right) == 0:
                    right_no_merge = True
                else:
                    aquila_right_read_seq, rightmost_read_seq, right_concat_index,right_queryname, right_reverse  = merged_right
                    if right_reverse:
                        merged_right = rightmost_read_seq + aquila_right_read_seq[right_concat_index:]
                    else:
                        merged_right = rightmost_read_seq[:right_concat_index] + aquila_right_read_seq
                    sequences[rightmost_index] = merged_right

                # adjust the case left and right are the same contig
                if len(merged_left) > 0 and len(merged_right) > 0:
                    if leftmost_read == rightmost_read:
                        if left_reverse == False and right_reverse == False:
                            merged_seq = aquila_left_read_seq + leftmost_read_seq[left_concat_index:]
                            right_concat_index = len(leftmost_read_seq) - right_concat_index
                            merged_seq = merged_seq[:-right_concat_index] + aquila_right_read_seq
                            sequences[leftmost_index] = merged_seq
                            assert leftmost_index == rightmost_index
                        else:
                            if left_reverse:
                                merged_seq = aquila_left_read_seq[:left_concat_index] + leftmost_read_seq
                            else:
                                merged_seq = aquila_left_read_seq +  leftmost_read_seq[left_concat_index:]

                            if right_reverse:
                                merged_seq = merged_seq + aquila_right_read_seq[right_concat_index:]
                            else:
                                right_concat_index = len(leftmost_read_seq) - right_concat_index
                                merged_seq = merged_seq[:-right_concat_index] + aquila_right_read_seq
                            sequences[leftmost_index] = merged_seq
                            assert leftmost_index == rightmost_index


                # replace reads in the original fasta
    #             aquila_need_replaced_reads,aquila_need_replaced_reads_before_supplementary = [],[]
                aquila_need_replaced_reads = []
                if chromo == 23:
                    reads = aquila_bamFile.fetch('chrX',start,end)
                else:
                    reads = aquila_bamFile.fetch('chr' + str(chromo),start,end)

                for read in reads:
                    # ignore secondary reads

                    if read.is_secondary or read.is_supplementary:
                        continue

                    if read.reference_start >= start-500 and read.reference_end <= end+500:
                        aquila_need_replaced_reads += [aquila_header_dict[int(read.query_name)]]

    #             for read in aquila_need_replaced_reads_before_supplementary:
    #                 if check_read_supplemntary(read,aquila_bamFile,chromo):
    #                     continue
    #                 else:
    #                 aquila_need_replaced_reads += [aquila_header_dict[int(read.query_name)]]


                if left_queryname != "" and left_no_merge == False:
                    aquila_need_replaced_reads +=  [aquila_header_dict[int(left_queryname)]]
                if right_queryname != "" and right_no_merge == False:
                    aquila_need_replaced_reads +=  [aquila_header_dict[int(right_queryname)]]
                aquila_need_replaced_reads = list(set(aquila_need_replaced_reads))

                aquila_need_replaced_indices = []
                for i in range(len(aquila_header)):
                    h = aquila_header[i]
                    if h in aquila_need_replaced_reads:
                        aquila_need_replaced_indices += [i]

                for i in sorted(aquila_need_replaced_indices,reverse=True):
                    assert aquila_header[i] in aquila_need_replaced_reads
                    aquila_header.pop(i)
                    aquila_sequence.pop(i)

                new_headers = []
                for i in range(len(headers)):
                    new_headers += ["%s_PS%s:%s_%s"%(i+1,start,end,hp)]

                aquila_header.extend(new_headers)
                aquila_sequence.extend(sequences)

                merged_fasta = ""
                for i in range(len(aquila_header)):
                    merged_fasta += (">" + aquila_header[i] + "\n")
                    merged_fasta += (aquila_sequence[i] + "\n")

                new_aquila_header = []
                for i in range(len(aquila_header)):
                    new_aquila_header += [str(i)]

                rewrite_merged_fasta = ""
                for i in range(len(new_aquila_header)):
                    rewrite_merged_fasta += (">" + new_aquila_header[i] + "\n")
                    rewrite_merged_fasta += (aquila_sequence[i] + "\n")



                with open("./%s/merged_assembly/Aquila_Contig_chr%s_%s.fasta"%(output_dir,chromo,hp), "w") as f:
                    f.write(merged_fasta)

                with open("./%s/merged_assembly/Aquila_Contig_chr%s_%s_rewrite.fasta"%(output_dir,chromo,hp), "w") as f:
                    f.write(rewrite_merged_fasta)

                merged_dir = "./%s/merged_assembly"%output_dir
                contig_file = "%s/Aquila_Contig_chr%s_%s_rewrite.fasta"%(merged_dir,chromo,hp)
                sam_file = "%s/Aquila_Contig_chr%s_%s.bam"%(merged_dir,chromo,hp)
                sorted_sam_file = "%s/Aquila_Contig_chr%s_%s_sorted.bam"%(merged_dir,chromo,hp)

                if chromo == 23:
                    indexFile = "./bam_index/chrX.mmi"
                else:
                    indexFile = "./bam_index/chr%s.mmi"%chromo
                map_cmd = "minimap2 -a %s %s > %s"%(indexFile,contig_file,sam_file)
                sort_cmd = "samtools sort %s -o %s"%(sam_file,sorted_sam_file)
                process = subprocess.Popen(map_cmd,shell=True)
                process.wait()
                
                process = subprocess.Popen(sort_cmd,shell=True)
                process.wait()
                process = subprocess.Popen("samtools index %s"%sorted_sam_file,shell=True)
                process.wait()

process = subprocess.Popen("mkdir ./%s/merged_assembly/"%output_dir,shell=True)
process.wait()

def helper(i):
    return concat_with_original_assembly(i, num_threads)

pool = ThreadPool(num_threads) 
results = pool.map(helper, range(num_threads))
