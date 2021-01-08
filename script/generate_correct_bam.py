import subprocess
import argparse
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool 

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

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', type=str, help='input directory that contains aquila assembly')
parser.add_argument('--correct_bam_output_dir', type=str, help='output directory for correct bam file')
parser.add_argument('--num_threads', type=int, default="6", help='number of threads for the process')

args = parser.parse_args()
correct_bam_output_dir = args.correct_bam_output_dir
input_dir = args.input_dir
num_threads = args.num_threads

process = subprocess.Popen("mkdir %s"%correct_bam_output_dir, shell=True)
process.wait()
chromoList = range(1,24)
for chromo in chromoList:
    hp1_fasta_file = "./%s/Assembly_Contigs_files/Aquila_Contig_chr%s_hp1.fasta"%(input_dir,chromo)
    hp1_headers, hp1_sequences = loadFasta(hp1_fasta_file)

    hp1_header_dict = {}
    for i in range(len(hp1_headers)):
        h = hp1_headers[i]
        hp1_header_dict[h] = i

    hp1_rewrite_headers = [hp1_header_dict[h] for h in hp1_headers]
    hp1_fasta = ''
    for i in range(len(hp1_rewrite_headers)):
        hp1_fasta += '>'
        hp1_fasta += str(hp1_rewrite_headers[i])
        hp1_fasta += '\n'
        hp1_fasta += hp1_sequences[i]
        hp1_fasta += '\n'
        
    hp1_fasta_file = open("./%s/Aquila_Contig_chr%s_hp1.fasta"%(correct_bam_output_dir,chromo),"w")
    hp1_fasta_file.write(hp1_fasta)

    hp2_fasta_file = "./%s/Assembly_Contigs_files/Aquila_Contig_chr%s_hp2.fasta"%(input_dir,chromo)
    hp2_headers, hp2_sequences = loadFasta(hp2_fasta_file)

    hp2_header_dict = {}
    for i in range(len(hp2_headers)):
        h = hp2_headers[i]
        hp2_header_dict[h] = i

    hp2_rewrite_headers = [hp2_header_dict[h] for h in hp2_headers]
    hp2_fasta = ''
    for i in range(len(hp2_rewrite_headers)):
        hp2_fasta += '>'
        hp2_fasta += str(hp2_rewrite_headers[i])
        hp2_fasta += '\n'
        hp2_fasta += hp2_sequences[i]
        hp2_fasta += '\n'
    hp2_fasta_file = open("./%s/Aquila_Contig_chr%s_hp2.fasta"%(correct_bam_output_dir,chromo),"w")
    hp2_fasta_file.write(hp2_fasta)

def generate_bam_file(offset,stride):
    for i in range(offset, len(chromoList), stride):
        chromo = chromoList[i]
        if chromo == 23:
            indexFile = "./bam_index/chrX.mmi"
        else:
            indexFile = "./bam_index/chr%s.mmi"%chromo
        fasta_file = "./%s/Aquila_Contig_chr%s_hp1.fasta"%(correct_bam_output_dir,chromo)
        samFile = "./%s/Aquila_Contig_chr%s_hp1.bam"%(correct_bam_output_dir,chromo)
        minimap_cmd = "minimap2 -a %s %s > %s"%(indexFile,fasta_file,samFile)
        process = subprocess.Popen(minimap_cmd,shell=True)
        process.wait()
        sorted_samFile = "./%s/Aquila_Contig_chr%s_hp1_sorted.bam"%(correct_bam_output_dir,chromo)
        sort_cmd = "samtools sort %s -o %s"%(samFile,sorted_samFile)
        process = subprocess.Popen(sort_cmd,shell=True)
        process.wait()
        index_cmd = "samtools index %s"%sorted_samFile
        process = subprocess.Popen(index_cmd,shell=True)
        process.wait()

        fasta_file = "./%s/Aquila_Contig_chr%s_hp2.fasta"%(correct_bam_output_dir,chromo)
        samFile = "./%s/Aquila_Contig_chr%s_hp2.bam"%(correct_bam_output_dir,chromo)
        minimap_cmd = "minimap2 -a %s %s > %s"%(indexFile,fasta_file,samFile)
        process = subprocess.Popen(minimap_cmd,shell=True)
        process.wait()
        sorted_samFile = "./%s/Aquila_Contig_chr%s_hp2_sorted.bam"%(correct_bam_output_dir,chromo)
        sort_cmd = "samtools sort %s -o %s"%(samFile,sorted_samFile)
        process = subprocess.Popen(sort_cmd,shell=True)
        process.wait()
        index_cmd = "samtools index %s"%sorted_samFile
        process = subprocess.Popen(index_cmd,shell=True)
        process.wait()
        
def helper(i):
    return generate_bam_file(i, num_threads)

pool = ThreadPool(num_threads) 
results = pool.map(helper, range(num_threads))
        
    