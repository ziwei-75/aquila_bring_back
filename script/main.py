import argparse
import subprocess

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

parser.add_argument('--input_dir', type=str, help='input directory that contains aquila assembly')
parser.add_argument('--correct_bam_output_dir', type=str, help='output directory for correct bam file')
parser.add_argument('--output_dir', type=str, default="less_covered_regions", help='output directory for assembly and merging results for less covered regions')
parser.add_argument('--merge_supplementary', type=str2bool, default="True", help='whether to merge on supplementary reads')
parser.add_argument('--num_threads', type=int, default="6", help='number of threads for the process')

args = parser.parse_args()
correct_bam_output_dir = args.correct_bam_output_dir
output_dir = args.output_dir
input_dir = args.input_dir
merge_supplementary = args.merge_supplementary
num_threads = args.num_threads

if __name__ == "__main__": 
    p1 = subprocess.Popen("python ./script/generate_correct_bam.py --input_dir=%s --correct_bam_output_dir=%s --num_threads=%s"%(input_dir,correct_bam_output_dir,num_threads),shell=True)
    p1.wait()
    print("finish generating correct bam files by renaming the headers")
    
    p2 = subprocess.Popen("python ./script/detect_low_coverage_regions.py --input_dir=%s --output_dir=%s --correct_bam_output_dir=%s  --num_threads=%s"%(input_dir,output_dir,correct_bam_output_dir,num_threads),shell=True)
    p2.wait()
    print("finish detecting less covered regions")
    
    p3 = subprocess.Popen("python ./script/merge_region_and_generate_bams.py --input_dir=%s --output_dir=%s --num_threads=%s"%(input_dir,output_dir,num_threads),shell=True)
    p3.wait()
    print("finish generating bams")
    
    p4 = subprocess.Popen("python ./script/partition_reads.py --input_dir=%s --output_dir=%s --num_threads=%s"%(input_dir,output_dir,num_threads),shell=True)
    p4.wait()
    print("finish partitioning reads")
    
    p5 = subprocess.Popen("python ./script/concat_and_assembly.py --output_dir=%s --num_threads=%s"%(output_dir,num_threads),shell=True)
    p5.wait()
    print("finish assembly")
    
    p6 = subprocess.Popen("python ./script/decide_better_assembly.py --output_dir=%s --correct_bam_output_dir=%s --num_threads=%s"%(output_dir,correct_bam_output_dir,num_threads),shell=True)
    p6.wait()
    print("finish deciding")
    
    p7 = subprocess.Popen("python ./script/concat_new_assembly_with_original_aquila.py --input_dir=%s --correct_bam_output_dir=%s --output_dir=%s --merge_supplementary=%s --num_threads=%s"%(input_dir,correct_bam_output_dir,output_dir,str(merge_supplementary),num_threads),shell=True)
    p7.wait()
