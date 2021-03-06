# aquila_bring_back

Bring back discarded reads in regions with low coverage in Aquila for better assembly. 

The scripts have five steps:

1. detect 1-kb regions with low coverage in Aquila (< 90% bases are covered relative to the reference).
2. merge adjacent regions into larger intervals and fetch raw reads in those regions. 
3. partition reads into hp1, hp2 and discarded reads
4. concat discarded reads with hp1 or hp2 reads and perform local assembly 
5. replace the original assembly in aquila with the new assembly including discarded reads if the new assembly is better


Usage 

```bash
create_bam.sh
python ./script/main.py --input_dir=L5_NA24385 --correct_bam_output_dir=correct_bam_files \
--output_dir=less_covered_regions --merge_supplementary=True --num_threads=8

```


Argument:

 - ` --input_dir:` the directory that contains the original aquila assembly and reads
 - `--correct_bam_output_dir:` the directory to output the correct sorted bam files for aquila (the correct bam files are generated by renaming the contig headers)
 - `--output_dir:` the output directory for the aquila assembly after merging with new assembly including discarded reads
 - `--merge_supplementary:` whether to merge contigs with supplementary mappings


Outputs:

The final fasta files after meging with new assembly are in the $output\_dir/merged_assembly/ Aquila\_Contig\_chr\*\_hp\*.fasta

```bash
│── chromo*_less_covered_regions.txt
├── chromo*_read_bams
├── chromo*_phased_reads
├── chromo*_concat_reads
├── chromo*_assembly
├── improved_assemblies
└── merged_assembly (final outputs)
  	 ├── Aquila_Contig_chr*_hp*.fasta 
  	 ├── Aquila_Contig_chr*_hp*_sorted.bam 
  	 ├── Aquila_Contig_chr*_hp*_rewrite.fasta 
     	 (intermediate fasta files used to generate the correct bam files by renaming the headers)
```





