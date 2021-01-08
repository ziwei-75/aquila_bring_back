# aquila_bring_back

Bring back discarded reads in regions with low coverage in Aquila for better assembly. 

Usage 

```python
python ./script/main.py --input_dir=L5_NA24385 --correct_bam_output_dir=correct_bam_files --output_dir=less_covered_regions --merge_supplementary=True --num_threads=8

```

The scripts have Five steps:

1. detect 1-kb regions with low coverage in Aquila (< 90% bases are covered relative to the reference).
2. merged adjacent regions into larger intervals and fetch raw reads in those regions. 
3. partition reads into hp1, hp2 and discarded reads
4. concat discarded reads with hp1 or hp2 reads and perform local assembly 
5. replace the original assembly in aquila with the new assembly using discarded reads if the new assembly is better





