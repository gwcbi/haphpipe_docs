### Demo Module

After successful installation, the demo dataset can be run to ensure HAPHPIPE is installed and set up correctly. 


These are the steps completed in the demo:


### Running Demo Automatically
Running the demo is simple and requires a single command:
`hp_demo` or `haphpipe demo`

A specific outdirectory can be specified by:
`hp_demo --outdir $outdir_name`

The output of the entire demo is as such
[put fig here]()

If running the entire demo is not desired, this command can be executed to just pull the references included in HAPHPIPE into the directory that is specified (default outdir is `.`).

`hp_demo --refonly`

Output on the terminal is as such, and the three HIV reference files are located in the subdirectory `refs`. See the [User Guide](https://gwcbi.github.io/haphpipe_docs/install/#reference-files) for more information regarding these reference files.

```
/base/directory/path/of/haphpipe)
Demo was run with --refonly. References are now in outdirectory: $outdir_name/haphpipe_demo/refs.

```
---

#### Output files
If the entire demo is run (i.e., no use of `--refonly`) then these are *all* the files for output **_if_ PredictHaplo is *NOT* installed**. The files are listed for a single sample as example, not all 5 samples. The other files from `multiple_align`, `model_test`, and `build_tree` are all listed.


```
haphpipe_demo
├── SRR8525886
|   ├── SRR8525886_1.fastq
|   ├── SRR8525886_2.fastq
|   └── haphpipe_assemble_02
|       ├── corrected_1.fastq
|       ├── corrected_2.fastq
|       ├── corrected_U.fastq
|       ├── final.bam
|       ├── final.bam.bai
|       ├── final_bt2.out
|       ├── final.fna
|       ├── final.vcf.gz
|       ├── final.vcf.gz.tbi
|       ├── haphpipe.out
|       ├── refined.01.fna
|       ├── refined.02.fna
|       ├── refined.03.fna
|       ├── refined.fna
|       ├── refined_bt2.01.out
|       ├── refined_bt2.02.out
|       ├── refined_bt2.03.out
|       ├── refined_bt2.out
|       ├── refined_summary.out
|       ├── trimmed_1.fastq
|       ├── trimmed_2.fastq
|       ├── trimmed_U.fastq
|       └── trimmomatic_summary.out
|
├── dir_list.txt
├── haphpipe.out
├── hp_multiple_align
|       ├── alignment_region00.fasta
|       ├── alignment_region00.phy
|       ├── alignment_region01.fasta
|       ├── alignment_region01.phy
|       ├── alignment_region02.fasta
|       ├── alignment_region02.phy
|       ├── all_sequences.fasta
|       ├── all_sequences_region00.fasta
|       ├── all_sequences_region01.fasta
|       └── all_sequences_region02.fasta
|
├── alignment_region00_modeltest_results.out
├── alignment_region00_modeltest_results_summary.tsv
├── alignment_region01_modeltest_results.out
├── alignment_region01_modeltest_results_summary.tsv
├── alignment_region02_modeltest_results.out
├── alignment_region02_modeltest_results_summary.tsv
|
├── hp_build_tree
|       ├── RAxML_bestTree.alignment_region00
|       ├── RAxML_bestTree.alignment_region01
|       ├── RAxML_bipartitions.alignment_region00
|       ├── RAxML_bipartitions.alignment_region01
|       ├── RAxML_bipartitionsBranchLabels.alignment_region00
|       ├── RAxML_bipartitionsBranchLabels.alignment_region01
|       ├── RAxML_bootstrap.alignment_region00
|       ├── RAxML_bootstrap.alignment_region01
|       ├── RAxML_info.alignment_region00
|       └── RAxML_info.alignment_region01
|
├── sample#
|   ├── sample#_1.fastq
|   ├── sample#_2.fastq
|   └── haphpipe_assemble_02
....
```

If the entire demo is run (i.e., no use of `--refonly`) then these are *all* the files for output **_if_ PredictHaplo is installed.**


```
haphpipe_demo
├── SRR8525886
|   ├── SRR8525886_1.fastq
|   ├── SRR8525886_2.fastq
|   └── haphpipe_assemble_02
|       ├── corrected_1.fastq
|       ├── corrected_2.fastq
|       ├── corrected_U.fastq
|       ├── final.bam
|       ├── final.bam.bai
|       ├── final_bt2.out
|       ├── final.fna
|       ├── final.vcf.gz
|       ├── final.vcf.gz.tbi
|       ├── haphpipe.out
|       ├── refined.01.fna
|       ├── refined.02.fna
|       ├── refined.03.fna
|       ├── refined.fna
|       ├── refined_bt2.01.out
|       ├── refined_bt2.02.out
|       ├── refined_bt2.03.out
|       ├── refined_bt2.out
|       ├── refined_summary.out
|       ├── trimmed_1.fastq
|       ├── trimmed_2.fastq
|       ├── trimmed_U.fastq
|       ├── trimmomatic_summary.out
|       ├── ph_haplotypes_comb.fna
|       ├── PH01_PRRT
|       |	  ├── PH01_PRRT.best_1_1197.fas
|       |	  ├── PH01_PRRT.best_1_1197.html
|       |	  ├── PH01_PRRT.config.log
|       |	  ├── ph_haplotypes.fna
|       |	  └── ph_summary.txt 
|       ├── PH02_INT
|       |	  ├── PH02_INT.best_1_964.fas
|       |	  ├── PH02_INT.best_1_964.html
|       |	  ├── PH02_INT.config.log
|       |	  ├── ph_haplotypes.fna
|       |	  └── ph_summary.txt 
|       └── PH03_gp120
|       	  ├── PH03_gp120.best_1_1641.fas
|       	  ├── PH03_gp120.best_1_1641.html
|       	  ├── PH03_gp120.config.log
|       	  ├── ph_haplotypes.fna
|       	  └── ph_summary.txt 
|
├── dir_list.txt
├── haphpipe.out
├── hp_multiple_align
|       ├── alignment_region00.fasta
|       ├── alignment_region00.phy
|       ├── alignment_region01.fasta
|       ├── alignment_region01.phy
|       ├── alignment_region02.fasta
|       ├── alignment_region02.phy
|       ├── all_sequences.fasta
|       ├── all_sequences_region00.fasta
|       ├── all_sequences_region01.fasta
|       └── all_sequences_region02.fasta
|
├── alignment_region00_modeltest_results.out
├── alignment_region00_modeltest_results_summary.tsv
├── alignment_region01_modeltest_results.out
├── alignment_region01_modeltest_results_summary.tsv
├── alignment_region02_modeltest_results.out
├── alignment_region02_modeltest_results_summary.tsv
|
├── hp_build_tree
|       ├── RAxML_bestTree.alignment_region00
|       ├── RAxML_bestTree.alignment_region01
|       ├── RAxML_bestTree.alignment_region01
|       ├── RAxML_bipartitions.alignment_region00
|       ├── RAxML_bipartitions.alignment_region01
|       ├── RAxML_bipartitions.alignment_region02
|       ├── RAxML_bipartitionsBranchLabels.alignment_region00
|       ├── RAxML_bipartitionsBranchLabels.alignment_region01
|       ├── RAxML_bipartitionsBranchLabels.alignment_region02
|       ├── RAxML_bootstrap.alignment_region00
|       ├── RAxML_bootstrap.alignment_region01
|       ├── RAxML_bootstrap.alignment_region02
|       ├── RAxML_info.alignment_region00
|       ├── RAxML_info.alignment_region01
|       └── RAxML_info.alignment_region02
|
├── sample#
|   ├── sample#_1.fastq
|   ├── sample#_2.fastq
|   └── haphpipe_assemble_02
....
```

