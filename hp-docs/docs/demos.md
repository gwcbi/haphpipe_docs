Data for both demos is included in the repository. Both samples should run quickly (less than five minutes using eight CPUs). _Note: ensure PredictHaplo is installed prior to running the demos_.<br>

## Demo01

We provide a subsample of reads for this sample. We obtained this sample using [fastq-dump](https://ncbi.github.io/sra-tools/fastq-dump.html).
```
fastq-dump --origfmt --defline-qual "+" --split-files --accession SRR8525886

condo activate haphpipe
hp_sample_reads --fq1 SRR8525886_1.fastq --fq2 SRR8525886_2.fastq --nreads 10000
```

Demo01 uses real HIV-1 data (SRA: [SAMN10820316](https://www.ncbi.nlm.nih.gov/biosample/SAMN10820316)). This demo runs both the _haphpipe_assemble_01_ de novo assembly and _haphpipe_assemble_02_ amplicon assembly pipelines (see [here](https://github.com/gwcbi/haphpipe/wiki/Example-Pipelines) for details). Then, the final output of each pipeline (_final.fna_) is run through the haplotype stages (see [here](https://github.com/gwcbi/haphpipe/wiki/hp_haplotype) for details). To run, activate your HAPHPIPE environment and execute the following code:
```
code block here
```

After running, you should see the following output files:
```
final.fna
final.vcf.gz
final.bam

# For all three amplicons:

ph_summary.txt
ph_haplotypes.fna
```

## Demo02

Demo02 uses a mixture of three HIV-1 subtype B references (Genbank: [AY795905](https://www.ncbi.nlm.nih.gov/nuccore/AY795905), [KP411829](https://www.ncbi.nlm.nih.gov/nuccore/KP411829), [MF373129](https://www.ncbi.nlm.nih.gov/nuccore/MF373129)) in a ratio of 60%, 30%, and 10%, respectively. The provided FASTA file contains the protease and reverse transcriptase gene regions. This demo runs the _hp_assemble_01_ de novo assembly pipeline (see [here](https://github.com/gwcbi/haphpipe/wiki/Example-Pipelines) for details) and uses the final output (_final.fna_) as input for the haplotype stages (see [here](https://github.com/gwcbi/haphpipe/wiki/hp_haplotype) for details).
To run, activate your HAPHPIPE environment and execute the following code:
```
code block here
```
The output for demo02 includes three 1,617bp reconstructed haplotypes and a haplotype diversity estimate of 0.5481. After running, you should see the following output files:
```
final.fna
final.vcf.gz
final.bam
ph_summary.txt
ph_haplotypes.fna
```