# Haphpipe User Guide
_**HA**plotype and **PH**ylodynamics pipeline for viral assembly, population genetics, and phylodynamics._


We assume basic familiarity with conda environments, basic bash knowledge and knowledge of general next-generation sequencing (NGS) concepts. For more information regarding all of these, see helpful links and FAQ. **(insert links to markdown pages once these are split up)** Remember that a directory is a simply a folder nesting system, similar to what you see on your computer. Where you see "folder" below you can also say "directory."

![directory structre](/Users/keyliegibson/Desktop/Working_projects/haphpipe/Things to add to github/593d447b94e14c2d4418f0f58a604169.png)

HAPHPIPE is intended only for Linux and Mac OS X platforms. 

This User Guide was developed by undergrads and graduate students to be accessible for users at all stages. Maggie Steiner, Keylie M Gibson, Matthew L Bendall, .... all contributed to developing this User Guide and testing HAPHPIPE.

See **paper** for more information and **paper** for validation study.
Citing HAPHPIPE:

```
insert eventually
```

***I want each of these a separate document, for now I'll just make a giant one.***



---

# Installation
For now, the installation goes like this. We are hoping to make it into a stand alone bioconda recipe soon.

HAPHIPE depends on more than a dozen different programs, each of which may itself depend on other programs and libraries. Installing everything separately would be a nightmare, so you should use the package manager "Bioconda" to install everything. Bioconda is a popular package manager in the bioinformatics world. See helpful links for more information and resources for Bioconda. This User Guide describes where to get Bioconda and how to install it, then how to use Bioconda to install the necessary programs for HAPHPIPE. This User Guide details the acquisition and installation of one program, GTAK, that is not handled by Bioconda, and finally the acquisition and installation of HAPHPIPE itself.

Here, we describe the procedure for installing HAPHPIPE using the package manager Bioconda (Grüning et al. 2018) on the command line. If you are unfamiliar with Bioconda, please see https://bioconda.github.io for installation and channel setup. HAPHPIPE is available at https://github.com/gwcbi/haphpipe and is written in Python 3.7.2 coding language. An online documentation is available at https://github.com/gwcbi/haphpipe/wiki. The installation process begins with the creation of a conda environment that installs the necessary dependencies required by HAPHPIPE. Once the conda environment has been created, it can be activated for use with the command conda activate haphpipe. Due to license restrictions, Bioconda cannot distribute and install GATK (McKenna et al. 2010; Van der Auwera et al. 2013; Poplin et al. 2018) directly. To fully install GATK, you must download a licensed copy of GATK (version 3.8-0) from the Broad Institute (https://software.broadinstitute.org/gatk/download/archive). You can then install HAPHPIPE using the single command pip install git+git://github.com/gwcbi/haphpipe.git, which pulls the repository from Github.


<need homebrew then wget>

__1. Install [conda](https://bioconda.github.io/user/install.html#set-up-channels)__

Download the conda package:

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh`

__2. Set up [conda channels](https://bioconda.github.io/user/install.html#set-up-channels)__

These need to be put in as a command in the order they come, as this sets the priority for where packages are pulled from. Therefore, in this order, conda-forge is top priority and explored first when downloading a program.

```
conda config --add channels R
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```


__3. Create a conda environment with the following dependencies__

```bash
conda create -n haphpipe \
    python \
    future \
    pyyaml \
    biopython \
    seqtk \
    bowtie2 \
    bwa \
    flash \
    freebayes \
    mummer \
    picard \
    trimmomatic \
    samtools=1.9 \
    gatk=3.8 \
    spades \
    blast \
    sierrapy

```

__4. Activate the environment__

```
conda activate haphpipe
```

__5. Install GATK__

GATK can be dowloaded before or after creating the haphpipe conda environment, but it *must be registered after* creating the conda environment.

Due to license restrictions, bioconda cannot distribute
and install GATK directly. To fully install GATK, you must
download a licensed copy of GATK (version 3.8-0) from the Broad Institute:
[https://software.broadinstitute.org/gatk/download/archive](https://software.broadinstitute.org/gatk/download/archive).

Once you download a copy, *register the package using gatk3-register*:

```
gatk3-register /path/to/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2
```

This will copy GATK into your conda environment.

NOTE: HAPHPIPE was developed and tested using GATK 3.8.

__6. Install HAPHPIPE__

```
pip install git+git://github.com/gwcbi/haphpipe.git
```

Upon completion of the installation, you can test it to ensure the repository has been installed completely and correctly by running haphpipe -h (in quick start). Once HAPHPIPE is installed and performing correctly, there is no need to install it again; simply activate the conda environment when needed by executing `conda activate haphpipe`. If a new version is released in the future, HAPHPIPE can be updated by running the command line `pip install --upgrade git+git://github.com/gwcbi/haphpipe.git`. At any point, the `-h` option that follows any HAPHPIPE stage will output a help message that provides a description of the stage and the desired input(s) and output(s).

__7. Installing PredictHaplo__
The user is required to install PredictHaplo on their own, as there are system-dependent variables within the installation process. We were unable to successfully incorporate PredictHaplo into HAPHPIPE so this step was unnecessary for you to install.

Users are required to download PredictHaplo on their own computing system prior to running any of the haplotype stages (_hp_predict_haplo_ and _hp_ph_parser_).

Here is how the GW CBI team installed PredictHaplo onto our HPC, which has a slurm scheduling system and uses Lua module files.

We cannot help with the installation of this software, but have provided the code that we used here to install PredictHaplo onto our system. Please see [their website](https://bmda.dmi.unibas.ch/software.html) for contact information if you need help installing their program.

You can download the new version for paired-end reads of the PredictHaplo program from [their website](https://bmda.dmi.unibas.ch/software.html). The instructions are in the README in the folder downloaded. We have also provided them here for your access. **provide the readme file for them**


This module loads predicthaplo onto GWU's HPC - Colonial One. See https://bmda.dmi.unibas.ch/software.html

```bash
cd /path/to/modules/predicthaplo

# use gcc 4.9.4, add blas and lapack to library path
module load gcc/4.9.4
module load blas/gcc/64
module load lapack/gcc/64

# download source
cd archive
wget https://bmda.dmi.unibas.ch/software/PredictHaplo-Paired-0.4.tgz
cd ..

# unzip source and change directory name
tar xfvz archive/PredictHaplo-Paired-0.4.tgz
mv PredictHaplo-Paired-0.4 0.4
cd 0.4

# install scythestat
tar xfvz scythestat-1.0.3.tar.gz
cd scythestat-1.0.3
./configure --prefix=/path/to/modules/predicthaplo/0.4/NEWSCYTHE
make install
cd ..

# compile predicthaplo
make
```


If a `segfault` error occurs during the `hp_predict_haplo` stage, this is **not** a characteristic of HAPHPIPE but rather that of PredictHaplo. Sometimes, we have luck if we just rerun the code again or move to an interactive CPU node. We are unsure what causes this error, and we only see it between the local and global reconstruction phases in PredictHaplo.



</br>

---
---

# Quick start
##	Activate HAPHPIPE

__1. Activate haphpipe__

Make sure you have conda running.
For students at GW using Colonial One, you need to load the `miniconda3` module like such prior to activating the haphpipe conda environemnt: 
`module load miniconda3`.

```
conda activate haphpipe
```

**2. Test that it is loaded correctly**

```
haphpipe -h
```

should produce:

```
Program: haphpipe (haplotype and phylodynamics pipeline)
Version: 0.8.1

Commands:
 -- Reads
    sample_reads             subsample reads using seqtk
    trim_reads               trim reads using Trimmomatic
    join_reads               join reads using FLASh
    ec_reads                 error correct reads using SPAdes

 -- Assemble
    assemble_denovo          assemble reads denovo
    assemble_amplicons       assemble contigs to amplicon regions
    assemble_scaffold        assemble contigs to genome
    align_reads              align reads to reference
    call_variants            call variants
    vcf_to_consensus         create consensus sequence from VCF
    refine_assembly          iterative refinement: align - variants - consensus
    finalize_assembly        finalize consensus sequence

 -- Haplotype
    predict_haplo            assemble haplotypes with PredictHaplo
    ph_parser                parse output from PredictHaplo.

 -- Annotate
    pairwise_align           align consensus to an annotated reference
    extract_pairwise         extract sequence regions from pairwise alignment
    annotate_from_ref        annotate consensus from reference annotation

 -- Miscellaneous
    demo                     setup demo directory and test data
```


##	Demo
Demo1, de novo assembly, tell them exactly what fastq read files they need to download, exactly where to find those files, and exactly where to download them.  Provide URLs and what to expect to see in the SRL archives.  These are beginners, they don't know how to find anything.  Then tell them that Demo2 you will lead them through a reference guided assembly and that they will need to download the following finished genomes sequences and exactly where to get them

We will now go over two ways to complete these demos: (i) automatically, which is to ultimately test if HAPHPIPE has been installed correctly and (ii) interactively, to get a feel and experience running HAPHPIPE yourself.

Data for both demos is included in the repository. Both samples should run quickly (less than five minutes using eight CPUs). Note: ensure PredictHaplo is installed prior to running the demos.

1. Automatically to test the installation of HAPHPIPE

`insert the code to run the demo data`

2. Interactively to familiarize yourself with running HAPHPIPE



##	Example pipelines
The example pipelines are written in bash scripting language. The reference files used in both examples are included in the demo data. To run in haphpipe, execute one of the following lines:

Pipeline 1 implements amplicon assembly using a *de novo* approach. Reads are
error-corrected and used to refine the initial assembly, with up to 5
refinement steps.

Pipeline 2 implements amplicon assembly using a reference-based mapping 
approach. Reads are error-corrected and used to refine the initial assembly,
with up to 5 refinement steps.


### Pipeline 1: __`haphpipe_assemble_01`__

This pipeline implements *de novo* assembly. Reads are first trimmed (*trim_reads*) and used as input for denovo assembly (*assemble_denovo*). The *de novo* assembly stage automatically performs error correction on the trimmed reads. The assembled contigs are used as input for amplicon assembly (*assemble_amplicons*) along with reference FASTA and GTF files. The assembly is then iteratively refined up to five times (*refine_assembly*) by mapping corrected reads to the assembled FASTA file and lastly finalized (*finalize_assembly*), resulting in a FASTA file with final consensus sequences, final VCF, and aligned BAM file.


![haphpipe01](/Users/keyliegibson/Desktop/Working_projects/haphpipe/haphpipe.pdf)

To see the input information for Pipeline 1, use the `-h` option again like so:
`haphpipe_assemble_01 -h`, and it will show the output:


```
USAGE:
haphpipe_assemble_01 [read1] [read2] [reference_fasta] [reference_gtf] [samp_id] <outdir>

----- HAPHPIPE assembly pipeline 01 -----

This pipeline implements amplicon assembly using a denovo approach. Reads are
error-corrected and used to refine the initial assembly, with up to 5
refinement steps.

Input:
read1:             Fastq file for read 1. May be compressed (.gz)
read2:             Fastq file for read 2. May be compressed (.gz)
reference_fasta:   Reference sequence (fasta)
reference_gtf:     Amplicon regions (GTF)
samp_id:           Sample ID
outdir:            Output directory (default is [sample_dir]/haphpipe_assemble_01)
```

General command to execute pipeline 1:
```
haphpipe_assemble_01 samp/read1.fq.gz samp/read2.fq.gz refs/ref.fasta refs/ref.gtf samp
```

Example command to run with demo01:
```
insert
```


### Pipeline 2: __`haphpipe_assemble_02`__

This pipeline implements reference-based mapping assembly. Reads are first trimmed (*trim_reads*) and error-corrected (*ec_reads*). The corrected reads are used as input for reference-based mapping assembly (*refine_assembly*) for up to five iterations. Lastly, the assembly is finalized (*finalize_assembly*) by mapping reads onto the refined reference sequence. The final output is a FASTA file with final consensus sequences, final VCF, and aligned BAM file.

**(insert Figure 2 here)**

To see the input information for Pipeline 1, use the `-h` option again like so:
`haphpipe_assemble_02 -h`, and it will show the output:

```
USAGE:
haphpipe_assemble_02 [read1] [read2] [amplicons_fasta] [samp_id] <outdir>

----- HAPHPIPE assembly pipeline 02 -----

This pipeline implements amplicon assembly using a reference-based approach.
Reads are error-corrected and aligned to provided amplicon reference with up to
five refinement steps.

Input:
read1:             Fastq file for read 1. May be compressed (.gz)
read2:             Fastq file for read 2. May be compressed (.gz)
amplicons_fasta:   Amplicon reference sequence (fasta)
samp_id:           Sample ID
outdir:            Output directory (default is sample_dir/haphpipe_assemble_02)
```

General command to execute pipeline 1:
```
haphpipe_assemble_02 samp/read1.fq.gz samp/read2.fq.gz refs/ref.fasta samp
```

Example command to run with demo01:
```
insert
```

</br>

---
---

# HAPHPIPE suite
Each stage can be run on its own. Stages are grouped into 4 categories: `hp_reads`, `hp_assemble`, `hp_haplotype`, and `hp_annotate`.
More detailed description of command line options for each stage are available in the [wiki](https://github.com/gwcbi/haphpipe/wiki). To view all available stages in HAPHPIPE, run: 

```
haphpipe -h
```

Output will look like:

```
Program: haphpipe (haplotype and phylodynamics pipeline)
Version: 0.8.1

Commands:
 -- Reads
    sample_reads             subsample reads using seqtk
    trim_reads               trim reads using Trimmomatic
    join_reads               join reads using FLASh
    ec_reads                 error correct reads using SPAdes

 -- Assemble
    assemble_denovo          assemble reads denovo
    assemble_amplicons       assemble contigs to amplicon regions
    assemble_scaffold        assemble contigs to genome
    align_reads              align reads to reference
    call_variants            call variants
    vcf_to_consensus         create consensus sequence from VCF
    refine_assembly          iterative refinement: align - variants - consensus
    finalize_assembly        finalize consensus sequence

 -- Haplotype
    predict_haplo            assemble haplotypes with PredictHaplo
    ph_parser                parse output from PredictHaplo.

 -- Annotate
    pairwise_align           align consensus to an annotated reference
    extract_pairwise         extract sequence regions from pairwise alignment
    annotate_from_ref        annotate consensus from reference annotation

 -- Miscellaneous
    demo                     setup demo directory and test data

```

HAPHPIPE consists of a suite of sub-commands under each stage that are invoked as follows:

`haphpipe [stage] [sub-command] [options]`

For example, to join paired end reads, one would invoke the following:

`haphpipe join_reads --fq1 trimmed_1.fastq --fq2 trimmed_2.fastq`


I want all stages to have a table like this detailing the options like this.
![example](/Users/keyliegibson/Downloads/Screen Shot 2019-09-28 at 5.51.07 PM.png)

##	hp_reads

*hp_reads* involves cleaning up the raw read sequences, as well as other processing steps. Modules to manipulate reads. 

### *sample_reads*
Subsample reads using seqtk ([documentation](https://github.com/lh3/seqtk)). Input is reads in FASTQ format. Output is sampled reads in FASTQ format. You do not have to have all read options (i.e., read1, read2 AND unpaired reads). You can have a combination of any of those.

#### Usage and option summary

**Usage:**

`haphpipe sample_reads [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> [--outdir]`

**(or):**

`hp_sample_reads [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> [--outdir]`

*Output files:* <br> `sample_1.fastq` <br> `sample_2.fastq`

*Input/Output Arguments:* 

Option   | Description 
---------|-------------
--fq1    | Fastq file with read 1.
--fq2    | Fastq file with read 2.
--fqU    | Fastq file with unpaired reads.
--outdir | Output directory (default: current directory).

*Options:*

Option   | Description |
---------|-------------|
--nreads | Number of reads to sample. If greater than the number of reads in file, all reads will be sampled. |
--frac   | Fraction of reads to sample, between 0 and 1. Each read has [frac] |probability of being sampled, so the number of sampled reads is not precisely [frac * number of reads]. |
--seed   | Seed for random number generator.|

*Settings:*

Option    | Description |
----------|-------------|
--quiet   | Do not write output to console (silence stdout and stderr), default is False. |
--logfile | Append console output to this file. |
--debug   | Print commands but do not run, default is False. |

#### Example usage
This pulls 1000 reads from these paired end files with a starting seed of 1234.
```
haphpipe sample_reads --fq1 read_1.fastq --fq2 read_2.fastq --nreads 1000 --seed 1234
```

--
	
### *trim_reads*
Trim reads using Trimmomatic ([documentation](http://www.usadellab.org/cms/?page=trimmomatic)). Input is reads in FASTQ format. Output is trimmed reads in FASTQ format.

#### Usage and option summary

**Usage:**

`haphpipe trim_reads [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> [--outdir]`

**(or):**

`hp_trim_reads [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> [--outdir]`

*Output files:* <br> 
`trimmed_1.fastq` <br> 
`trimmed_2.fastq` <br> 
`trimmed_U.fastq` <br> 
`trimmomatic_summary.out`

*Input/Output Arguments:* 

Option    | Description
----------|-------------
--fq1     | Fastq file with read 1.
--fq2     | Fastq file with read 2.
--fqU     | Fastq file with unpaired reads.
--outdir  | Output directory (default: current directory).

*Options:*

Option         | Description
---------------|-------------
--adapter_file | Adapter file.
--trimmers     | Trim commands for trimmomatic (default: ['LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36']).
--encoding     | Quality score encoding.

*Settings:*

Option     | Description
-----------|-------------
--ncpu     | Number of CPU to use (default: 1).
--quiet    | Do not write output to console (silence stdout and stderr) (default: False).
--logfile  | Append console output to this file.
--debug    | Print commands but do not run (default: False).

#### Example usage
This trims paired end read files 1 and 2.

```
haphpipe trim_reads --fq1 read_1.fastq --fq2 read_2.fastq 
```

--

### *join_reads*
Join reads using FLASH ([paper](https://www.ncbi.nlm.nih.gov/pubmed/21903629)). Input is reads in FASTQ format. Output is joined reads in FASTQ format.

#### Usage and option summary

**Usage:**

`haphpipe join_reads [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> [--outdir]`

**(or):**

`hp_join_reads [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> [--outdir]`

*Output files:* <br> `joined.fastq` <br> `notjoined_1.fastq` <br> `notjoined_2.fastq`

*Input/Output Arguments:* 

Option   | Description
---------|-------------
--fq1    | Fastq file with read 1.
--fq2    | Fastq file with read 2.
--outdir | Output directory (default: current directory).

*Settings:*

Option          | Description
----------------|-------------
--min_overlap   | The minimum required overlap length between two reads to provide a confident overlap (default: 10).
--max_overlap   | Maximum overlap length expected in approximately 90% of read pairs, longer overlaps are penalized.
--allow_outies  | Also try combining read pairs in the "outie" orientation (default: False).
--encoding      | Quality score encoding.

*Settings:*

Option     | Description
-----------|-------------
--ncpu     | Number of CPU to use (default: 1).
--keep_tmp | Keep temporary directory (default: False).
--quiet    | Do not write output to console (silence stdout and stderr) (default: False).
--logfile  | Append console output to this file.
--debug    | Print commands but do not run (default: False).

#### Example usage

--

### *ec_reads*
Error correction using SPAdes ([documentation](http://cab.spbu.ru/software/spades/)). Input is reads in FASTQ format. Output is error-corrected reads in FASTQ format. Remember that HAPHPIPE is intended for Illumina reads, therefore the error correction is based on Illumina sequencing errors.

#### Usage and option summary

**Usage:**

`haphpipe ec_reads [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> [--outdir]`

**(or):**

`hp_ec_reads [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> [--outdir]`

*Output files:* <br> `corrected_1.fastq` <br> `corrected_2.fastq` <br> `corrected_U.fastq`

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         | Fastq file with read 1.
--fq2         | Fastq file with read 2.
--fqU         | Fastq file with unpaired reads.
--outdir      | Output directory (default: current directory).

*Settings:*

Option        | Description
------------- |-------------
--ncpu        | Number of CPU to use (default: 1).
--keep_tmp    | Keep temporary directory (default: False).
--quiet       | Do not write output to console (silence stdout and stderr) (default: False).
--logfile     | Append console output to this file.
--debug       | Print commands but do not run, default is False.

#### Example usage

--
	
##	hp_assemble
Assemble consensus sequence(s). Input reads (in FASTQ format) are assembled 
using either denovo assembly or reference-based alignment. 
Resulting consensus can be further refined.

### *assemble_denovo*
Assemble reads via de novo assembly using SPAdes ([documentation](http://cab.spbu.ru/software/spades/)). Input is reads in FASTQ format. Output is contigs in FNA format.

#### Usage and option summary

**Usage:**

`haphpipe assemble_denovo [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> [--outdir]`

**(or):**

`hp_assemble_denovo [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> [--outdir]`

*Output files:* <br> `denovo_contigs.fna` <br> `denovo_summary.txt`

*Input/Output Arguments:* 

Option      | Description
------------|-------------
--fq1       | Fastq file with read 1.
--fq2       | Fastq file with read 2.
--fqU       | Fastq file with unpaired reads.
--outdir    | Output directory (default: current directory).

*Options:*

Option                | Description
----------------------|-------------
--no_error_correction | Do not perform error correction (default: False)
--subsample           | Use a subsample of reads for assembly
--seed                | Seed for random number generator (ignored if not subsampling)

*Settings:*

Option      | Description
------------|-------------
--ncpu      | Number of CPU to use (default: 1).
--keep_tmp  | Keep temporary directory (default: False).
--quiet     | Do not write output to console (silence stdout and stderr) (default: False).
--logfile   | Append console output to this file.
--debug     | Print commands but do not run (default: False).

#### Example usage


### *assemble_amplicons*
Assemble contigs from de novo assembly using both a reference sequence and amplicon regions with MUMMER 3+ ([documentation](http://mummer.sourceforge.net/manual/)). Input is contigs and reference sequence in FASTA format and amplicon regions in GTF format.

#### Usage and option summary

**Usage:**

`haphpipe assemble_amplicons [OPTIONS] [SETTINGS] --contigs_fa <FASTA> --ref_fa <FASTA> --ref_gtf <GTF> [--outdir]`

**(or):**

`hp_assemble_amplicons [OPTIONS] [SETTINGS] --contigs_fa <FASTA> --ref_fa <FASTA> --ref_gtf <GTF> [--outdir]`

*Output files:* <br> `amplicon_assembly.fna`

*Input/Output Arguments:* 

Option        | Description
--------------|-------------
--contigs_fa  | Fasta file with assembled contigs.
--ref_fa      | Fasta file with reference genome to scaffold against.
--ref_gtf     | GTF format file containing amplicon regions.
--outdir      | Output directory (default: current directory).

*Scaffold Options:*

Option      | Description
------------|-------------
--sample_id | Sample ID (default: sampleXX).
--padding   | Bases to include outside reference annotation (default: 50).

*Settings:*

Option      | Description
------------|-------------
--keep_tmp  | Keep temporary directory (default: False).
--quiet     | Do not write output to console (silence stdout and stderr) (default: False).
--logfile   | Append console output to this file.
--debug     | Print commands but do not run (default: False).


#### Example usage

### *assemble_scaffold*
Scaffold contigs against a reference sequence with MUMMER 3+ ([documentation](http://mummer.sourceforge.net/manual/)). Input is contigs in FASTA format and reference sequence in FASTA format. Output is scaffold assembly, alligned scaffold, imputed scaffold, and padded scaffold in FASTA format.

#### Usage and option summary

**Usage:**

`haphpipe assemble_scaffold [OPTIONS] [SETTINGS] --contigs_fa <FASTA> --ref_fa <FASTA> [--outdir]`

**(or):**

`hp_assemble_scaffold [OPTIONS] [SETTINGS] --contigs_fa <FASTA> --ref_fa <FASTA> [--outdir]`

*Output files:* <br> `scaffold_aligned.fa` <br> `scaffold_assembly.fa` <br> `scaffold_imputed.fa` <br> `scaffold_padded.out`

*Input/Output Arguments:* 

Option        | Description
--------------|-------------
--contigs_fa  | Fasta file with assembled contigs.
--ref_fa      | Fasta file with reference genome to scaffold against.
--outdir      | Output directory (default: current directory).

*Options:*

Option      | Description
------------|-------------
--seqname   | Name to append to scaffold sequence (default: sample01).

*Settings:*

Option      | Description
------------|-------------
--keep_tmp  | Additional options (default: False).
--quiet     | Do not write output to console (silence stdout and stderr) (default: False).
--logfile   | Append console output to this file.
--debug     | Print commands but do not run (default: False).

#### Example usage

### *align_reads*
Map reads to reference sequence (instead of running de novo assembly) using Bowtie2 ([documentation](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)) and Picard ([documentation](https://broadinstitute.github.io/picard/)). Input is reads in FASTQ format and reference sequence in FASTA format. 

#### Usage and option summary

**Usage:**

`haphpipe align_reads [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> --ref_fa <FASTA> [--outdir]`

**(or):**

`hp_align_reads [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> --ref_fa <FASTA> [--outdir]`

*Output files:* <br> *aligned.bam* <br> *aligned.bt2.out*

*Input/Output Arguments:* 

Option      | Description
------------|-------------
--fq1       | Fastq file with read 1.
--fq2       | Fastq file with read 2.
--fqU       | Fastq file with unpaired reads.
--ref_fa    | Reference fasta file.
--outdir    | Output directory (default: current directory).

*Options:*

Option              | Description
--------------------|-------------
--bt2_preset        | {very-fast, fast, sensitive,very-sensitive,very-fast-local,fast-local,sensitive-local,very-sensitive-local}
--sample_id         | Sample ID. Used as read group ID in BAM (default: sampleXX).
--no_realign        | Do not realign indels (default: False).
--remove_duplicates | Remove duplicates from final alignment. Otherwise duplicates are marked but not removed (default: False).
--encoding          | {Phred+33,Phred+64} Quality score encoding.

*Settings:*

Option      | Description
------------|-------------
--ncpu      | Number of CPUs to use (default: 1).
--xmx       | Maximum heap size for Java VM, in GB (default: 32).
--keep_tmp  | Additional options (default: False).
--quiet     | Do not write output to console (silence stdout and stderr) (default: False).
--logfile   | Append console output to this file.
--debug     | Print commands but do not run (default: False).

#### Example usage

### *call_variants*
Variant calling from alignment using GATK ([documentation](https://software.broadinstitute.org/gatk/download/archive)). Input is alignment file in BAM format and reference sequence in FASTA format (either reference from reference-based assembly or consensus final sequence from de novo assembly). Output is a Variant Call File (VCF) format file. 

#### Usage and option summary

**Usage:**

`haphpipe call_variants [OPTIONS] [SETTINGS] --aln_bam <BAM> --ref_fa <FASTA> [--outdir]`

**(or):**

`hp_call_variants [OPTIONS] [SETTINGS] --aln_bam <BAM> --ref_fa <FASTA> [--outdir]`

*Output files:* <br> *variants.vcf.gz*

*Input/Output Arguments:* 

Option    | Description
----------|-------------
--aln_bam | Alignment file.
--ref_fa  | Reference fasta file.
--outdir  | Output directory (default: False).

*Options:*

Option          | Description
----------------|-------------
--emit_all      | Output calls for all site (default: False).
--min_base_qual | Minimum base quality required to consider a base for calling (default: 15).

*Settings:*

Option       | Description
-------------|-------------
--ncpu       | Number of CPUs to use (default: 1).
--xmx        | Maximum heap size for Java VM, in GB (default: 32).
--keep_tmp   | Additional options (default: False).
--quiet      | Do not write output to console (silence stdout and stderr) (default: False).
--logfile    | Append console output to this file.
--debug      | Print commands but do not run (default: False).

#### Example usage

### *vcf_to_consensus*
Generate a consensus sequence from a VCF file. Input is a VCF file. Output is the consensus sequence in FASTA format. 

#### Usage and option summary

**Usage:**

`haphpipe vcf_to_consensus [OPTIONS] [SETTINGS] --vcf <FASTQ> [--outdir] [--sampidx]`

**(or):**

`hp_vcf_to_consensus [OPTIONS] [SETTINGS] --vcf <FASTQ> [--outdir] [--sampidx]`

*Output files:* <br> *consensus.fna*

*Input/Output Arguments:* 

Option    | Description
----------|-------------
--vcf     | VCF file (created with all sites).
--outdir  | Output directory (default: False).
--sampidx | Index for sample if multi-sample VCF (default: 0).

*Options:*

Option    | Description
----------|-------------
--min_DP  | Minimum depth to call site (default: 1).
--major   | Allele fraction to make unambiguous call (default: 0.5).
--minor   | Allele fraction to make ambiguous call (default: 0.2).

*Settings:*

Option       | Description
-------------|-------------
--keep_tmp   | Additional options (default: False).
--quiet      | Do not write output to console (silence stdout and stderr) (default: False).
--logfile    | Append console output to this file.


#### Example usage

### *refine_assembly*
Map reads to a denovo assembly or reference alignment. Assembly or alignment is iteratively updated. Input is reads in FASTQ format and reference sequence (assembly or reference alignment) in FASTA format. Output is refined assembly in FASTA format.

#### Usage and option summary

**Usage:**

`haphpipe refine_assembly [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> --ref_fa <FASTA> [--outdir]`

**(or):**

`hp_refine_assembly [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> --ref_fa <FASTA> [--outdir]`

*Output files:* <br> *refined.fna*

*Input/Output Arguments:* 

Option     | Description
-----------|-------------
--fq1      | Fastq file with read 1.
--fq2      | Fastq file with read 2.
--fqU      | Fastq file with unpaired reads.
--ref_fa   | Reference fasta file.
--outdir   | Output directory (default: False).

*Options:*

Option      | Description
------------|-------------
--max_step  | Maximum number of refinement steps (default: 1).
--subsample | Use a subsample of reads for refinement.
--seed      | Seed for random number generator (ignored if not subsampling).
--sample_id | Sample ID. Used as read group ID in BAM (default: sampleXX).

*Settings:*

Option     | Description
-----------|-------------
--ncpu     | Number of CPUs to use (default: 1).
--xmx      | Maximum heap size for Java VM, in GB (default: 32).
--keep_tmp | Additional options (default: False).
--quiet    | Do not write output to console (silence stdout and stderr) (default: False).
--logfile  | Append console output to this file.
--debug    | Print commands but do not run (default: False).

#### Example usage

### *finalize_assembly*
Finalize consensus, map reads to consensus, and call variants. Input is reads in FASTQ format and reference sequence in FASTA format. Output is finalized reference sequence, alignment, and variants (in FASTA, BAM, and VCF formats, respectively).

#### Usage and option summary

**Usage:**

`haphpipe finalize_assembly [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> --ref_fa <FASTA> [--outdir]`

**(or):**

`hp_finalize_assembly [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --fqU <FASTQ> --ref_fa <FASTA> [--outdir]`

*Output files:* <br> final.fna <br> final.ban <br> final.vcf.gz

*Input/Output Arguments:* 

Option     | Description
-----------|-------------
--fq1      | Fastq file with read 1.
--fq2      | Fastq file with read 2.
--fqU      | Fastq file with unpaired reads.
--ref_fa   | Consensus fasta file.
--outdir   | Output directory (default: current directory).

*Options:*

Option       | Description
-------------|-------------
--bt2_preset | {very-fast,fast,sensitive,very-sensitive,very-fast-local,fast-local,sensitive-local,very-sensitive-local} Bowtie2 preset to use (default: very-sensitive).
--sample_id  | Sample ID (default: sampleXX).

*Settings:*

Option     | Description
-----------|-------------
--ncpu     | Number of CPU to use (default: 1).
--keep_tmp | Keep temporary directory (default: False).
--quiet    | Do not write output to console (silence stdout and stderr) (default: False).
--logfile  | Append console output to this file.
--debug    | Print commands but do not run (default: False).


#### Example usage


##	hp_haplotype
Haplotype assembly stages. HAPHPIPE implements PredictHaplo ([paper](https://www.ncbi.nlm.nih.gov/pubmed/26355517)), although other haplotype reconstruction programs can be utilized outside of HAPHPIPE using the final output of HAPHPIPE, typically with the final consensus sequence (FASTA) file, reads (raw, trimmed, and/or corrected), and/or final alignment (BAM) file as input.

### *predict_haplo*
Haplotype identification with PredictHaplo. Input is reads in FASTQ format and and reference sequence in FASTA format. Output is the longest global haplotype file and corresponding HTML file. _Note: PredictHaplo must be installed separately before running this stage._ 

#### Usage and option summary

**Usage:**

`haphpipe predict_haplo [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --ref_fa <FASTA> --interval_txt [TXT] [--outdir]`

**(or):**

`hp_ predict_haplo [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --ref_fa <FASTA> --interval_txt [TXT] [--outdir]`

*Output files:* <br> *best.fa*

*Input/Output Arguments:* 

Option         | Description
---------------|-------------
--fq1          | Fastq file with read 1.
--fq2          | Fastq file with read 2.
--ref_fa       | Reference sequence used to align reads (Fasta).
--interval_txt | File with intervals to perform haplotype reconstruction.
--outdir       | Output directory (default: current directory).

*Options:*

Option           | Description
-----------------|-------------
--min_readlength | Minimum read length passed to PredictHaplo (default: 36).

*Settings:*

Option      | Description
------------|-------------
--keep_tmp  | Keep temporary directory (default: False).
--quiet     | Do not write output to console (silence stdout and stderr) (default: False).
--logfile   | Append console output to this file.
--debug     | Print commands but do not run (default: False).


#### Example usage

### *ph_parser*
Returns PredictHaplo output as a correctly formatted FASTA file. Input is the output file from predict_haplo (longest global .fas file). Output is a correctly formatted FASTA file.

#### Usage and option summary

**Usage:**

`haphpipe ph_parser [OPTIONS] [SETTINGS] --haplotypes_fa <best.fas> [--outdir]`

**(or):**

`hp_ph_parser [OPTIONS] [SETTINGS] --haplotypes_fa <best.fas> [--outdir]`


*Output files:* <br> *ph_summary.txt* <br> *ph_haplotypes.fna*

*Input/Output Arguments:* 

Option          | Description
----------------|-------------
--haplotypes_fa | Haplotype file created by PredictHaplo.
--outdir        | Output directory (default: current directory).

*Options:*

Option      | Description
------------|-------------
--prefix    | Prefix to add to sequence names.
--keep_gaps | Do not remove gaps from alignment (default: False).

*Settings:*

Option    | Description
----------|-------------
--quiet   | Do not write output to console (silence stdout and stderr) (default: False).
--logfile |Append console output to this file.



#### Example usage
By default, PredictHaplo outputs their own unique version of a fasta file. It includes the frequency, some information regarding their unqiue overlapping scores, and their unique confidence scores. This file is always named `PH#.best_#_#.fas`, where the first number is the reconstructed haplotype number and the next numbers are the start and end of the longest haplotype reconstructed by PredictHaplo.

```
$ cat PH01.best_1_1208.fas
>reconstructed_0
;Freq:0.190638
;Overlap quality scores (5,10):1,1
;Confidence scores
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~n7~~~y~~~~~t~~~i~jjk~~zz
;~~~~{~~~~~N{~{~Sx~~~~~~z~K~F~~~~~~~y~~~~~~~|~~F~wx|~~{~~|~~s<|]~~~~~~
;~m~kj|{|v~{|_`~~~z~~~~~~~jy{y~~~~~a~__~~|~~~~~{wXZ~}~~~qm~xV~~~~~~~}~
;Q~}~y||~~~}}~~~z~~~~~~{}A~}b|~~u~~|}}|~}}z~}bx~~n||~~||{~}~~d}~bz~~~}
;|~}}~~~}~~~{|}}g{~~~}~r~}}~~~~u~~{kx{~}~}~|~}~}{}~}~}||~~~~[~}~}}~~~~
;~~}~~~U||U}~|}}~~}}~~~~}u~b|}~w~~~~~{}|wv~}Dxzp{|}~@~~P~}~}~V~~z}~|ry
;q~}|}~}~~}t~o~}~f}~~}{~~~}~{~~~~~~~}~~}~~|~~~M~~}~}x~}c~}^v~~~yzA~}~}
;y}~z}~~~~~~~~~{z~~}~~~}~{}~~~~~}~~~~~~~|~~v~~}~~|y~]|{~||~~~~~~~~||~~
;Y||~~|Q~|~~~|~~~~~|~~~z~~z~{{{y~~~~~~~~~~w{{~wz|~Z|~z|~~}p|~~|}}~~x}}
;z}~}}|a|}}}}{}|~~}}~}}~}~{~|~}}~}{}|}|}}}~|}}}}{}}}|}}|}|}}}|}}}}~}}}
;||}}}}{}}~}~}}}}}}}}}}}~}}}}}}}}}}}}}}}}}y}~}}}}}}}}}}~|}}}}|}}}}|}}}
;}}}}}}}}|}}}}}}}}~}}}}}}|}}}}}}}}|}}}}}|}}}}|~|}}}}{}}}}}}}}}}}}}}}}}
;}}=}}}{}}}}}}}}}}}}}}}}}}}}}}}}~||}}|~}}}}}{}}}}}|}}|}}}}}}}}}}}}}|}}
;|}~}}}}|}}}}}}}}|}|~J}}}}}}~}}}}}}}}}}}}}}}}}}}}|}}`~}}}}|}}}}}}}}}}~
;|}}}}}}}}}}}~}}}|}}}}}}}}}}|}}}}}}}}}}}t}}{}`}}}}}}|}~}~}|~}}}|k}}}}}
;|}|~}}|}}~}}~}|}}z}}}}}}}}}}}}}}~}}~}}}~~}|}}}~}}}}}}}}~}~}}||~~}~z}}
;~~}~}||~|~}|{}||~|z~~||}~|~}}|}~~}|}}~}}|}z~~~~{~}{}~y~~~~~{|}}~~|y~~
;~|~~||~~~~~~~|n~~{~~~~~~~~~~~~~~~
;EndOfComments
ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAATATCTGGCCTTCCCGCAAGGGAAG
GCCAGGGAACTTTCCTCAGAGCAGGCTAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGCTTGGGGAA
GCAACAACAAAGCCCTCTCAGAAACAGGAGCCAATAGACAAGGAAATGTATCCTTTAACTTCCCTCAGAT
CACTCTTTGGCAACGACCCCTTGTCACAATAAAGATAGGGGAGCAACTAAAGGAAGCCCTGTTAGATACA
GGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAGAATGATAGGGGGAA
TCGGGGGTTTTATCAAAGTAAGACAGTATGATCAGATAGCCATAGACATCTGTGGCCATAAAGCTATAGG
TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGCTTGGTTGCACT
TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTGAAGCCAGGAATGGATGGCCCAAAGG
TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAAGA
AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGGAAAAAAGAC
AGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAGGACTTCTGGGAAGTCC
AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGCGGGTGATGC
ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTATACTGCATTTACCATACCTAGTGTAAACAAT
GAGACACTAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
AAGCTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACAT
GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
GTTATGAACTCCATCCTT
>reconstructed_1
;Freq:0.294104
;Overlap quality scores (5,10):1,1
;Confidence scores
;~~~~~~~~z~~~~~~~~~|~y~~~~u}|~~~~~}}~~~~~~~ya~|~~}}~~~y~~~j~YXH~~s~
;~~~~z~~|~~b|}^}xZ}}~}z~t~r~{}~}x}}yb}~~~}~}u~~}~gu|}}{~|}}~ozsw}|}}h~
;|l}v|^][t}zz}vw}}s}}}}~v}|~|~~}}}~~}}}}}~}}}}}}zo\}}}~}vn}iv}}}}}}}}}
;v}}~|u}}}}}}}}}|}}}}}}|}^}}[}~}r}}}{|}}|t|}}{z}}pz}}{}}}}}}}v~}y|~}}}
;}}}}~}}~}}}tz}}`}}}~}}q}}}}}~}y}}|r||}}}}~}}~}}{}}}}}}}}}}}y~}}}}t}}}
;}}}}}}{|}|}}{r}}}}tv~}}}t}m|~}h}}}}}v}}qt}|y|{|u~}}|}~_}~|~~|}~|}~}|e
;Y}~}}~}}~}v}t}}}j}}~}}}~}}~|}}}}|}}|}}}~}}}}};}}~}}x}}|}}yu}}}wy`}|~~
;{}}x}}}}}~}}}~{v}}}}}{}~}}~}~}}|}}}}}}}}}}u}}y~~{g}L}}}}}|~}{~|}|}|}}
;l}}}}{_}}}}~z~}|}~~|~~x|~}~z~~|~}~~~}~~~}y}~~|k}}Pf}w~~~~T}}}F~~}~z}{
;{~|}|zZ~s~|}r~}}~~|}}}}}~|}}}}n}}}~}}}}~}}}}}}}}}}}|~}}}}}}~~}}~~~~~}
;N}}}}}|~}}~}}}}{}}}~}}|~~|}|}}}~|}~}}}}}}D}~}}}}}}}}}}}}}}}}|}}z}|}}}
;}k}}}}}}}}}}}}}}}}}}}|}}}}}|}}}}}|}}}}}|}}}}|~{~}}}G}~}}}}}~}}}}}}}}}
;}}h}}~|}|}}}}}}}}}}}}}}}}}}}}}}}}{}}|}}}}|~|}}}}}}}M}}}}}}}~}}}}}}}}}
;}}~|}}}}}~~~}}}~}}|}g}}}}}}}}~|}}}}}}}}}}}}}}}}}|}~t~}|}}}}}}}||~}}}}
;|}}}}~}}}}}}~}}}gi}}~}}}}}}}}}}}}|}}~}}g}~}}v}}}}}|{}}}}}}~}}~|e}}}}}
;}}|}}~}~}~}}~}|~}|}~|||}}}}|~}}~}}}}}}}}|}}~}}~}|~}|y~}~}}}}}|}}~~}~}
;}}~}~}|~}}}~}}}}~}{}}~}}~}}}||~~}||~}x|}{~|}|~~|~}}}|{~}}~}{}}||~}}}~
;|~}~}}}~|}}~~|,~~~~~~~~~~~~~~m~~~
;EndOfComments
ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCCGCGGCGGAAG
GCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAG
GAGACAACAACTCCCTCTCAGAAGCAGGAGCCGAAGGACAAGGAAACATATCCTTTAGCTTCCCTCAAAT
CACTCTTTGGCAACGACCCCTTGTTACAATAAAGATAGGGGGGCAGCTAAAGGAAGCTCTATTAGATACA
GGAGCAGATGACACAGTATTAGAAGAAATGGATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAA
TTGGAGGTTTTATCAAAGTAAAACAGTATGATCAGATACCCATAGAAATTTGTGGACATAAAGTGATAGG
TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGCTTGGTTGCACT
TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTGAAGCCAGGAATGGATGGCCCAAAGG
TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAATAGAAATCTGTGCAGAAATGGAAAAGGA
AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGAAAAAAAGAC
AGTACTAAATGGAGAAAATTAGTAGATTTCAAAGAACTTAATAAGAGAACTCAGGACTTCTGGGAAGTCC
AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAGTCAGTAACAGTACTGGATGTGGGTGATGC
ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTACACTGCATTTACCATACCTAGTGTAAACAAT
GAGACACCAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
AATGTAGCATGACAAAAATCTTAGAACCTTTTAGAAAACAAAATCCAGATATAGTTATCTATCAATACAT
GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
GTTATGAACTCCATCCTT
>reconstructed_2
;Freq:0.515259
;Overlap quality scores (5,10):4,4
;Confidence scores
;~~~~~~~w~~zz~~~~z~~w~|wy~|~~{|~~|~~~~~~~~}~{~hD~~~}|~|}}t~~_;~tue}}S}
;~}~~p}}l}~P}}ZlvL~z}}h~g~A~x}}}s}}{S~}~}}}}l}~|}SOq}}q}w}}~i\pz~w~~z~
;}7~|vzUu4|smtgm}~W}}~~}r}n`f`}{PZLSI:Vs_V_}}}}_<db}}}}}CQ}YJ}t}}}~}}}
;T}}~|}}}}}}}~}}k}}}}}}q}L}}]}}}L}}}tt|~}{q}}{M}lAIu}}}s{}}}~y~}st}|}}
;z}}}}|}}}}}{q|}Z}}}}}}P~}|}|}}E}}cLbi}vu}}|}}|}g}}}}}}|}~}~{}~}|}{}}}
;}~}|}}}{{|~}pp}~~~{|}}~}[yEx}}\}}}{}N}uoK}}pauSy}}}v~~J}}q}}v~}c}~}\i
;]|}}~}|}~}r}Y}}}B}}~||~}}}}{}}}|}|~|}~~~~}}~{w}{}~|R|}s|}k{~}}}ra}q}}
;o}}d}~}}}}}~}}|X~}}~}e}}}}}}}|}y|}}}}~}|~}c~~b~~nL~a}}}}}~}}m}}~{}}}}
;g}}}}uLm~l}}`~|n}}~}}|d~~y~u}|tz~~~~|}}~~t}}}tx}|Tu~v}}~}^|~~w~|~}m{}
;u}}}|w;}{}{}r}||}}z|}|~|}z}||}y}}}~}}}|}}~}}}|}y}}}}}}}~}||}}}}}}~}}|
;v{}|}}{}}}|~|}}|}}}}}}}}}}|}}}}}|}~}}}}}}p}}||}}}}}~}}}{}}}}}}}z}}}}}
;}|}}~}}}|}}}}}}}}}}}}}}}}}}}}}}|}||}}|}|}}}}|}|}}}}u}}}}}}}}}}}}}}}}}
;}}g}}}{}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}y}}}}}|}m}}}}}}}}}}}}}}}}}
;}}}}}}}}}}}}}}}}}}|}i}}}}}}}}}}|}}}}}}}}}}}}}}}}}}}H}}}}}}}}}|}{}}}}}
;|}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}y}}|}H}}}}}}}}}}}}}}}}}}u}}}}}
;}}}}}}}}}}}}}}}}}z}}}}}}}}}}}}}}}}|}}}|}}}}}}}}}|}}}}}}}}}}}}|}}}}|}}
;}}}}}}}}}}}}||}|}}|}}}}}}|}}|}}}~}|}|}}}|}{}}}}|}}}}|{|}}}}x||}}}|{}}
;}}}}}}}~}}~}~|n}~}}}}}~}}~~|~~~y~,
;EndOfComments
ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAAATCTGGCCTTCCCACAAGGGAAG
GCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAA
GAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGAAATGTATCCTTTAGCTTCCCTCAAAT
CACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACA
GGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAA
TTGGAGGTTTTATCAAAGTAAAACAGTATGATCAGATACCCATAGAAATCTGTGGACATAAAGCTATAGG
TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACT
TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAG
TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGA
AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGGAAAAAAGAC
AGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCCGGGAAGTCC
AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGC
ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTATACTGCATTTACCATACCTAGTGTAAACAAT
GAGACACCAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
AAGCTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACAT
GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
GTTATGAACTCCATCCCC
```

*ph_parser* takes this output and creates a proper fasta file with each resontructed haplotype and a text file that has the hpalotype diversity estimate.

```
$ cat ph_summary.txt
PH_num_hap 3
PH_hap_diversity 0.611668153059
PH_seq_len 1208


$ cat ph_haplotypes.fna
>reconstructed_0 Freq=0.190638
ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAATATCTGGCCTTCCCGCAAGGGAAG
GCCAGGGAACTTTCCTCAGAGCAGGCTAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGCTTGGGGAA
GCAACAACAAAGCCCTCTCAGAAACAGGAGCCAATAGACAAGGAAATGTATCCTTTAACTTCCCTCAGAT
CACTCTTTGGCAACGACCCCTTGTCACAATAAAGATAGGGGAGCAACTAAAGGAAGCCCTGTTAGATACA
GGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAGAATGATAGGGGGAA
TCGGGGGTTTTATCAAAGTAAGACAGTATGATCAGATAGCCATAGACATCTGTGGCCATAAAGCTATAGG
TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGCTTGGTTGCACT
TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTGAAGCCAGGAATGGATGGCCCAAAGG
TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAAGA
AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGGAAAAAAGAC
AGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAGGACTTCTGGGAAGTCC
AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGCGGGTGATGC
ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTATACTGCATTTACCATACCTAGTGTAAACAAT
GAGACACTAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
AAGCTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACAT
GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
GTTATGAACTCCATCCTT
>reconstructed_1 Freq=0.294104
ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCCGCGGCGGAAG
GCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAG
GAGACAACAACTCCCTCTCAGAAGCAGGAGCCGAAGGACAAGGAAACATATCCTTTAGCTTCCCTCAAAT
CACTCTTTGGCAACGACCCCTTGTTACAATAAAGATAGGGGGGCAGCTAAAGGAAGCTCTATTAGATACA
GGAGCAGATGACACAGTATTAGAAGAAATGGATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAA
TTGGAGGTTTTATCAAAGTAAAACAGTATGATCAGATACCCATAGAAATTTGTGGACATAAAGTGATAGG
TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGCTTGGTTGCACT
TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTGAAGCCAGGAATGGATGGCCCAAAGG
TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAATAGAAATCTGTGCAGAAATGGAAAAGGA
AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGAAAAAAAGAC
AGTACTAAATGGAGAAAATTAGTAGATTTCAAAGAACTTAATAAGAGAACTCAGGACTTCTGGGAAGTCC
AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAGTCAGTAACAGTACTGGATGTGGGTGATGC
ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTACACTGCATTTACCATACCTAGTGTAAACAAT
GAGACACCAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
AATGTAGCATGACAAAAATCTTAGAACCTTTTAGAAAACAAAATCCAGATATAGTTATCTATCAATACAT
GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
GTTATGAACTCCATCCTT
>reconstructed_2 Freq=0.515259
ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAAATCTGGCCTTCCCACAAGGGAAG
GCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAA
GAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGAAATGTATCCTTTAGCTTCCCTCAAAT
CACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACA
GGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAA
TTGGAGGTTTTATCAAAGTAAAACAGTATGATCAGATACCCATAGAAATCTGTGGACATAAAGCTATAGG
TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACT
TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAG
TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGA
AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGGAAAAAAGAC
AGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCCGGGAAGTCC
AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGC
ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTATACTGCATTTACCATACCTAGTGTAAACAAT
GAGACACCAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
AAGCTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACAT
GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
GTTATGAACTCCATCCCC
```


			
			
##	hp_annotate
### *pairwise_align*
Apply correct coordinate system to final sequence(s) to facilitate downstream analyses. Input is the final sequence file in FASTA format, a reference sequence in FASTA format, and a reference GFT file. Output is a JSON file to be used in _extract_pairwise_.

#### Usage and option summary

**Usage:**

`haphpipe pairwise_align [SETTINGS] --amplicons_fa <FASTA> --ref_fa <FASTA> --ref_gtf <GTF> [--outdir]`

**(or):**

`hp_pairwise_align [SETTINGS] --amplicons_fa <FASTA> --ref_fa <FASTA> --ref_gtf <GTF> [--outdir]`

*Output files:* <br> 
pairwise_aligned.json

*Input/Output Arguments:* 

Option          | Description
----------------|-------------
--amplicons_fa  | Fasta file with assembled amplicons.
--ref_fa        | Reference fasta file.
--ref_gtf       | GTF format file containing amplicon regions. Primary and alternate coding regions should be provided in the attribute field (for amino acid alignment).
--outdir        | Output directory (default: False).

*Settings:*

Option      | Description
------------|-------------
--keep_tmp  | Keep temporary directory (default: False).
--quiet     | Do not write output to console (silence stdout and stderr) (default: False).
--logfile   | Append console output to this file.
--debug     | Print commands but do not run (default: False).


#### Example usage

### *extract_pairwise*
Extract sequence regions from the pairwise alignment produced in _pairwise_align_. Input is the JSON file from _pairwise_align_. Output is either an unaligned nucleotide FASTA file, an aligned nucleotide FASTA file, an amino acid FASTA file, an amplicon GTF file, or a tab-separated values (TSV) file (default: nucleotide FASTA with regions of interest from GTF file used in _pairwise_align_). 

#### Usage and option summary

**Usage:**

`haphpipe extract_pairwise [OPTIONS] [SETTINGS] --align_json <JSON> [--outdir]`

**(or):**

`hp_extract_pairwise [OPTIONS] [SETTINGS] --align_json <JSON> [--outdir]`

*Output files:* <br> 
stdout.fasta

*Input/Output Arguments:* 

Option        | Description
--------------|-------------
--align_json  | JSON file describing alignment (output of _pairwise_align_ stage).
--outfile     | Output file (default: stdout).

*Options:*

Option    | Description
----------|-------------
--outfmt  | Format for output: nuc_fa, aln_fa, amp_gtf, ost, or prot_fa (default: nuc_fa).
--refreg  | Reference region. String format is ref:start-stop. For example, the region string to extract pol when aligned to HXB2 is HIV_B.K03455.HXB2:2085-5096.

*Settings:*

Option    | Description
----------|-------------
--debug   | Print commands but do not run (default: False).

#### Example usage

### * annotate_ from_ref*
Annotate consensus sequence from reference annotation. Input is JSON file from _pairwise_align_ and reference GTF file. 

#### Usage and option summary

**Usage:**

`haphpipe annotate_from_ref [OPTIONS] [SETTINGS] --haplotypes_fa <best.fas> [--outdir]`

**(or):**

`hp_annotate_from_ref [OPTIONS] [SETTINGS] --haplotypes_fa <best.fas> [--outdir]`

#### Example usage


</br>

---
---

# Example usage
## `hp_assemble_01`

## `hp_assemble_02`

## HIV
Data from [NCBI bioproject: PRJNA486832](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA486832) with [associated article](https://www.biorxiv.org/content/biorxiv/early/2018/09/12/414995.full.pdf). Brief description of the paper and data: **insert from students**

### Quick start
#### Miniconda and Bioconda
Step 1. Installing [conda](https://bioconda.github.io/user/install.html#set-up-channels).

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh`

Step 2. Setting up [conda channels](https://bioconda.github.io/user/install.html#set-up-channels).

```
conda config --add channels R
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### Downloading data
This is showing how to load miniconda on GW's HPC, Colonial One, which is how undergrads in CBI work on data. 

```
module use /groups/cbi/shared/modulefiles
module load miniconda3/4.3.27.1
```

Creating a new conda environment to pull data from SRA.

```
conda create -n fastqdump sra-tools
conda activate fastqdump
```

Downloading SRA data.
This pulls the first 20 accessions from the SRA list, and only pulls 20,000 reads from SRA data. We pulled the first 20 in order, because thier paper says that between 22-24 samples were put on a single MiSeq run.

This code says for each accession in the accession list, sort them in order, and then take the first 20 accession numbers. For each of these 20 accession numbers, pull 20,000 reads from the SRA data and put the output read files into the accession folder. Then deactive the conda environment used to pull these SRA data.

```
for acc in $(cat SraAccList.txt | sort | head -n 20); do
    mkdir -p $acc
    fastq-dump --outdir $acc --split-files -N 10001 -X 30000 $acc
done

# deactivate the conda environment now that you are done downloading data
conda deactivate
```

#### Installing HAPHPIPE
**insert from students**

#### Activating HAPHPIPE
**insert from students**


### HAPHPIPE suite
**insert from students**

#### hp_reads
        sample_reads
        trim_reads
        join_reads
        ec_reads
#### hp_assemble
        assemble_denovo
        assemble_amplicons
        assemble_scaffold
        align_reads
        call_variants
        vcf_to_consensus
        refine_assembly
        finalize_assembly
#### hp_haplotype
        predict_haplo
        ph_parser
#### hp_annotate
        pairwise_align
        extract_pairwise
        annotate_from_ref
### Example usage

#### Your virus with both pipelines. Document all code and explanation.

##### Starting directory structure
**insert from students**

##### Ending directory structure
**insert from students**

### Helpful sources
	bash scripting
	conda
### FAQ


## HCV
Data from [NCBI bioproject: PRJNA237978](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=651945) with [associated article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4067308). Brief description of the paper and data: **insert from students**

### Quick start
#### Miniconda and Bioconda
Step 1. Installing [conda](https://bioconda.github.io/user/install.html#set-up-channels).

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh`

Step 2. Setting up [conda channels](https://bioconda.github.io/user/install.html#set-up-channels).

```
conda config --add channels R
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### Downloading data
This is showing how to load miniconda on GW's HPC, Colonial One, which is how undergrads in CBI work on data. 

```
module use /groups/cbi/shared/modulefiles
module load miniconda3/4.3.27.1
```

Creating a new conda environment to pull data from SRA.

```
conda create -n fastqdump sra-tools
conda activate fastqdump
```

Downloading SRA data.
This pulls the all the accessions (23) from the SRA list, and only pulls 20,000 reads from SRA data. 

This code says for each accession in the accession list, pull 20,000 reads from the SRA data and put the output read files into the accession folder. Then deactive the conda environment used to pull these SRA data.

```
for acc in $(cat SraAccList.txt); do
    mkdir -p $acc
    fastq-dump --outdir $acc --split-files -N 10001 -X 30000 $acc
done

# deactivate the conda environment now that you are done downloading data
conda deactivate
```

#### Installing HAPHPIPE
**insert from students**

#### Activating HAPHPIPE
**insert from students**


### HAPHPIPE suite
**insert from students**

#### hp_reads
        sample_reads
        trim_reads
        join_reads
        ec_reads
#### hp_assemble
        assemble_denovo
        assemble_amplicons
        assemble_scaffold
        align_reads
        call_variants
        vcf_to_consensus
        refine_assembly
        finalize_assembly
#### hp_haplotype
        predict_haplo
        ph_parser
#### hp_annotate
        pairwise_align
        extract_pairwise
        annotate_from_ref
### Example usage

#### Your virus with both pipelines. Document all code and explanation.

##### Starting directory structure
**insert from students**

##### Ending directory structure
**insert from students**

### Helpful sources
	bash scripting
	conda
	**insert from students**
### FAQ
**insert from students**


## Influenza
Data from [NCBI bioproject: PRJNA506454](https://www.ncbi.nlm.nih.gov/bioproject/506454) with [associated article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731279). Brief description of the paper and data: **insert from students**

### Quick start
#### Miniconda and Bioconda
Step 1. Installing [conda](https://bioconda.github.io/user/install.html#set-up-channels).

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh`

Step 2. Setting up [conda channels](https://bioconda.github.io/user/install.html#set-up-channels).

```
conda config --add channels R
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### Downloading data
This is showing how to load miniconda on GW's HPC, Colonial One, which is how undergrads in CBI work on data. 

```
module use /groups/cbi/shared/modulefiles
module load miniconda3/4.3.27.1
```

Creating a new conda environment to pull data from SRA.

```
conda create -n fastqdump sra-tools
conda activate fastqdump
```

Downloading SRA data.
This pulls the 20 random accessions from the SRA list, and only pulls 20,000 reads from SRA data. We pulled randomly, because the study had two sets of influenza types sequenced, and we wanted to pull from both, but randomly

This code says for each accession in the accession list, sort them randomly (the `-R` option applied in `sort`, and then take the first 20 accession numbers from that randomly sorted list. For each of these 20 accession numbers, pull 20,000 reads from the SRA data and put the output read files into the accession folder. Then deactive the conda environment used to pull these SRA data.

```
for acc in $(cat SraAccList.txt | sort -R | head -n 20); do
    mkdir -p $acc
    fastq-dump --outdir $acc --split-files -N 10001 -X 30000 $acc
done

# deactivate the conda environment now that you are done downloading data
conda deactivate
```

#### Installing HAPHPIPE
**insert from students**

#### Activating HAPHPIPE
**insert from students**


### HAPHPIPE suite
**insert from students**

#### hp_reads
        sample_reads
        trim_reads
        join_reads
        ec_reads
#### hp_assemble
        assemble_denovo
        assemble_amplicons
        assemble_scaffold
        align_reads
        call_variants
        vcf_to_consensus
        refine_assembly
        finalize_assembly
#### hp_haplotype
        predict_haplo
        ph_parser
#### hp_annotate
        pairwise_align
        extract_pairwise
        annotate_from_ref
### Example usage

#### Your virus with both pipelines. Document all code and explanation.

##### Starting directory structure
**insert from students**

##### Ending directory structure
**insert from students**

### Helpful sources
	bash scripting
	conda
	**insert from students**
### FAQ
**insert from students**


</br>

---
---

# Advanced Usage
	individual modules + making your own pipelines w/ bash in advanced usage

</br>

---
---

# File extension table

| Module                 | Input                                | Input Format                                          | Output                                                                  | Output File Names                                                                                       |
|------------------------|--------------------------------------|-------------------------------------------------------|-------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
| **sample_reads**       | FASTQ file(s)                        | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_            | Subampled FASTQ file(s)                                                 | *sample_1.fastq* <br> *sample_2.fastq*                                                                  |
| **trim_reads**         | FASTQ file(s)                        | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_            | Trimmed FASTQ file(s) and <br> summary from Trimmomatic text file       | *trimmed_1.fastq* <br> *trimmed_2.fastq* <br> *trimmed_U.fastq* <br> *trimmomatic_summary.out*          |
| **join_reads**         | FASTQ file(s)                        | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_            | FASTQ files for joined reads and unjoined reads                         | _joined.fastq_ <br> *notjoined_1.fastq* <br> *notjoined_2.fastq*                                        |
| **ec_reads**           | FASTQ file(s)                        | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_            | Error corrected FASTQ file(s)                                           | *corrected_1.fastq* <br> *corrected_2.fastq* <br> *corrected_U.fastq*                                   |
| **assemble_denovo**    | FASTQ file(s)                        | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_            | De novo assembled contigs FASTA file and <br> de novo summary text file | *denovo_contigs.fna* <br> *denovo_summary.txt*                                                          |
| **assemble_amplicons** | FASTA files and <br> GTF file        | _.fna_ or <br> _.fasta_ or <br> _.fa_ and <br> _.gtf_ | FASTA file with assembled amplicons                                     | *amplicon_assembly.fna*                                                                                 |
| **assemble_scaffold**  | FASTA files                          | _.fna_ or <br> _.fasta_ or <br> _.fa_                 | FASTA files with scaffolded and aligned sequences and <br> FASTA file with assembled amplicons | *scaffold_aligned.fa* <br> *scaffold_assembly.fa* <br> *scaffold_imputed.fa* <br> *scaffold_padded.out* |
| **align_reads**        | FASTQ files and <br> FASTA file      | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_ and <br> _.fna_ or <br> _.fasta_ or <br> _.fa_ | Aligned BAM file and <br> Bowtie2 alignment output summary text file | _aligned.bam_ <br> _aligned.bt2.out_                                |
| **call_variants**      | BAM and <br> FASTA file              | _.bam_ and <br> _.fna_ or <br> _.fasta_ or <br> _.fa_ | VCF file with variants                                                  | _variants.vcf.gz_                                                                                       |
| **vcf_to_consensus**   | VCF file                             | _.vcf_                                                | FASTA file with consensus sequence                                      | _consensus.fna_                                                                                         |
| **refine_assembly**    | FASTQ files and <br> FASTA file      | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_ and <br> _.fna_ or <br> _.fasta_ or <br> _.fa_ | FASTA file with refined consensus sequence | _refined.fna_                                                       |
| **finalize_assembly**  | FASTQ files and <br> FASTA file      | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_ and <br> _.fna_ or <br> _.fasta_ or <br> _.fa_ | FASTA file with final refined consensus sequence and <br> final BAM file with reads aligned to final FASTA file and <br> VCF file with variants relative to *final.fna* | _final.fna_ and <br> _final.bam_ and <br> _final.vcf.gz_            |
| **predict_haplo**      | FASTQ files and <br> FASTA file      | _.fastq.gz_ or <br> _.fastq_ or <br> _.fq_ and <br> _.fna_ or <br> _.fasta_ or <br> _.fa_ | PredictHaplo's fasta-like output file| _best.fa_                                                           |
| **ph_parser**          | PredictHaplo's FAS output file       | _.fas_                                                | Summary text file with haplotype diversity statistics and <br> FASTA file with haplotype sequences | *ph_summary.txt* and <br> *ph_haplotypes.fna*                                                           |
| **pairwise_align**     | FASTA files and <br> GTF file        | _.fna_ or <br> _.fasta_ or <br> _.fa_ and <br> _.gtf_ | JSON file                                                               | *pairwise_aligned.json*                                                                                 |
| **extract_pairwise**   | JSON file from pairwise_align output | *.json*                                               | FASTA output with region extracted to standard out                      | *stdout.fasta*                                                                                          |
| **annotate_from_ref**  | JSON file from pairwise_align output and GTF file | _.json_ and <br> _.gtf_              | | |

</br>

---
---

# Helpful resources
- Conda 
	- [Bioconda](https://bioconda.github.io)
	- [Bioconda channels](https://bioconda.github.io/user/install.html#set-up-channels)
- Bash 
	- [Command line tutorial](https://www.codecademy.com/learn/learn-the-command-line)
	- [List of bash commands](https://www.codecademy.com/articles/command-line-commands)
	- [Bash scripting tutorial](https://www.codecademy.com/learn/learn-the-command-line/modules/bash-scripting)
	- Additional scripting help:  
		- [linuxconfig.org](https://linuxconfig.org/bash-scripting-tutorial-for-beginners)
		- [flaviocopesc.com](https://flaviocopes.com/bash-scripting/)
		- [ryanstutorials.net](https://ryanstutorials.net/bash-scripting-tutorial/)
- NGS overview:
	-  [NYU resource](https://learn.gencore.bio.nyu.edu)
- Accessing on a PC using a virtual machine
	- [Guide to options](https://blog.storagecraft.com/the-dead-simple-guide-to-installing-a-linux-virtual-machine-on-windows/)
	-  [VirtualBox](https://www.virtualbox.org)
</br>

---
---

# FAQ
- what is directory structure?


</br>

---
---

# Related software