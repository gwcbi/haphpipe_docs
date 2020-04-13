The Reads stage involves cleaning up the raw read sequences, as well as other processing steps. Modules to manipulate reads. Use -h after any command for a list of options. 


### *sample_reads*
Subsample reads using seqtk ([documentation](https://github.com/lh3/seqtk)). Input is reads in FASTQ format. Output is sampled reads in FASTQ format. You do not have to have all read options (i.e., read1, read2 AND unpaired reads). You can have a combination of any of those.

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

_Example usage:_

This pulls 1000 reads from these paired end files with a starting seed of 1234.
```
haphpipe sample_reads --fq1 read_1.fastq --fq2 read_2.fastq --nreads 1000 --seed 1234
```

--
	
### *trim_reads*
Trim reads using Trimmomatic ([documentation](http://www.usadellab.org/cms/?page=trimmomatic)). Input is reads in FASTQ format. Output is trimmed reads in FASTQ format.

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

_Example usage:_

This trims paired end read files 1 and 2.

```
haphpipe trim_reads --fq1 read_1.fastq --fq2 read_2.fastq 
```

--

### *join_reads*
Join reads using FLASH ([paper](https://www.ncbi.nlm.nih.gov/pubmed/21903629)). Input is reads in FASTQ format. Output is joined reads in FASTQ format.

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

_Example usage:_

```
haphpipe join_reads --fq1 trimmed_1.fastq --fq2 trimmed_2.fastq
```

### *ec_reads*
Error correction using SPAdes ([documentation](http://cab.spbu.ru/software/spades/)). Input is reads in FASTQ format. Output is error-corrected reads in FASTQ format. Remember that HAPHPIPE is intended for Illumina reads, therefore the error correction is based on Illumina sequencing errors.

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

_Example usage:_

```
haphpipe ec_reads --fq1 trimmed_1.fastq --fq2 trimmed_2.fastq
```
