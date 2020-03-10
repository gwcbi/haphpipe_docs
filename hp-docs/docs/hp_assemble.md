Stages in *hp_assemble* are designed to construct consensus sequence(s). Input reads (in FASTQ format) are assembled 
using either denovo assembly or reference-based alignment. Resulting consensus can be further refined. Use -h after any command for a list of options.

### *assemble_denovo*
Assemble reads via de novo assembly using SPAdes ([documentation](http://cab.spbu.ru/software/spades/)). Input is reads in FASTQ format. Output is contigs in FNA format.

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

_Example usage:_


### *assemble_amplicons*
Assemble contigs from de novo assembly using both a reference sequence and amplicon regions with MUMMER 3+ ([documentation](http://mummer.sourceforge.net/manual/)). Input is contigs and reference sequence in FASTA format and amplicon regions in GTF format.

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


_Example usage:_

### *assemble_scaffold*
Scaffold contigs against a reference sequence with MUMMER 3+ ([documentation](http://mummer.sourceforge.net/manual/)). Input is contigs in FASTA format and reference sequence in FASTA format. Output is scaffold assembly, alligned scaffold, imputed scaffold, and padded scaffold in FASTA format.



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

_Example usage:_

### *align_reads*
Map reads to reference sequence (instead of running de novo assembly) using Bowtie2 ([documentation](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)) and Picard ([documentation](https://broadinstitute.github.io/picard/)). Input is reads in FASTQ format and reference sequence in FASTA format. 

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

_Example usage:_

### *call_variants*
Variant calling from alignment using GATK ([documentation](https://software.broadinstitute.org/gatk/download/archive)). Input is alignment file in BAM format and reference sequence in FASTA format (either reference from reference-based assembly or consensus final sequence from de novo assembly). Output is a Variant Call File (VCF) format file. 



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

_Example usage:_

### *vcf_to_consensus*
Generate a consensus sequence from a VCF file. Input is a VCF file. Output is the consensus sequence in FASTA format. 

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


_Example usage:_

### *refine_assembly*
Map reads to a denovo assembly or reference alignment. Assembly or alignment is iteratively updated. Input is reads in FASTQ format and reference sequence (assembly or reference alignment) in FASTA format. Output is refined assembly in FASTA format.

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

_Example usage:_

### *finalize_assembly*
Finalize consensus, map reads to consensus, and call variants. Input is reads in FASTQ format and reference sequence in FASTA format. Output is finalized reference sequence, alignment, and variants (in FASTA, BAM, and VCF formats, respectively).

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


_Example usage:_

