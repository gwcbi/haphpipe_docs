## hp_reads
Stages to manipulate reads. Use -h after any command for a list of options. 

### sample_reads
Subsample reads using seqtk ([seqtk documentation](https://github.com/lh3/seqtk)).

*Output files:* <br> sample_1.fastq <br> sample_2.fastq

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         |Fastq file with read 1.
--fq2         |Fastq file with read 2.
--fqU         |Fastq file with unpaired reads.
--outdir      |Output directory (default: current directory).

*Sample Options:*

Option                   | Description
-------------------------|-------------
--nreads                 |Number of reads to sample. If greater than the number of reads in file, all reads will be sampled.
--frac                   |Fraction of reads to sample, between 0 and 1. Each read has [frac] probability of being sampled, so the number of sampled reads is not precisely [frac * number of reads].
--seed                   |Seed for random number generator.

*Settings:*

Option        | Description
------------- |-------------
--quiet       |Do not write output to console (silence stdout and stderr), default is False.
--logfile     |Append console output to this file.
--debug       |Print commands but do not run, default is False.

### trim_reads
Trim reads using Trimmomatic ([Trimmomatic documentation](https://github.com/timflutre/trimmomatic)).

*Output files:* <br> trimmed_1.fastq <br> trimmed_2.fastq <br> trimmed_U.fastq <br> trimmomatic_summary.out

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         |Fastq file with read 1.
--fq2         |Fastq file with read 2.
--fqU         |Fastq file with unpaired reads.
--outdir      |Output directory (default: current directory).

*Trimmomatic Options:*

Option                   | Description
-------------------------|-------------
--adapter_file           |Adapter file.
--trimmers               |Trim commands for trimmomatic (default: ['LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36']).
--encoding               |Quality score encoding.

*Settings:*

Option        | Description
------------- |-------------
--ncpu        |Number of CPU to use (default: 1).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

### join_reads
Join reads using FLASH ([FLASH documentation](https://github.com/dstreett/FLASH2)).

*Output files:* <br> joined.fastq <br> notjoined_1.fastq <br> notjoined_2.fastq

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         |Fastq file with read 1.
--fq2         |Fastq file with read 2.
--outdir      |Output directory (default: current directory).

*FLASH Settings:*

Option                   | Description
-------------------------|-------------
--min_overlap            |The minimum required overlap length between two reads to provide a confident overlap (default: 10).
--max_overlap            |Maximum overlap length expected in approximately 90% of read pairs, longer overlaps are penalized.
--allow_outies           |Also try combining read pairs in the "outie" orientation (default: False).
--encoding               |Quality score encoding.

*Settings:*

Option        | Description
------------- |-------------
--ncpu        |Number of CPU to use (default: 1).
--keep_tmp    |Keep temporary directory (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

### ec_reads
Error correct reads using spades ([Spades documentation](https://github.com/ablab/spades)).

*Output files:* <br> corrected_1.fastq <br> corrected_2.fastq <br> corrected_U.fastq

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         |Fastq file with read 1.
--fq2         |Fastq file with read 2.
--fqU         |Fastq file with unpaired reads.
--outdir      |Output directory (default: current directory).

*Settings:*

Option        | Description
------------- |-------------
--ncpu        |Number of CPU to use (default: 1).
--keep_tmp    |Keep temporary directory (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run, default is False.

## hp_assemble

Stages to assemble consensus sequence(s). Use -h after any command for a list of options.

### assemble_denovo
Assemble reads using denovo assembly.

*Output files:* <br> denovo_contigs.fna <br> denovo_summary.txt

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         |Fastq file with read 1.
--fq2         |Fastq file with read 2.
--fqU         |Fastq file with unpaired reads.
--outdir      |Output directory (default: current directory).

*Assembly Options:*

Option                   | Description
-------------------------|-------------
--no_error_correction| Do not perform error correction (default: False)
--subsample |Use a subsample of reads for assembly
--seed                   |Seed for random number generator (ignored if not subsampling)

*Settings:*

Option        | Description
------------- |-------------
--ncpu | Number of CPU to use (default: 1).
--keep_tmp | Keep temporary directory (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

### assemble_amplicons

Assemble contigs using reference and amplicon regions.

*Output files:* <br> amplicon_assembly.fna

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--contigs_fa         |Fasta file with assembled contigs.
--ref_fa         |Fasta file with reference genome to scaffold against.
--ref_gtf         |GTF format file containing amplicon regions.
--outdir      |Output directory (default: current directory).

*Scaffold Options:*

Option                   | Description
-------------------------|-------------
--sample_id| Sample ID (default: sampleXX).
--padding |Bases to include outside reference annotation (default: 50).

*Settings:*

Option        | Description
------------- |-------------
--keep_tmp | Keep temporary directory (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

### assemble_scaffold

Scaffold contigs using reference sequence.

*Output files:* <br> scaffold_aligned.fa <br> scaffold_assembly.fa <br> scaffold_imputed.fa <br> scaffold_padded.out

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--contigs_fa         |Fasta file with assembled contigs.
--ref_fa       |Fasta file with reference genome to scaffold against.
--outdir      |Output directory (default: current directory).

*Scaffold Options:*

Option                   | Description
-------------------------|-------------
--seqname           |Name to append to scaffold sequence (default: sample01).

*Settings:*

Option        | Description
------------- |-------------
--keep_tmp        |Additional options (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

### align_reads

Align reads to reference.

*Output files:* <br> aligned.bam <br> aligned.bt2.out

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         |Fastq file with read 1.
--fq2         |Fastq file with read 2.
--fqU         |Fastq file with unpaired reads.
--ref_fa.     |Reference fasta file.
--outdir      |Output directory (default: current directory).

*Alignment Options:*

Option                   | Description
-------------------------|-------------
--bt2_preset           |{very-fast, fast, sensitive,very-sensitive,very-fast-local,fast-local,sensitive-local,very-sensitive-local}
--sample_id |Sample ID. Used as read group ID in BAM (default: sampleXX).
--no_realign |Do not realign indels (default: False).
--remove_duplicates |Remove duplicates from final alignment. Otherwise duplicates are marked but not removed (default: False).
--encoding |{Phred+33,Phred+64} Quality score encoding.

*Settings:*

Option        | Description
------------- |-------------
--ncpu    |Number of CPUs to use (default: 1).
--xmx | Maximum heap size for Java VM, in GB (default: 32).
--keep_tmp        |Additional options (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

### call_variants

Call variants.

*Output files:* <br> variants.vcf.gz

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--aln_bam         |Alignment file.
--ref_fa       |Reference fasta file.
--outdir      |Output directory (default: False).

*Variant Calling Options:*

Option                   | Description
-------------------------|-------------
--emit_all           |Output calls for all site (default: False).
--min_base_qual | Minimum base quality required to consider a base for calling (default: 15).

*Settings:*

Option        | Description
------------- |-------------
--ncpu    |Number of CPUs to use (default: 1).
--xmx | Maximum heap size for Java VM, in GB (default: 32).
--keep_tmp        |Additional options (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

### vcf_to_consensus

Create consensus sequence from VCF.

*Output files:* <br> consensus.fna

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--vcf         |VCF file (created with all sites).
--outdir      |Output directory (default: False).
--sampidx | Index for sample if multi-sample VCF (default: 0).

*Variant Options:*

Option                   | Description
-------------------------|-------------
--min_DP           | Minimum depth to call site (default: 1).
--major |Allele fraction to make unambiguous call (default: 0.5).
--minor | Allele fraction to make ambiguous call (default: 0.2).

*Settings:*

Option        | Description
------------- |-------------
--keep_tmp        |Additional options (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.

### refine_assembly

Three step assembly refinement: align reads, call variants, and update reference.

*Output files:* <br> refined.fna

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         |Fastq file with read 1.
--fq2         |Fastq file with read 2.
--fqU         |Fastq file with unpaired reads.
--ref_fa.     |Reference fasta file.
--outdir      |Output directory (default: False).

*Refinement Options:*

Option                   | Description
-------------------------|-------------
--max_step           |Maximum number of refinement steps (default: 1).
--subsample |Use a subsample of reads for refinement.
--seed |Seed for random number generator (ignored if not subsampling).
--sample_id |Sample ID. Used as read group ID in BAM (default: sampleXX).


*Settings:*

Option        | Description
------------- |-------------
--ncpu    |Number of CPUs to use (default: 1).
--xmx | Maximum heap size for Java VM, in GB (default: 32).
--keep_tmp        |Additional options (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

### finalize_assembly
Finalize consensus sequence, align all reads to consensus, and call variants in dataset.

*Output files:* <br> final.fna <br> final.ban <br> final.vcf.gz

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--fq1         |Fastq file with read 1.
--fq2         |Fastq file with read 2.
--fqU         |Fastq file with unpaired reads.
--ref_fa |Consensus fasta file.
--outdir      |Output directory (default: current directory).

*Consensus Options:*

Option                   | Description
-------------------------|-------------
--bt2_preset|{very-fast,fast,sensitive,very-sensitive,very-fast-local,fast-local,sensitive-local,very-sensitive-local} Bowtie2 preset to use (default: very-sensitive).
--sample_id |Sample ID (default: sampleXX).

*Settings:*

Option        | Description
------------- |-------------
--ncpu | Number of CPU to use (default: 1).
--keep_tmp | Keep temporary directory (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).

## hp_annotate

Stages to annotate consensus sequence(s). Use -h after any command for a list of options.

### pairwise_align
Align amplicons to reference.

*Output files:* <br> 
pairwise_aligned.json

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--amplicons_fa         |Fasta file with assembled amplicons.
--ref_fa         |Reference fasta file.
--ref_gtf         |GTF format file containing amplicon regions. Primary and alternate coding regions should be provided in the attribute field (for amino acid alignment).
--outdir      |Output directory (default: False).

*Settings:*

Option        | Description
------------- |-------------
--keep_tmp | Keep temporary directory (default: False).
--quiet       |Do not write output to console (silence stdout and stderr) (default: False).
--logfile     |Append console output to this file.
--debug       |Print commands but do not run (default: False).


### extract_pairwise
Extract sequence regions from pairwise alignment.

*Output files:* <br> 
stdout.fasta

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--align_json         |JSON file describing alignment (output of _pairwise_align_ stage).
--outfile         |Output file (default: stdout).

*Options:*

Option        | Description
------------- |-------------
--outfmt       |Format for output: nuc_fa, aln_fa, amp_gtf, ost, or prot_fa (default: nuc_fa).
--refreg |Reference region. String format is ref:start-stop. For example, the region string to extract pol when aligned to HXB2 is HIV_B.K03455.HXB2:2085-5096.

*Settings:*

Option        | Description
------------- |-------------
--debug       |Print commands but do not run (default: False).


### hp_annotate_from_ref
Extract sequence regions from pairwise alignment.

*Output files:* <br> 

*Input/Output Arguments:* 

Option        | Description
------------- |-------------
--align_json         |JSON file describing alignment (output of _pairwise_align_ stage).
--ref_gtf         |GTF format file containing amplicon regions. 
--outfile         |Output file (default: stdout).

*Settings:*

Option        | Description
------------- |-------------
--debug       |Print commands but do not run (default: False).


