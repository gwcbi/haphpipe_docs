*hp_annotate* includes stages to annotate consensus sequence(s). Use -h after any command for a list of options.

### *pairwise_align*
Apply correct coordinate system to final sequence(s) to facilitate downstream analyses. Input is the final sequence file in FASTA format, a reference sequence in FASTA format, and a reference GFT file. Output is a JSON file to be used in _extract_pairwise_.

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


_Example usage:_
```
haphpipe pairwise_align --amplicons_fa final.fna --ref_fa HIV_B.K03455.HXB2.fasta --ref_gtf HIV_B.K03455.HXB2.gtf
```

### *extract_pairwise*
Extract sequence regions from the pairwise alignment produced in _pairwise_align_. Input is the JSON file from _pairwise_align_. Output is either an unaligned nucleotide FASTA file, an aligned nucleotide FASTA file, an amino acid FASTA file, an amplicon GTF file, or a tab-separated values (TSV) file (default: nucleotide FASTA with regions of interest from GTF file used in _pairwise_align_). 

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

_Example usage:_
```
haphpipe extract_pairwise --align_json pairwise_aligned.json --refreg HIV_B.K03455.HXB2:2085-5096
```

### *summary_stats*
Report summary statistics from an alignment and/or haplotype calling as TXT and TSV files. Input is a list of paths to directories (TXT format, one per line), each of which contain the following files: `final_bt2.out`, `trimmomatic_summary.out`, `final.bam`, `final.fna`, and `final.vcf.gz`.
If applicable, also input a list of directories containing PredictHaplo summary files (`ph_summary.txt`). If amplicons were used in assembly, use the `--amplicons` option to report statistics per amplicon.

**Usage:**

`haphpipe summary_stats [SETTINGS] --dir_list <TXT> [--ph_list <TXT>] [--amplicons] [--outdir] `

**(or):**

`hp_summary_stats [SETTINGS] --dir_list <TXT> [--ph_list <TXT>] [--amplicons] [--outdir] `

*Output files:* <br> 
summary_stats.txt, summary_stats.tsv, PH_summary_stats.tsv

*Input/Output Arguments:* 

Option        | Description
--------------|-------------
--dir_list  | List of directories which include the required files, one on each line.
--ph_list   | List of directories which include haplotype summary files, one on each line.
--amplicons | Amplicons used in assembly (default: False).


*Settings:*

Option    | Description
----------|-------------
--quiet | Do not write output to console (silence stdout and stderr) (default: False).
--logfile | Name for log file.
--debug   | Print commands but do not run (default: False).

_Example usage:_
```
haphpipe summary_stats --dir_list demo_sra_list.txt --ph_list demo_sra_ph_list.txt --amplicons
```

