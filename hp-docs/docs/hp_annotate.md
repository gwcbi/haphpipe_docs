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

### *annotate_ from_ref*
Annotate consensus sequence from reference annotation. Input is JSON file from _pairwise_align_ and reference GTF file. 

**Usage:**

`haphpipe annotate_from_ref [OPTIONS] [SETTINGS] --haplotypes_fa <best.fas> [--outdir]`

**(or):**

`hp_annotate_from_ref [OPTIONS] [SETTINGS] --haplotypes_fa <best.fas> [--outdir]`

_Example usage:_

