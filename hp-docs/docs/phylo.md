**Phylo** includes phylogenomics stages.

### *multiple_align_*
Align consensus sequences using MAFFT ([documentation](https://mafft.cbrc.jp/alignment/software/manual/manual.html)). Input can be a list of directories which contain `final.fna` files or a fasta file, or both (in which case the sequences in the FASTA file are combined with the `final.fna` files retreived before the alignment.
Sequences will be separated by amplicons using a supplied GTF file before alignment (unless the `--alignall` option is specified). This module may also be used to separate files by amplicons (without aligning) by specifying the `--fastaonly` option.
Alignments are by default outputted as FASTA files, although PHYLIP (`--phylipout`) or CLUSTAL (`--clustalout`) output options are also available.
Many options from MAFFT are available in this module. Please refer to the MAFFT documentation above for information about these options.

**Usage:**

`haphpipe multiple_align [MAFFT OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> --dir_list <TXT> --ref_gtf <GTF> [--outdir]`

**(or):**

`hp_multiple_align [MAFFT OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> --dir_list <TXT> --ref_gtf <GTF> [--outdir]`

*Output files:* alignment files in FASTA format (default), one per amplicon (or one `alignment.fasta` file if using `--alignall` option)

*Input/Output Arguments:* 

Option   | Description 
---------|-------------
--seqs SEQS |           FASTA file with sequences to be aligned
  --dir_list DIR_LIST |   List of directories which include a final.fna file,one on each line
  --ref_gtf REF_GTF |    Reference GTF file
  --out_align OUT_ALIGN | Name for alignment file
  --nuc |                Assume nucleotide (default: False)
  --amino  |             Assume amino (default: False)
  --clustalout   |       Clustal output format (default: False)
  --phylipout      |     PHYLIP output format (default: False)
  --inputorder      |    Output order same as input (default: False)
  --reorder         |    Output order aligned (default: False)
  --treeout          |   Guide tree is output to the input.tree file (default:False)
  --quiet_mafft      |   Do not report progress (default: False)
  --outdir OUTDIR    |   Output directory


*MAFFT Options:*

Option   | Description |
---------|-------------|
--algo ALGO           |Use different algorithm in command: linsi, ginsi,   einsi, fftnsi, fftns, nwns, nwnsi
  --auto              |  Automatically select algorithm (default: False)
  --sixmerpair         | Calculate distance based on shared 6mers, on by default (default: False)
  --globalpair         | Use Needleman-Wunsch algorithm (default: False)
  --localpair           | Use Smith-Waterman algorithm (default: False)
  --genafpair           | Use local algorithm with generalized affine gap cost (default: False)
  --fastapair         |  Use FASTA for pairwise alignment (default: False)
  --weighti WEIGHTI    | Weighting factor for consistency term
  --retree RETREE     |  Number of times to build guide tree
  --maxiterate MAXITERATE | Number of cycles for iterative refinement
  --noscore          |   Do not check alignment score in iterative alignment (default: False)
  --memsave           |  Use Myers-Miller algorithm (default: False)
  --parttree         |   Use fast tree-building method with 6mer distance (default: False)
  --dpparttree        |  Use PartTree algorithm with distances based on DP (default: False)
  --fastaparttree     |  Use PartTree algorithm with distances based on FASTA (default: False)
  --partsize PARTSIZE  | Number of partitions for PartTree
  --groupsize GROUPSIZE | Max number of sequences for PartTree
  
*MAFFT Parameters:*

Option   | Description |
---------|-------------|
 --lop LOP            | Gap opening penalty
  --lep LEP            | Offset value
  --lexp LEXP          | Gap extension penalty
  --LOP LOP            | Gap opening penalty to skip alignment
  --LEXP LEXP          | Gap extension penalty to skip alignment
  --bl BL              | BLOSUM matrix: 30, 45, 62, or 80
  --jtt JTT            | JTT PAM number >0
  --tm TM              | Transmembrane PAM number >0
  --aamatrix AAMATRIX  | Path to user-defined AA scoring matrix
  --fmodel             | Incorporate AA/nuc composition info into scoring matrix (default: False)


*Options:*

Option   | Description |
---------|-------------|
--ncpu NCPU      |     Number of CPU to use (default: 1)
  --quiet          |     Do not write output to console (silence stdout and  stderr) (default: False)
  --logfile LOGFILE  |   Name for log file (output)
  --debug             |  Print commands but do not run (default: False)
  --fastaonly         |  Output fasta files separated by region but do not     align (default: False)
  --alignall         |   Do not separate files by region, align entire file (default: False)


_Example usage:_


```
haphpipe multiple_align --dir_list demo_sra_list.txt --ref_gtf HIV_B.K03455.HXB2.gtf --phylipout --logfile demo_multiple_align.log
```
