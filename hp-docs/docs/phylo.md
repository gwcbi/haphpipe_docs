**hp_phylo** includes phylogenomics stages.

### *multiple_align*
Align consensus sequences using MAFFT ([documentation](https://mafft.cbrc.jp/alignment/software/manual/manual.html)). Input can be a list of directories which contain `final.fna` files or a fasta file, or both (in which case the sequences in the FASTA file are combined with the `final.fna` files retreived before the alignment.
Sequences will be separated by amplicons using a supplied GTF file before alignment (unless the `--alignall` option is specified). This module may also be used to separate files by amplicons (without aligning) by specifying the `--fastaonly` option.
Alignments are by default outputted as FASTA files, although PHYLIP (`--phylipout`) or CLUSTAL (`--clustalout`) output options are also available.
Many options from MAFFT are available in this module. Please refer to the MAFFT documentation above for information about these options.

**Usage:**

`haphpipe multiple_align [MAFFT OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> --dir_list <TXT> --ref_gtf <GTF> [--outdir]`

**(or):**

`hp_multiple_align [MAFFT OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> --dir_list <TXT> --ref_gtf <GTF> [--outdir]`

*Output files:* alignment files in FASTA format (default), one per amplicon (or one `alignment.fasta` file if using `--alignall` option)

Note: MAFFT stores intermediate files in a temporary directory located in /tmp. More information is available [here](https://mafft.cbrc.jp/alignment/software/mpi.html#TMPDIR).

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

### *build_tree*
Phylogeny reconstruction with RAxML ([documentation](https://cme.h-its.org/exelixis/resource/download/NewManual.pdf)). Input is an alignment (FASTA or PHYLIP format). Output is a tree file (TRE format).


**Usage:**

`haphpipe build_tree [RAxML OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> --output_name <TXT> [--outdir]`

**(or):**

`hp build_tree [RAxML OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> --output_name <TXT> [--outdir]`

*Output files:*

File   | Description |
---------|-------------|
RaxML_info.exampleRun              | This file is written regardless of the command line option. It contains information about the model and algorithm used. If --InitialRearrangement is called, it will indicate the rearrangement setting used 
RAxML_log.exampleRun               | This file prints out the time, likelihood value of the current tree, and number of the checkpoint file (if called) after each iteration of the search algorithm. Not generated in the case of multiple bootstraps. 
RAxML_result.exampleRun            | Unless multiple bootstraps are executes, this file is written after each iteration of the search algorithm. It contians the final tree topology of the current run
RAxML_parsimonyTree.exampleRun     | If a starting tree is not specified by --UserStartingTree, this file will contain the randomized parsimony starting tree
RAxML_randomTree.exampleRun        | If --rand_starting_tree if called, this file will contain the completely random starting tree
RAxML_checkpoint.exampleRun.checkpointNumber          | Generated if --print_intermediate_trees is called
RAxML_bootstrap.exampleRun         | Consolidates final bootstrapped trees if called with --NumberofRuns
RAxML_bipartitions.exampleRun      | Contain the input tree with confidence values at nodes if --algo_option is called
RAxML_bipartitionsBranchLabels.exampleRun             | Support values are displayed as Newick branch labels rather than node labels
AxML_bipartitionFrequencies.exampleRun                | If --algo_optio m is called, this file contains the pair-wise bipartition frequencies of all trees 
RAxML_perSiteLLs.exampleRun        | This file contains contains the per–site log likelihood scores. Only generated if --algo_option g is called
RAxML_bestTree.exampleRun          | Outputs the best-scoring ML tree
RAxML_distances.exampleRun         | Contains the pair-wise ML-based distances between taxonpairs. This file is only generated when --algo_option x option is called.

*Input/Output Arguments:* 

Option   | Description 
---------|-------------

  --seqs SEQS                  |  Input alignment in PHYLIP or FASTA format
  --output_name build_tree.tre |  Run name for trees
  --model GTRGAMMA             |  Substitution Model
  --outdir treedir             |  Output directory
 
*RAxML Options:*

Option   | Description 
---------|-------------

  --outgroup                  |  outgrpup for tree
  --parsimony seed            |  Parsimony Random Seed
  --wgtFile                   |  Column weight file name to assign individual weights to each column of the alignment
  --secsub                    |  Specify secondary structure substitution models, must also include file defining the secondary structure 
  --bootstrap                 |  bootstrapRandomNumberSeed for non-parametric bootstrapping
  --parsimony seed            |  Parsimony Random Seed
  --wgtFile                   |  Column weight file name to assign individual weights to each column of the alignment
  --secsub                    |  Specify secondary structure substitution models, must also include file defining the secondary structure 
  --bootstrap_threshold       | Threshold for bootstopping criteria
  --numCat                    | Number of distinct rate categories for RAxML when model of rate heterogeneity is set to CAT
  --rand_starting_tree        | ML optimization from random starting tree
  --convergence_criterion     |ML search convergence criterion
  --likelihoodEpsilon         | Set model optimization precision in log likelihood units for final optimization of tree topology
  --excludeFileName           | File contains specifications of alignment positions to be excluded
  --algo_option               | Select what kind of algorithm RAxML shall execute
  --cat_model                 | Enable ML tree searches under CAT model for very large trees
  --groupingFile              | File name of a multifurcating constraint tree
  --placementThreshold        | Threshold value for ML­based evolutionary placement algorithm heuristics
  --disable_pattern_compression | Disable pattern compression
  --InitialRearrangement      | Radius for pruned sub-tree re-insertion
  --posteriori                | posteriori bootstopping analysis
  --print_intermediate_trees  | Print out a couple of intermediate trees
  --majorityrule              | Compute majority rule consensus tree
  --print_branch_length       | Bootstrapped trees should be printed with branch lengths
  --ICTCmetrics               | Compute the TC and IC metrics on a consensus tree
  --partition_branch_length   | Switch on estimation of individual per­partition branch lengths
  --disable_check             | Disable check for completely undetermined sequence in alignment
  --AAmodel                   | Specify the file name of a user­defined AA (Protein) substitution model
  --multiplemodelFile         | Specify the file name which contains the assignment of models to alignment partitions for multiple models of substitution
  --binarytree                | Specify the file name of a binary constraint tree
  --BinaryParameterFile       | Specify the file name of a binary model parameter file that has previously been generated with RAxML using the ­f e tree evaluation option.
  --SecondaryStructure        | Specify the name of a secondary structure file
  --UserStartingTree          | Specifies a user starting tree file name which must be in Newick format
  --median_GAMMA              | Use the median for the discrete approximation of the GAMMA model of rateheterogeneity
  --implement_SEV             | Try to save memory by using SEV­based implementation
  --version_info              | Display version information
  --rate_heterogeneity        | Disable rate heterogeneity among site model and use one without rate heterogeneity instead
  --directory                 | Full directory of output file
  --window                    | Sliding window size for leave­one­out site­specific placement bias algorithm
  --RapidBootstrapNumSeed     | Specify an integer number (random seed) and turn on rapid bootstrapping
  --random_addition           | RAxML will only do a randomized stepwise addition order parsimony tree reconstruction without performing any additional SPRs
  --starting_tree             | Only for computing parsimony
  --quartetGroupingFileName   | Pass a quartet grouping file name defining four groups from which to draw quartets
  --multipleTreeFile          | Specify the file name of a file containing multiple trees e.g. from a bootstrap that shall be used to draw bipartition values onto a tree provided with ­t.
  --NumberofRuns              | Specify the number of alternative runs on distinct starting trees
  --mesquite                  | Print output files that can be parsed by Mesquite
  --silent                    | Disables printout of warnings related to identical sequences and entirely undetermined sites in the alignment
  --noseqcheck                | Disables checking the input MSA for identical sequences and entirely undetermined sites
  --nobfgs                    | Disables automatic usage of BFGS method to optimize GTR rates on unpartitioned DNA datasets
  --asccorrlewis              | The standard correction by Paul Lewis
  --asccorrfelsenstein        | A correction introduced by Joe Felsenstein
  --asccorrstamatakis         | A correction introduced by Stamatakis
  --flagcheck                 | RAxML will only check if all command line flags specifed are available
  --autoprotml                | When using automatic protein model selection you can chose the criterion for selecting these models
  --autoprotbic               | When using automatic protein model selection you can chose the criterion for selecting these models
  --autoprotaic'              | When using automatic protein model selection you can chose the criterion for selecting these models
  --autoprotaicc              | When using automatic protein model selection you can chose the criterion for selecting these models
  --epaPlaceNum               | Specify the number of potential placements you want to keep for each read in the EPA algorithm
  --epaProbThreshold          | Specify a percent threshold for including potential placements of a read depending on the maximum placement weight for this read
  --epaLikelihood             | Specify an accumulated likelihood weight threshold
  --JC69                      | Specify that all DNA partitions will evolve under the Jukes­Cantor model
  --K80                       | Specify that all DNA partitions will evolve under the K80 model
  --HKY85                     | Specify that all DNA partitions will evolve under the HKY85 model
  --BootstrapPerm             | Specify the number of permutations to be conducted for the bootstopping/bootstrap convergence test; minimum 100
  --quartetswithoutreplacement | Specify that quartets are randomly subsampled, but without replacement
  --printidenticalsequences    | Specify that RAxML shall automatically generate a .reduced alignment with all undetermined columns removed')
  --option_help                | Display Help

*Options:*

Option   | Description |
---------|-------------|

 --keep_tmp tempdir | Keep temporary directory  
 --quiet            |     Do not write output to console (silence stdout and  stderr) (default: False)
 --logfile LOGFILE  |   Name for log file (output)
 --debug            |  Print commands but do not run (default: False)
 
