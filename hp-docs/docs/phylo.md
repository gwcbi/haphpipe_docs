**hp_phylo** includes phylogenomics stages.

### *multiple_align*
Align consensus sequences using MAFFT ([documentation](https://mafft.cbrc.jp/alignment/software/manual/manual.html)). Input can be a list of directories which contain `final.fna` and/or `ph_haplotypes.fna` files or a fasta file, or both (in which case the sequences in the FASTA file are combined with the `final.fna` and/or `ph_haplotypes.fna` files retreived before the alignment.
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
  --dir_list DIR_LIST |   List of directories which include either a final.fna or ph_haplotypes.fna file, one on each line
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
--ncpu NCPU      |     Number of CPU to use (default: 1)|
  --quiet          |     Do not write output to console (silence stdout and  stderr) (default: False)|
  --logfile LOGFILE  |   Name for log file (output)|
  --debug             |  Print commands but do not run (default: False)|
  --fastaonly         |  Output fasta files separated by region but do not     align (default: False)|
  --alignall         |   Do not separate files by region, align entire file (default: False)|
  
_Example usage:_


```
haphpipe multiple_align --dir_list demo_sra_list.txt --ref_gtf HIV_B.K03455.HXB2.gtf --phylipout --logfile demo_multiple_align.log
```


### *model_test*
Select the best-fit model of evolution from an alignment file using ModelTest-NG ([documentation](https://github.com/ddarriba/modeltest/wiki)). Input is an alignment in FASTA or PHYLIP format. Output is ModelTest-NG results (text file) containing information for the best performing models.

**Usage:**

`haphpipe model_test [ModelTest-NG OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> [--outdir]`

**(or):**

`hp_model_test [ModelTest-NG OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> [--outdir]`

*Output files:* ModelTest-NG output file (`modeltest_results.out`).

*Input/Output Arguments:* 

Option   | Description |
---------|-------------|
--seqs SEQS |   Alignment in FASTA or PHYLIP format
--outname | Name for output file
--outdir OUTDIR    |   Output directory


*ModelTest-NG Options:*

Option   | Description |
---------|-------------|
--data_type|Data type: nt or aa (default: nt)
--partitions|Partitions file
--seed|Seed for random number generator
--topology TOPOLOGY|Starting topology: ml, mp, fixed-ml-jc, fixed-ml-gtr, fixed-mp, random, or user (default: ml)
--utree|User-defined starting tree
--force| Force output overriding (default: False)
--asc_bias ASC_BIAS|Ascertainment bias correction: lewis, felsenstein, or stamatakis
--frequencies FREQUENCIES|Candidate model frequencies: e (estimated) or f (fixed)
--het HET| Set rate heterogeneity: u (uniform), i (invariant sites +I), g (gamma +G), or f (bothinvariant sites and gamma +I+G)
--models MODELS |  Text file with candidate models, one per line
--schemes SCHEMES|Number of predefined DNA substitution schemes evaluated: 3, 5, 7, 11, or 203
--template TEMPLATE|Set candidate models according to a specified tool: raxml, phyml, mrbayes, or paup

*Options:*

Option   | Description |
---------|-------------|
--ncpu NCPU      |     Number of CPU to use (default: 1)|
  --quiet          |     Do not write output to console (silence stdout and  stderr) (default: False)|
  --logfile LOGFILE  |   Name for log file (output)|
  --debug             |  Print commands but do not run (default: False)|
--keep_tmp|Keep temporary directory (default: False)  

_Example usage:_


```
haphpipe model_test --seqs multiple_align/alignment.fasta
```



### *build_tree*
Phylogeny reconstruction with RAxML ([documentation](https://cme.h-its.org/exelixis/resource/download/NewManual.pdf)). Input is an alignment (FASTA or PHYLIP format). Output is a tree file (TRE format).
Please see the RAxML documentation for a full description of RAxML options. For convenience, we have included an option `--run_full_analysis` which will automatically find the best maximum likelihood tree, complete bootstrapping, and then merge output together for a final tree.

**Usage:**

`haphpipe build_tree [RAxML OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> --output_name <TXT> [--outdir]`

**(or):**

`hp build_tree [RAxML OPTIONS] [HAPHPIPE OPTIONS] --seqs <FASTA> --output_name <TXT> [--outdir]`

*Output files:*

File   | Description |
---------|-------------|
RaxML_info.build_tree.tre              | This file is written regardless of the command line option. It contains information about the model and algorithm used. If --InitialRearrangement is called, it will indicate the rearrangement setting used 
RAxML_log.build_tree.tre               | This file prints out the time, likelihood value of the current tree, and number of the checkpoint file (if called) after each iteration of the search algorithm. Not generated in the case of multiple bootstraps. 
RAxML_result.build_tree.tre            | Unless multiple bootstraps are executes, this file is written after each iteration of the search algorithm. It contians the final tree topology of the current run
RAxML_parsimonyTree.build_tree.tre     | If a starting tree is not specified by --UserStartingTree, this file will contain the randomized parsimony starting tree
RAxML_randomTree.build_tree.tre        | If --rand_starting_tree if called, this file will contain the completely random starting tree
RAxML_checkpoint.build_tree.tre.checkpointNumber          | Generated if --print_intermediate_trees is called
RAxML_bootstrap.build_tree.tre         | Consolidates final bootstrapped trees if called with --NumberofRuns
RAxML_bipartitions.build_tree.tre      | Contain the input tree with confidence values at nodes if --algo_option is called
RAxML_bipartitionsBranchLabels.build_tree.tre             | Support values are displayed as Newick branch labels rather than node labels
AxML_bipartitionFrequencies.build_tree.tre                | If --algo_optio m is called, this file contains the pair-wise bipartition frequencies of all trees 
RAxML_perSiteLLs.build_tree.tre        | This file contains contains the per–site log likelihood scores. Only generated if --algo_option g is called
RAxML_bestTree.build_tree.tre          | Outputs the best-scoring ML tree
RAxML_distances.build_tree.tre         | Contains the pair-wise ML-based distances between taxonpairs. This file is only generated when --algo_option x option is called.

*Input/Output Arguments:* 

Option   | Description |
---------|-------------|
  --seqs SEQS                  |  Input alignment in PHYLIP or FASTA format|
  --output_name NAME |  Run name for trees (default: build_tree.tre)|
  --model MODEL             |  Substitution Model (default: GTRGAMMAIX)|
  --outdir OUTDIR              |  Output directory (default: .)|
 
*RAxML Options:*

Option   | Description 
---------|-------------
--run_full_analysis | Run bootstrap search and find best ML tree |
  --outgroup                  |  outgrpup for tree|
  --parsimony_seed            |  Parsimony Random Seed|
  --wgtFile                   |  Column weight file name to assign individual weights to each column of the alignment|
  --secsub                    |  Specify secondary structure substitution models, must also include file defining the secondary structure |
  --bootstrap                 |  bootstrapRandomNumberSeed for non-parametric bootstrapping|
  --wgtFile                   |  Column weight file name to assign individual weights to each column of the alignment|
  --secsub                    |  Specify secondary structure substitution models, must also include file defining the secondary structure |
  --bootstrap_threshold       | Threshold for bootstopping criteria|
  --numCat                    | Number of distinct rate categories for RAxML when model of rate heterogeneity is set to CAT|
  --rand_starting_tree        | ML optimization from random starting tree|
  --convergence_criterion     |ML search convergence criterion|
  --likelihoodEpsilon         | Set model optimization precision in log likelihood units for final optimization of tree topology|
  --excludeFileName           | File contains specifications of alignment positions to be excluded|
  --algo_option               | Select what kind of algorithm RAxML shall execute|
  --cat_model                 | Enable ML tree searches under CAT model for very large trees|
  --groupingFile              | File name of a multifurcating constraint tree|
  --placementThreshold        | Threshold value for ML­based evolutionary placement algorithm heuristics|
  --disable_pattern_compression | Disable pattern compression|
  --InitialRearrangement      | Radius for pruned sub-tree re-insertion|
  --posteriori                | posteriori bootstopping analysis|
  --print_intermediate_trees  | Print out a couple of intermediate trees|
  --majorityrule              | Compute majority rule consensus tree|
  --print_branch_length       | Bootstrapped trees should be printed with branch lengths|
  --ICTCmetrics               | Compute the TC and IC metrics on a consensus tree|
  --partition_branch_length   | Switch on estimation of individual per­partition branch lengths|
  --disable_check             | Disable check for completely undetermined sequence in alignment|
  --AAmodel                   | Specify the file name of a user­defined AA (Protein) substitution model|
  --multiplemodelFile         | Specify the file name which contains the assignment of models to alignment partitions for multiple models of substitution|
  --binarytree                | Specify the file name of a binary constraint tree|
  --BinaryParameterFile       | Specify the file name of a binary model parameter file that has previously been generated with RAxML using the ­f e tree evaluation option.|
  --SecondaryStructure        | Specify the name of a secondary structure file|
  --UserStartingTree          | Specifies a user starting tree file name which must be in Newick format|
  --median_GAMMA              | Use the median for the discrete approximation of the GAMMA model of rateheterogeneity|
  --implement_SEV             | Try to save memory by using SEV­based implementation|
  --version_info              | Display version information|
  --rate_heterogeneity        | Disable rate heterogeneity among site model and use one without rate heterogeneity instead|
  --directory                 | Full directory of output file|
  --window                    | Sliding window size for leave­one­out site­specific placement bias algorithm|
  --RapidBootstrapNumSeed     | Specify an integer number (random seed) and turn on rapid bootstrapping|
  --random_addition           | RAxML will only do a randomized stepwise addition order parsimony tree reconstruction without performing any additional SPRs|
  --starting_tree             | Only for computing parsimony|
  --quartetGroupingFileName   | Pass a quartet grouping file name defining four groups from which to draw quartets|
  --multipleTreeFile          | Specify the file name of a file containing multiple trees e.g. from a bootstrap that shall be used to draw bipartition values onto a tree provided with ­t.|
  --NumberofRuns              | Specify the number of alternative runs on distinct starting trees|
  --mesquite                  | Print output files that can be parsed by Mesquite|
  --silent                    | Disables printout of warnings related to identical sequences and entirely undetermined sites in the alignment|
  --noseqcheck                | Disables checking the input MSA for identical sequences and entirely undetermined sites|
  --nobfgs                    | Disables automatic usage of BFGS method to optimize GTR rates on unpartitioned DNA datasets|
  --asccorrlewis              | The standard correction by Paul Lewis|
  --asccorrfelsenstein        | A correction introduced by Joe Felsenstein|
  --asccorrstamatakis         | A correction introduced by Stamatakis|
  --flagcheck                 | RAxML will only check if all command line flags specifed are available|
  --autoprotml                | When using automatic protein model selection you can chose the criterion for selecting these models|
  --autoprotbic               | When using automatic protein model selection you can chose the criterion for selecting these models|
  --autoprotaic'              | When using automatic protein model selection you can chose the criterion for selecting these models|
  --autoprotaicc              | When using automatic protein model selection you can chose the criterion for selecting these models|
  --epaPlaceNum               | Specify the number of potential placements you want to keep for each read in the EPA algorithm|
  --epaProbThreshold          | Specify a percent threshold for including potential placements of a read depending on the maximum placement weight for this read|
  --epaLikelihood             | Specify an accumulated likelihood weight threshold|
  --JC69                      | Specify that all DNA partitions will evolve under the Jukes­Cantor model|
  --K80                       | Specify that all DNA partitions will evolve under the K80 model|
  --HKY85                     | Specify that all DNA partitions will evolve under the HKY85 model|
  --BootstrapPerm             | Specify the number of permutations to be conducted for the bootstopping/bootstrap convergence test; minimum 100|
  --quartetswithoutreplacement | Specify that quartets are randomly subsampled, but without replacement|
  --printidenticalsequences    | Specify that RAxML shall automatically generate a .reduced alignment with all undetermined columns removed')|
  --option_help                | Display Help|

*Options:*

Option   | Description |
---------|-------------|
 --keep_tmp  | Keep temporary directory  |
 --quiet            |     Do not write output to console (silence stdout and  stderr) (default: False)|
 --logfile LOGFILE  |   Name for log file (output)|
 --debug            |  Print commands but do not run (default: False)|
 
_Example usage:_


```
haphpipe build_tree --seqs multiple_align/alignment.fasta --run_full_analysis
```
