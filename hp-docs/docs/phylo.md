**hp_phylo** includes phylogenomics stages.

### Phylo Quick-Start 

HAPHPIPE includes three stages for phylogenomics: *multiple_align*, *model_test*, and *build_tree*. These three stages are sufficient to turn your consensus and/or haplotype sequences from the other stages into a phylogenetic tree! For purposes of this quick-start guide, we will demonstrate the stages to create a tree from HIV _pol_ consensus sequences. 

**Step 1: Alignment**

After running either of the assembly pipelines, `final.fna` files will be located in directories named `./<SampleID>/haphpipe_assemble_0[1|2]`. For the *multiple_align* stage, we need to create a list of all of these directories. We can do so easily with one command (shown for `haphpipe_assemble_01` output:
```
ls -d ./SRR*/haphpipe_assemble_01 > ./dir_list.txt
```

Now, we will align all of these `final.fna` files:
```
haphpipe multiple_align --dir_list dir_list.txt --ref_gtf refs/HIV_B.K03455.HXB2.gtf 
```

The output will be located in a new directory, `hp_multiple_align`. The alignment of _pol_ sequences is the file `alignment_region00.fasta`.

**Step 2: Model Selection**

Now, we will use the *model_test* stage to determine the best-fit evolutionary model for our data. This is an input to the tree building stage. We will use this command to generate best-fit models available in RAxML:
```
haphpipe model_test --seqs hp_multiple_align/alignment_region00.fasta --run_id alignment_region00 --template raxml
```

The ModelTest output will be written to a file called `modeltest_results.out` and a summary of all the best models will be written to `modeltest_results_summary.tsv`. Examples of both are below.

<details>
  <summary>ModelTest-NG Output</summary>
  
  ```
  --------------------------------------------------------------------------------
ModelTest-NG vx.y.z

Input data:
  MSA:        multiple_align/alignment.fasta
  Tree:       Maximum likelihood
    file:           -
  #taxa:            6
  #sites:           1975
  #patterns:        114
  Max. thread mem:  0 MB

Output:
  Log:           /var/folders/lv/dqdkd8957_3fv6yxsyfsvn0r0000gn/T/tmpHP_model_testud4sc4ya/samp12_modeltest_results.log
  Starting tree: /var/folders/lv/dqdkd8957_3fv6yxsyfsvn0r0000gn/T/tmpHP_model_testud4sc4ya/samp12_modeltest_results.tree
  Results:       /var/folders/lv/dqdkd8957_3fv6yxsyfsvn0r0000gn/T/tmpHP_model_testud4sc4ya/samp12_modeltest_results.out

Selection options:
  # dna schemes:      11
  # dna models:       88
  include model parameters:
    Uniform:         true
    p-inv (+I):      true
    gamma (+G):      true
    both (+I+G):     true
    free rates (+R): false
    fixed freqs:     true
    estimated freqs: true
    #categories:     4
  gamma rates mode:   mean
  asc bias:           none
  epsilon (opt):      0.01
  epsilon (par):      0.05
  keep branches:      false

Additional options:
  verbosity:        very low
  threads:          1/6
  RNG seed:         12345
  subtree repeats:  enabled
--------------------------------------------------------------------------------

BIC       model              K            lnL          score          delta    weight
--------------------------------------------------------------------------------
       1  HKY                4     -5380.7535     10860.1551         0.0000    0.7049
       2  TrN                5     -5378.2013     10862.6391         2.4839    0.2036
       3  TPM1uf             5     -5379.8699     10865.9763         5.8211    0.0384
       4  TPM3uf             5     -5380.7175     10867.6716         7.5164    0.0164
       5  HKY+G4             5     -5380.8440     10867.9246         7.7694    0.0145
       6  HKY+I              5     -5381.4917     10869.2200         9.0648    0.0076
       7  TIM3               6     -5378.1637     10870.1523         9.9971    0.0048
       8  TrN+G4             6     -5378.3069     10870.4387        10.2836    0.0041
       9  TPM2uf+G4          6     -5378.9908     10871.8064        11.6512    0.0021
      10  TrN+I              6     -5379.0181     10871.8610        11.7059    0.0020
--------------------------------------------------------------------------------
Best model according to BIC
---------------------------
Model:              HKY
lnL:                -5380.7535
Frequencies:        0.3695 0.1709 0.2219 0.2377
Subst. Rates:       1.0000 2.4741 1.0000 1.0000 2.4741 1.0000 
Inv. sites prop:    -
Gamma shape:        -
Score:              10860.1551
Weight:             0.7049
---------------------------
Parameter importances
---------------------------
P.Inv:              0.0100
Gamma:              0.0216
Gamma-Inv:          0.0002
Frequencies:        1.0000
---------------------------
Model averaged estimates
---------------------------
P.Inv:              0.0215
Alpha:              94.2337
Alpha-P.Inv:        91.7960
P.Inv-Alpha:        0.0214
Frequencies:        0.3700 0.1702 0.2226 0.2372 

Commands:
  > phyml  -i multiple_align/alignment.fasta -m 010010 -f m -v 0 -a 0 -c 1 -o tlr
  > raxmlHPC-SSE3 -s multiple_align/alignment.fasta -c 1 -m GTRCATX -n EXEC_NAME -p PARSIMONY_SEED
  > raxml-ng --msa multiple_align/alignment.fasta --model HKY
  > paup -s multiple_align/alignment.fasta
  > iqtree -s multiple_align/alignment.fasta -m HKY

AIC       model              K            lnL          score          delta    weight
--------------------------------------------------------------------------------
       1  TrN                5     -5378.2013     10784.4026         0.0000    0.2246
       2  TIM2+G4            7     -5376.4972     10784.9944         0.5919    0.1671
       3  TIM3               6     -5378.1637     10786.3274         1.9249    0.0858
       4  TIM2+I             7     -5377.2570     10786.5141         2.1115    0.0782
       5  TrN+G4             6     -5378.3069     10786.6138         2.2113    0.0744
       6  TIM1+G4            7     -5377.3549     10786.7098         2.3073    0.0709
       7  HKY                4     -5380.7535     10787.5069         3.1044    0.0476
       8  TPM1uf             5     -5379.8699     10787.7398         3.3372    0.0423
       9  TPM2uf+G4          6     -5378.9908     10787.9815         3.5790    0.0375
      10  TrN+I              6     -5379.0181     10788.0362         3.6336    0.0365
--------------------------------------------------------------------------------
Best model according to AIC
---------------------------
Model:              TrN
lnL:                -5378.2013
Frequencies:        0.3720 0.1679 0.2249 0.2352
Subst. Rates:       1.0000 2.2008 1.0000 1.0000 3.1065 1.0000 
Inv. sites prop:    -
Gamma shape:        -
Score:              10784.4026
Weight:             0.2246
---------------------------
Parameter importances
---------------------------
P.Inv:              0.1262
Gamma:              0.3943
Gamma-Inv:          0.0384
Frequencies:        1.0000
---------------------------
Model averaged estimates
---------------------------
P.Inv:              0.0216
Alpha:              93.2595
Alpha-P.Inv:        94.4193
P.Inv-Alpha:        0.0216
Frequencies:        0.3718 0.1684 0.2238 0.2359 

Commands:
  > phyml  -i multiple_align/alignment.fasta -m 010020 -f m -v 0 -a 0 -c 1 -o tlr
  > raxmlHPC-SSE3 -s multiple_align/alignment.fasta -c 1 -m GTRCATX -n EXEC_NAME -p PARSIMONY_SEED
  > raxml-ng --msa multiple_align/alignment.fasta --model TrN
  > paup -s multiple_align/alignment.fasta
  > iqtree -s multiple_align/alignment.fasta -m TrN

AICc      model              K            lnL          score          delta    weight
--------------------------------------------------------------------------------
       1  TrN                5     -5378.2013     10784.4026         0.0000    0.2246
       2  TIM2+G4            7     -5376.4972     10784.9944         0.5919    0.1671
       3  TIM3               6     -5378.1637     10786.3274         1.9249    0.0858
       4  TIM2+I             7     -5377.2570     10786.5141         2.1115    0.0782
       5  TrN+G4             6     -5378.3069     10786.6138         2.2113    0.0744
       6  TIM1+G4            7     -5377.3549     10786.7098         2.3073    0.0709
       7  HKY                4     -5380.7535     10787.5069         3.1044    0.0476
       8  TPM1uf             5     -5379.8699     10787.7398         3.3372    0.0423
       9  TPM2uf+G4          6     -5378.9908     10787.9815         3.5790    0.0375
      10  TrN+I              6     -5379.0181     10788.0362         3.6336    0.0365
--------------------------------------------------------------------------------
Best model according to AICc
---------------------------
Model:              TrN
lnL:                -5378.2013
Frequencies:        0.3720 0.1679 0.2249 0.2352
Subst. Rates:       1.0000 2.2008 1.0000 1.0000 3.1065 1.0000 
Inv. sites prop:    -
Gamma shape:        -
Score:              10784.4026
Weight:             0.2246
---------------------------
Parameter importances
---------------------------
P.Inv:              0.1262
Gamma:              0.3943
Gamma-Inv:          0.0384
Frequencies:        1.0000
---------------------------
Model averaged estimates
---------------------------
P.Inv:              0.0216
Alpha:              93.2595
Alpha-P.Inv:        94.4193
P.Inv-Alpha:        0.0216
Frequencies:        0.3718 0.1684 0.2238 0.2359 

Commands:
  > phyml  -i multiple_align/alignment.fasta -m 010020 -f m -v 0 -a 0 -c 1 -o tlr
  > raxmlHPC-SSE3 -s multiple_align/alignment.fasta -c 1 -m GTRCATX -n EXEC_NAME -p PARSIMONY_SEED
  > raxml-ng --msa multiple_align/alignment.fasta --model TrN
  > paup -s multiple_align/alignment.fasta
  > iqtree -s multiple_align/alignment.fasta -m TrN
Done

  ```
  
</details>

<details>
  <summary>ModelTest-NG Output</summary>
  
```
File	Criteria	Best Model
multiple_align/alignment.fasta	BIC	HKY
multiple_align/alignment.fasta	AIC	TrN
multiple_align/alignment.fasta	AICc	TrN
```
  
</details>

**Step 3: Build a Tree**

Now, we will use *build_tree* to build our tree! You should use the best model outputted in *model_test* for the `--model` argument (here we are using GTRGAMMAX). The `--run_full_analysis` option will automatically run a full maximum likelihood & bootstrapping analysis for us:
```
haphpipe build_tree --seqs hp_multiple_align/alignment_region00.fasta --run_full_analysis --model GTRGAMMAX
```

The output will be written to a new directory, `hp_build_tree`. The best tree file from RAxML will be outputted as `RAxML_bipartitionsBranchLabels.build_tree.tre`. This tree can then be annotated in programs such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or [iTOL](https://itol.embl.de).

**Phylogenomics Pipelines**

For users who would like to build a full pipeline to run assembly and phylogenetics stages in one go, we recommend adapting the demo pipeline (`haphpipe_demo`) for this purpose. See the [demo page](https://gwcbi.github.io/haphpipe_docs/demos/) for more details.

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

Option   | RAxML Equivalent |Description |
---------|------------------|----        |
  --seqs SEQS        | -s          |  Input alignment in PHYLIP or FASTA format|
  --output_name NAME | -n          |  Run name for trees (default: build_tree.tre)|
  --model MODEL      | -m          |  Substitution Model (default: GTRGAMMAIX)|
  --outdir OUTDIR    | -w          |  Output directory (default: .)|
 
*RAxML Options:*

Option   | RAxML Equivalent |Description |
---------|------------------|----        |
--run_full_analysis | Run bootstrap search and find best ML tree |
  --outgroup                  | -o  | outgrpup for tree|
  --parsimony_seed            | -p  | Parsimony Random Seed|
  --wgtFile                   | -a  | Column weight file name to assign individual weights to each column of the alignment|
  --secsub                    | -A  | Specify secondary structure substitution models, must also include file defining the secondary structure |
  --bootstrap                 | -b  | bootstrapRandomNumberSeed for non-parametric bootstrapping|
  --bootstrap_threshold       | -B  | Threshold for bootstopping criteria|
  --numCat                    | -c  | Number of distinct rate categories for RAxML when model of rate heterogeneity is set to CAT|
  --rand_starting_tree        | -d  | ML optimization from random starting tree|
  --convergence_criterion     | -D  | ML search convergence criterion|
  --likelihoodEpsilon         | -e  | Set model optimization precision in log likelihood units for final optimization of tree topology|
  --excludeFileName           | -E  | File contains specifications of alignment positions to be excluded|
  --algo_option               | -f  | Select what kind of algorithm RAxML shall execute|
  --cat_model                 | -F  | Enable ML tree searches under CAT model for very large trees|
  --groupingFile              | -g  | File name of a multifurcating constraint tree|
  --placementThreshold        | -G  | Threshold value for ML­based evolutionary placement algorithm heuristics|
  --disable_pattern_compression | -H |Disable pattern compression|
  --InitialRearrangement      | -i  | Radius for pruned sub-tree re-insertion|
  --posteriori                | -I  | posteriori bootstopping analysis|
  --print_intermediate_trees  | -j  | Print out a couple of intermediate trees|
  --majorityrule              | -J  | Compute majority rule consensus tree|
  --print_branch_length       | -k  | Bootstrapped trees should be printed with branch lengths|
  --ICTCmetrics               | -L  | Compute the TC and IC metrics on a consensus tree|
  --partition_branch_length   | -M  | Switch on estimation of individual per­partition branch lengths|
  --disable_check             | -O  | Disable check for completely undetermined sequence in alignment|
  --AAmodel                   | -P  | Specify the file name of a user­defined AA (Protein) substitution model|
  --multiplemodelFile         | -q  | Specify the file name which contains the assignment of models to alignment partitions for multiple models of substitution|
  --binarytree                | -r  | Specify the file name of a binary constraint tree|
  --BinaryParameterFile       | -R  | Specify the file name of a binary model parameter file that has previously been generated with RAxML using the ­f e tree evaluation option.|
  --SecondaryStructure        | -S  | Specify the name of a secondary structure file|
  --UserStartingTree          | -t  | Specifies a user starting tree file name which must be in Newick format|
  --median_GAMMA              | -u  | Use the median for the discrete approximation of the GAMMA model of rateheterogeneity|
  --version_info              | -v  | Display version information|
  --rate_heterogeneity        | -V  | Disable rate heterogeneity among site model and use one without rate heterogeneity instead|
  --directory                 | -w  | Full directory of output file|
  --window                    | -W  | Sliding window size for leave­one­out site­specific placement bias algorithm|
  --RapidBootstrapNumSeed     | -x  | Specify an integer number (random seed) and turn on rapid bootstrapping|
  --random_addition           | -X  | RAxML will only do a randomized stepwise addition order parsimony tree reconstruction without performing any additional SPRs|
  --starting_tree             | -y  | Only for computing parsimony|
  --quartetGroupingFileName   | -Y  | Pass a quartet grouping file name defining four groups from which to draw quartets|
  --multipleTreeFile          | -z  | Specify the file name of a file containing multiple trees e.g. from a bootstrap that shall be used to draw bipartition values onto a tree provided with ­t.|
  --NumberofRuns              | -N  | Specify the number of alternative runs on distinct starting trees|
  --mesquite                  | --mesquite                   | Print output files that can be parsed by Mesquite|
  --silent                    | --silent                     | Disables printout of warnings related to identical sequences and entirely undetermined sites in the alignment|
  --noseqcheck                | --no-seq-check               | Disables checking the input MSA for identical sequences and entirely undetermined sites|
  --nobfgs                    | --no-bfgs                    | Disables automatic usage of BFGS method to optimize GTR rates on unpartitioned DNA datasets|
  --epaPlaceNum               | ­­epa­keep­placements=       | Specify the number of potential placements you want to keep for each read in the EPA algorithm|
  --epaProbThreshold          | ­­epa­prob­threshold=        | Specify a percent threshold for including potential placements of a read depending on the maximum placement weight for this read|
  --epaLikelihood             | ­­epa­accumulated­threshold= | Specify an accumulated likelihood weight threshold|
  --HKY85                     | --HKY85                      | Specify that all DNA partitions will evolve under the HKY85 model|
  --BootstrapPerm             | ­­bootstop­perms=            | Specify the number of permutations to be conducted for the bootstopping/bootstrap convergence test; minimum 100|
  --option_help               | -h                           | Display Help|

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
