### Making your own pipelines with Bash

This is an example of how to make your own pipeline with bash. This example uses four NGS samples for COVID-19 from NCBI and puts them through whole genome de novo assembly. 

SRA accession numbers:
[SRR11140744](https://www.ncbi.nlm.nih.gov/sra/SRR11140744)
[SRR11140746](https://www.ncbi.nlm.nih.gov/sra/SRR11140746)
[SRR11140748](https://www.ncbi.nlm.nih.gov/sra/SRR11140748)
[SRR11140750](https://www.ncbi.nlm.nih.gov/sra/SRR11140750)

<br/>

---

**Step 0 - Obtaining samples.**

We will assume that you have downloaded the reads from NCBI or know how to use fastq-dump. We will show the fastq-dump command here that we used to download the files. _**This part is not included in the pipeline.**_

```bash
for sra in SRR11140744 SRR11140746 SRR11140748 SRR11140750; do
    fastq-dump --outdir ${sra} --split-files --origfmt ${sra}
done
```
Output is:

```
Read 503344 spots for SRR11140744
Written 503344 spots for SRR11140744
Read 358971 spots for SRR11140746
Written 358971 spots for SRR11140746
Read 421395 spots for SRR11140748
Written 421395 spots for SRR11140748
Read 17657 spots for SRR11140750
Written 17657 spots for SRR11140750

```

You're starting directory should look like for this example, whether the reads were obtained manually or downloaded with fastq-dump:

```
.
├── SRR11140744
|   ├── SRR11140744_1.fastq
|   └── SRR11140744_2.fastq
├── SRR11140746
|   ├── SRR11140746_1.fastq
|   └── SRR11140746_2.fastq
├── SRR11140748
|   ├── SRR11140748_1.fastq
|   └── SRR11140748_2.fastq
└── SRR11140750
    ├── SRR11140750_1.fastq
    └── SRR11140750_2.fastq
```

We also know we will need a reference genome to help scaffold the contigs. We have downloaded the COVID19 reference genome [here](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) as a FASTA file.

<br/>

---

**Step 1 - Evaluate which modules you want to use.**

View all the module options using `haphpipe -h`. We have decided that we want to sample the reads (`sample_reads`), trim the reads (`trim_reads`), error correct the reads (`ec_reads`). We want to do genome de novo assembly, so we want to use `assemble_denovo` and `assemble_scaffold`. We then want to do refinement of the assembly (`refine_assembly`) and finalize the assembly (`finaliza_assembly`). 

<br/>

---

**Step 2 - Document files and options needed for each module.**

As we go through each module, we will take note of what we need and which input is specific for an individual sample (i.e., the input fastq reads will be different per sample, but we probably want the same number of reads for each sample.)

We know all modules have the option for a logfile. We'll include that in the bash scripting, but don't need to make notes of a logfile for each module.

_Part 2A - Sample reads._

Upon viewing the `haphpipe sample_reads -h` we see that we need to provide fastq reads 1 and 2, an output directory, and number of reads desired.

Sample specific options:

* fastq read 1
* fastq read 2
* output directory

Not sample specific options:

* number of reads desired 


_Part 2B - Trim reads._

Upon viewing the `haphpipe trim_reads -h` we see that we need to provide fastq reads 1 and 2, an output directory, number of reads, and number of CPUs desired. There is an option to change the timming commands, but we'll keep default. 

Sample specific options:

* fastq read 1 -- we will want these to be the subsampled read file
* fastq read 2 -- we will want these to be the subsampled read file
* output directory

Not sample specific options:

* number of reads desired 
* ncpu


_Part 2C - Error correct reads._

Upon viewing the `haphpipe ec_reads -h` we see that we need to provide fastq reads 1 and 2, an output directory, and number of CPUs desired.

Sample specific options:

* fastq read 1 -- we will want these to be the trimmed read file
* fastq read 2 -- we will want these to be the trimmed read file
* output directory

Not sample specific options:

* ncpu


_Part 2D - De novo assembly._

Upon viewing the `haphpipe assemble_denovo -h` we see that we need to provide fastq reads 1 and 2, an output directory, and number of CPUs desired.

Because we previously completed error correction. We will not want to include it here (i.e., we will invoke the option command `--no_error_correction`).

Sample specific options:

* fastq read 1 -- we will want these to be the error corrected reads
* fastq read 2 -- we will want these to be the error corrected reads
* output directory

Not sample specific options:

* ncpu
* no error correction option


_Part 2E - Assemble scaffold._

Upon viewing the `haphpipe assemble_scaffold -h` we see that we need to provide fasta file containing assembled contigs, an output directory, name to append to scaffold sequence and reference fasta desired. There is an option to change the timming commands, but we'll keep default. 

Sample specific options:

* assembled contigs file
* output directory
* sequence name

Not sample specific options:

* reference fasta


_Part 2F - Refine assembly._

Upon viewing the `haphpipe refine_assembly -h` we see that we need to provide fastq reads 1 and 2, an output directory, a reference sequence to refine, a sample ID and a maximum number of refinement steps, and number of CPUs desired.

Sample specific options:

* fastq read 1 -- we will want these to be the error corrected reads
* fastq read 2 -- we will want these to be the error corrected reads
* output directory
* reference fasta -- we will want this to be the sample's assembled scaffold
* Sample ID

Not sample specific options:

* ncpu
* maximum number of refinement steps (we'll do 3 for sake of simplicity and time)


_Part 2G - Finalize assembly._

Upon viewing the `haphpipe finalize_assembly -h` we see that we need to provide fastq reads 1 and 2, an output directory, a reference sequence to finalize, a sample ID and a maximum number of refinement steps, and number of CPUs desired. We could replace the preset bowtie2 option, but we will leave it as default for this sample pipeline.

Sample specific options:

* fastq read 1 -- we will want these to be the error corrected reads
* fastq read 2 -- we will want these to be the error corrected reads
* output directory
* reference fasta -- we will want this to be the sample's refined fasta sequence
* Sample ID

Not sample specific options:

* ncpu


_Part 2H - Gather the needed initial input for each sample._
We need to gather the necessary files that we will need to input for each sample so that we can make a script that takes in input files. We know from the User Guide that we can look at the output/input file types and names [here](https://gwcbi.github.io/haphpipe_docs/inout/). Because each module has a standard output file name, we will only need to be specific about the input for the raw fastq files for each sample. Therefore we need:

1. input raw fastq read 1
2. input raw fastq read 2
3. reference genome to scaffold against (this will be a covid19 reference) 
4. output directory for each sample
5. sample name

The other options we can code into the bash script for each module, since they will be the same for every sample in this analysis.

<br/>

---

**Step 3 - Create a bash script for each module.**

Now we will format a bash script for each module. 

_Part 3A - Input options._

Because we have a list of needed input options (specified by the user), we need to make a bash command to take in the inputs. 

First, we will specify the script name. We can do this one of two ways. i) explicitly or 2) through a command.

i) `SN='covid_genome'` # this sets the script name (SN variable) to covid_genome <br/>
ii) `SN=$(basename $0)` # this sets the script name (SN variable) to whatever the script filename is. If the script's file name is `covid.sh` then it is set as that. If the file name is `this_is_file` it will be set as that.

Because we are in charge of this script, we will explicitly set it.

--
Second, we want to set some input information for the user. 

We can do this as such:

```bash
read -r -d '' USAGE <<EOF
USAGE:
$SN [read1] [read2] [reference_fasta] [samp_id] <outdir>

----- COVID19 Genome Assembly Pipeline -----

This pipeline implements genome assembly using a denovo approach. Reads are
error-corrected and used to refine the scaffolded assembly, with up to 3
refinement steps. This pipeline is used as an example for advanced usage 
- making own pipeline in the User Guide. 

Input:
read1:             Fastq file for read 1. May be compressed (.gz)
read2:             Fastq file for read 2. May be compressed (.gz)
reference_fasta:   Reference sequence (fasta)
samp_id:           Sample ID
outdir:            Output directory (default is [sample_dir]/$SN)

EOF
```

--
Third, we want to provide help information for the script. We can do that using the avoce code, which has been saved in the variable `$USAGE`. Therefore,  

```bash
#--- Read command line args and if arg is -h, provide the usage information
[[ -n "$1" ]] && [[ "$1" == '-h' ]] && echo "$USAGE" && exit 0
```

If `pipeline.sh -h` is invoked here, then the output is:

```
USAGE:
covid_genome_assembly [read1] [read2] [reference_fasta] [samp_id] <outdir>

----- COVID19 Genome Assembly Pipeline -----

This pipeline implements genome assembly using a denovo approach. Reads are
error-corrected and used to refine the scaffolded assembly, with up to 3
refinement steps. This pipeline is used as an example for advanced usage
- making own pipeline in the User Guide.

Input:
read1:             Fastq file for read 1. May be compressed (.gz)
read2:             Fastq file for read 2. May be compressed (.gz)
reference_fasta:   Reference sequence (fasta)
samp_id:           Sample ID
outdir:            Output directory (default is [sample_dir]/covid_genome_assembly)
```

Note that this output looks identical to the information we put into `$USAGE` variable above in the second part of this section.

--
Fourth, we want to read in the command with the provided input options. Because we are using a bash script, the position of the input files is imparative.

```bash
#--- Read command line args
[[ -n "$1" ]] && raw1="$1"
[[ -n "$2" ]] && raw2="$2"
[[ -n "$3" ]] && refFA="$3"
[[ -n "$4" ]] && sampid="$4"
[[ -n "$5" ]] && outdir="$5"
```

--
Fifth, we want to check that the input files are provided and are not empty. We also want to set the outdirectory.

```bash
#--- Check that files are provided and exist
[[ -z ${raw1+x} ]] && echo "FAILED: read1 is not set" && echo "$USAGE" && exit 1
[[ ! -e "$raw1" ]] && echo "[---$SN---] ($(date)) FAILED: file $raw1 does not exist" && exit 1

[[ -z ${raw2+x} ]] && echo "FAILED: read2 is not set" && echo "$USAGE" && exit 1
[[ ! -e "$raw2" ]] && echo "[---$SN---] ($(date)) FAILED: file $raw2 does not exist" && exit 1

[[ -z ${refFA+x} ]] && echo "FAILED: refFA is not set" && echo "$USAGE" && exit 1
[[ ! -e "$refFA" ]] && echo "[---$SN---] ($(date)) FAILED: file $refFA does not exist" && exit 1

[[ -z ${sampid+x} ]] && echo "FAILED: sampid is not set" && echo "$USAGE" && exit 1

#--- Set outdirectory
[[ -z ${outdir+x} ]] && outdir=$(dirname $raw1)/$SN
mkdir -p $outdir
```

Sixth, we want to set the number of CPUs to use throughout the script.

```bash
#--- Determine CPUs to use
# First examines NCPU environment variable, then nproc, finally sets to  1
[[ -n "$NCPU" ]] && ncpu=$NCPU
[[ -z $ncpu ]] && ncpu=$(nproc 2> /dev/null)
[[ -z $ncpu ]] && ncpu=1
```

Seventh, we want to print out the variables and files. 

```bash
echo "[---$SN---] ($(date)) read1:             $raw1"
echo "[---$SN---] ($(date)) read2:             $raw2"
echo "[---$SN---] ($(date)) reference_fasta:   $refFA"
echo "[---$SN---] ($(date)) samp_id:           $sampid"
echo "[---$SN---] ($(date)) outdir:            $outdir"
echo "[---$SN---] ($(date)) num CPU:           $ncpu"
```

--
Finally, because here at GWU CBI, we like to time everything, we include a line of code to start a timer and end the timer:

```bash
#--- Start the timer
t1=$(date +"%s")

#### Put haphpipe module scripts here

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

```

--
After Part 3A, our pipeline file should read as such:

```bash
#!/usr/bin/env bash

###############################################################################
# This pipeline implements genome assembly using a denovo approach. Reads are
# error-corrected and used to refine the scaffolded assembly, with up to 3
# refinement steps. This pipeline is used as an example for advanced usage
# - making own pipeline in the User Guide.
###############################################################################
SN='covid_genome_assembly'

read -r -d '' USAGE <<EOF
USAGE:
$SN [read1] [read2] [reference_fasta] [samp_id] <outdir>

----- COVID19 Genome Assembly Pipeline -----

This pipeline implements genome assembly using a denovo approach. Reads are
error-corrected and used to refine the scaffolded assembly, with up to 3
refinement steps. This pipeline is used as an example for advanced usage
- making own pipeline in the User Guide.

Input:
read1:             Fastq file for read 1. May be compressed (.gz)
read2:             Fastq file for read 2. May be compressed (.gz)
reference_fasta:   Reference sequence (fasta)
samp_id:           Sample ID
outdir:            Output directory (default is [sample_dir]/$SN)

EOF

#--- Read command line args and if arg is -h, provide the usage information
[[ -n "$1" ]] && [[ "$1" == '-h' ]] && echo "$USAGE" && exit 0

#--- Read command line args
[[ -n "$1" ]] && raw1="$1"
[[ -n "$2" ]] && raw2="$2"
[[ -n "$3" ]] && refFA="$3"
[[ -n "$4" ]] && refGTF="$4"
[[ -n "$5" ]] && sampid="$5"
[[ -n "$6" ]] && outdir="$6"

#--- Check that files are provided and exist
[[ -z ${raw1+x} ]] && echo "FAILED: read1 is not set" && echo "$USAGE" && exit 1
[[ ! -e "$raw1" ]] && echo "[---$SN---] ($(date)) FAILED: file $raw1 does not exist" && exit 1

[[ -z ${raw2+x} ]] && echo "FAILED: read2 is not set" && echo "$USAGE" && exit 1
[[ ! -e "$raw2" ]] && echo "[---$SN---] ($(date)) FAILED: file $raw2 does not exist" && exit 1

[[ -z ${refFA+x} ]] && echo "FAILED: refFA is not set" && echo "$USAGE" && exit 1
[[ ! -e "$refFA" ]] && echo "[---$SN---] ($(date)) FAILED: file $refFA does not exist" && exit 1

[[ -z ${refGTF+x} ]] && echo "FAILED: refGTF is not set" && echo "$USAGE" && exit 1
[[ ! -e "$refGTF" ]] && echo "[---$SN---] ($(date)) FAILED: file $refGTF does not exist" && exit 1

[[ -z ${sampid+x} ]] && echo "FAILED: sampid is not set" && echo "$USAGE" && exit 1

#--- Set outdirectory
[[ -z ${outdir+x} ]] && outdir=$(dirname $raw1)/$SN
mkdir -p $outdir

#--- Determine CPUs to use
# First examines NCPU environment variable, then nproc, finally sets to  1
[[ -n "$NCPU" ]] && ncpu=$NCPU
[[ -z $ncpu ]] && ncpu=$(nproc 2> /dev/null)
[[ -z $ncpu ]] && ncpu=1


echo "[---$SN---] ($(date)) read1:             $raw1"
echo "[---$SN---] ($(date)) read2:             $raw2"
echo "[---$SN---] ($(date)) reference_fasta:   $refFA"
echo "[---$SN---] ($(date)) reference_gtf:     $refGTF"
echo "[---$SN---] ($(date)) samp_id:           $sampid"
echo "[---$SN---] ($(date)) outdir:            $outdir"
echo "[---$SN---] ($(date)) num CPU:           $ncpu"

#--- Start the timer
t1=$(date +"%s")

#### Put haphpipe module scripts here

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
```

_Part 3B - Sample reads._

Now we will begin constructing bash scripts for each module.
Each module will follow a similar trend. <br/>

1. We will set the stage name.
2. We will echo stage name to terminal. 
3. We will check to make sure the files that are the output of the stage are not already present. If they are present, we will skip the stage and continue on (no need to repeat the stage).
4. If the files are not present, we will complete the stage and call the command.
5. If the command completes, we will print to terminal. If command fails, we will print that to terminal and quit the script. <br/>

Remember, you can find the output file names [here](https://gwcbi.github.io/haphpipe_docs/inout/). The ouput names for this stage are: `sample_1.fastq` and `sample_2.fastq`. We have also decided to subsample the number of reads to 50,000 reads for each sample. <br/>

Now this in the input for the base haphpipe command for this stage:

```
haphpipe sample_reads\
 --fq1 $raw1\
 --fq2 $raw2\
 --nreads 50000\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
```

The entire stage's bash script is here:

```bash
stage="sample_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

# if the sampled files are present, skip this stage. Otherwise, call sample_reads
if [[ -e $outdir/sample_1.fastq && -e ${outdir}/sample_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage sample_1.fastq, sample_2.fastq"
else
	# this reads in the sample_reads command and saves it in the variable cmd
    read -r -d '' cmd <<EOF
haphpipe sample_reads\
 --fq1 $raw1\
 --fq2 $raw2\
 --nreads 50000\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

```

_Part 3C - Trim reads._

Remember, you can find the output file names [here](https://gwcbi.github.io/haphpipe_docs/inout/). The ouput names for this stage are: `trimmed_1.fastq` and `trimmed_2.fastq`. <br/>

Now, because we are doing trimming after the sampled reads module, we need to change the input fastq reads for this module to be the fastq reads output from the previous step (sample reads). Therefore, the input fastq files are named `sample_1.fastq` and `sample_2.fastq`. Also remember that these files are now contained in the outdirectory specified by the input, so we have to list the path to the input fastq files. <br/>

Now this in the input for the base haphpipe command for this stage:

```
haphpipe trim_reads\
 --ncpu $ncpu\
 --fq1 ${outdir}/sample_1.fastq\
 --fq2 ${outdir}/sample_2.fastq\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
```

The entire stage's bash script is here:

```bash
stage="trim_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/trimmed_1.fastq && -e ${outdir}/trimmed_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage trimmed_1.fastq,trimmed_2.fastq"
else
    read -r -d '' cmd <<EOF
haphpipe trim_reads\
 --ncpu $ncpu\
 --fq1 ${outdir}/sample_1.fastq\
 --fq2 ${outdir}/sample_2.fastq\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi
```


_Part 3D - Error correct reads._

Remember, you can find the output file names [here](https://gwcbi.github.io/haphpipe_docs/inout/). The ouput names for this stage are: `corrected_1.fastq` and `corrected_2.fastq`. <br/>

Now, because we are doing error correction after the trimming module, we need to change the input fastq reads for this module to be the fastq reads output from the previous step (trimmed reads). Therefore, the input fastq files are named `trimmed_1.fastq` and `trimmed_2.fastq`. Also remember that these files are now contained in the outdirectory specified by the input (just like the sampled reads in the previous step), so we have to list the path to the input fastq files. <br/>

Now this in the input for the base haphpipe command for this stage:

```
haphpipe ec_reads\
 --ncpu $ncpu\
 --fq1 ${outdir}/trimmed_1.fastq\
 --fq2 ${outdir}/trimmed_2.fastq\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
```

The entire stage's bash script is here:

```bash
stage="ec_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/corrected_1.fastq && -e $outdir/corrected_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage corrected_1.fastq,corrected_2.fastq"
else
    read -r -d '' cmd <<EOF
haphpipe ec_reads\
 --ncpu $ncpu\
 --fq1 ${outdir}/trimmed_1.fastq\
 --fq2 ${outdir}/trimmed_2.fastq\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi
```


_Part 3E - De novo assembly._

Remember, you can find the output file names [here](https://gwcbi.github.io/haphpipe_docs/inout/). The ouput name for this stage is `denovo_contigs.fna`. <br/>

Now, because we are doing denovo assembly after the error correction module, we need to change the input fastq reads for this module to be the fastq reads output from the previous step (error corrected reads). Therefore, the input fastq files are named `corrected_1.fastq` and `corrected_2.fastq`. Also remember that these files are now contained in the outdirectory specified by the input, so we have to list the path to the input fastq files. Finally, remember we want to specify that we do NOT want to do another round of error correction.  <br/>

Now this in the input for the base haphpipe command for this stage:

```
haphpipe assemble_denovo\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --no_error_correction\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
```

The entire stage's bash script is here:

```bash
stage="assemble_denovo"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/denovo_contigs.fna ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage denovo_contigs.fna"
else
    read -r -d '' cmd <<EOF
haphpipe assemble_denovo\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --no_error_correction\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi
```


_Part 3F - Assemble scaffold._

Remember, you can find the output file names [here](https://gwcbi.github.io/haphpipe_docs/inout/). The ouput name for this stage is `scaffold_assembly.fa`. There are other outputs, but we focus on this one for further use in the refinement and finalize modules. <br/>

Now, because we are doing scaffold assembly after the denovo assembly module, we need to specify that the input file is the output contig file (`denovo_contigs.fna`). Again, remember that this file is now contained in the outdirectory specified by the input, so we have to list the path too. We also have to use our input reference fasta file and input sampleID from the script command. <br/>

Now this in the input for the base haphpipe command for this stage:

```
haphpipe assemble_scaffold\
 --contigs_fa ${outdir}/denovo_contigs.fna\
 --ref_fa ${refFA}\
 --seqname ${sampid}\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
```

The entire stage's bash script is here:

```bash
stage="assemble_scaffold"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/scaffold_assembly.fa ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage scaffold_assembly.fa"
else
    read -r -d '' cmd <<EOF
haphpipe assemble_scaffold\
 --contigs_fa ${outdir}/denovo_contigs.fna\
 --ref_fa ${refFA}\
 --seqname ${sampid}\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi
```


_Part 3G - Refine assembly._

Remember, you can find the output file names [here](https://gwcbi.github.io/haphpipe_docs/inout/). The ouput name for this stage is `refined.fna`. <br/>

Now, because we are doing refining the assembly after the scaffold assembly module, we need to change the input fastq reads for this module to be the fastq reads output from the error correction step (error corrected reads). Therefore, the input fastq files are named `corrected_1.fastq` and `corrected_2.fastq`. Also remember that these files are now contained in the outdirectory specified by the input, so we have to list the path to the input fastq files. Finally, we need to specify that the input reference file is the `scaffold_assembly.fa` file from the previous scaffold assembly step (above), which is also located in the outdirectory.  <br/>

Now this in the input for the base haphpipe command for this stage:

```
hp_refine_assembly\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --ref_fa ${outdir}/scaffold_assembly.fa\
 --sample_id ${sampid}\
 --max_step 5\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
```

The entire stage's bash script is here:

```bash
stage="refine_assembly"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e ${outdir}/refined.fna ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage refined.fna"
else
    read -r -d '' cmd <<EOF
hp_refine_assembly\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --ref_fa ${outdir}/scaffold_assembly.fa\
 --sample_id ${sampid}\
 --max_step 5\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi
```

_Part 3H - Finalize assembly._

Remember, you can find the output file names [here](https://gwcbi.github.io/haphpipe_docs/inout/). The ouput files for this stage are: `final.fna`, `final.bam`, `final.vcf.gz`. <br/>

Now, because we are doing finalizing the assembly after the refinement module, we need to change the input fastq reads for this module to be the fastq reads output from the error correction step (error corrected reads). Therefore, the input fastq files are named `corrected_1.fastq` and `corrected_2.fastq`. Also remember that these files are now contained in the outdirectory specified by the input, so we have to list the path to the input fastq files. Finally, we need to specify that the input reference file is the `refined.fna` file from the previous refinement step (above), which is also located in the outdirectory.  <br/>

Now this in the input for the base haphpipe command for this stage:

```
hp_finalize_assembly\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --sample_id ${sampid}\
 --ref_fa ${outdir}/refined.fna\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
```

The entire stage's bash script is here:

```bash
stage="finalize_assembly"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e ${outdir}/final.fna && -e ${outdir}/final.bam && -e ${outdir}/final.vcf.gz ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage final.fna,final.bam,final.vcf.gz"
else
    read -r -d '' cmd <<EOF
hp_finalize_assembly\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --sample_id ${sampid}\
 --ref_fa ${outdir}/refined.fna\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi
```

_Part 3I - Gather all the individual module scripts into the final pipeline script._

For readable code, we will separate each module with `###` like so:

```bash
###############################################################################
# Step #: description here
###############################################################################

insert code here ..
```

Once we concatenate all our code, we end up with this:

```bash
###############################################################################
# Step 1: Sample Reads
###############################################################################
stage="sample_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

# if the sampled files are present, skip this stage. Otherwise, call sample_reads
if [[ -e $outdir/sample_1.fastq && -e ${outdir}/sample_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage sample_1.fastq, sample_2.fastq"
else
	# this reads in the sample_reads command and saves it in the variable cmd
    read -r -d '' cmd <<EOF
haphpipe sample_reads\
 --fq1 $raw1\
 --fq2 $raw2\
 --nreads 50000\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 2: Trim Reads
###############################################################################
stage="trim_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/trimmed_1.fastq && -e ${outdir}/trimmed_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage trimmed_1.fastq,trimmed_2.fastq"
else
    read -r -d '' cmd <<EOF
haphpipe trim_reads\
 --ncpu $ncpu\
 --fq1 $raw1\
 --fq2 $raw2\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 3: Error correction using Spades
###############################################################################
stage="ec_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/corrected_1.fastq && -e $outdir/corrected_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage corrected_1.fastq,corrected_2.fastq"
else
    read -r -d '' cmd <<EOF
haphpipe ec_reads\
 --ncpu $ncpu\
 --fq1 ${outdir}/trimmed_1.fastq\
 --fq2 ${outdir}/trimmed_2.fastq\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 4: Denovo assembly
###############################################################################
stage="assemble_denovo"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/denovo_contigs.fna ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage denovo_contigs.fna"
else
    read -r -d '' cmd <<EOF
haphpipe assemble_denovo\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --no_error_correction\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 5: Scaffold assembly
###############################################################################
stage="assemble_scaffold"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/scaffold_assembly.fa ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage scaffold_assembly.fa"
else
    read -r -d '' cmd <<EOF
haphpipe assemble_scaffold\
 --contigs_fa ${outdir}/denovo_contigs.fna\
 --ref_fa ${refFA}\
 --seqname ${sampid}\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi


###############################################################################
# Step 5: Refine assembly
###############################################################################
stage="refine_assembly"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e ${outdir}/refined.fna ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage refined.fna"
else
    read -r -d '' cmd <<EOF
hp_refine_assembly\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --ref_fa ${outdir}/scaffold_assembly.fa\
 --sample_id ${sampid}\
 --max_step 5\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 5: Finalize assembly
###############################################################################
stage="finalize_assembly"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e ${outdir}/final.fna && -e ${outdir}/final.bam && -e ${outdir}/final.vcf.gz ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage final.fna,final.bam,final.vcf.gz"
else
    read -r -d '' cmd <<EOF
hp_finalize_assembly\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --sample_id ${sampid}\
 --ref_fa ${outdir}/refined.fna\
 ${quiet} --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi
```


<br/>

---

**Step 4 - Combine both the input code and bash scripts for each module into a single script.**

We named this script `covid_genome_assembly.sh` and the entire code is:

```bash
#!/usr/bin/env bash

###############################################################################
# This pipeline implements genome assembly using a denovo approach. Reads are
# error-corrected and used to refine the scaffolded assembly, with up to 3
# refinement steps. This pipeline is used as an example for advanced usage
# - making own pipeline in the User Guide.
###############################################################################
SN='covid_genome_assembly'

read -r -d '' USAGE <<EOF
USAGE:
$SN [read1] [read2] [reference_fasta] [samp_id] <outdir>

----- COVID19 Genome Assembly Pipeline -----

This pipeline implements genome assembly using a denovo approach. Reads are
error-corrected and used to refine the scaffolded assembly, with up to 3
refinement steps. This pipeline is used as an example for advanced usage
- making own pipeline in the User Guide.

Input:
read1:             Fastq file for read 1. May be compressed (.gz)
read2:             Fastq file for read 2. May be compressed (.gz)
reference_fasta:   Reference sequence (fasta)
samp_id:           Sample ID
outdir:            Output directory (default is [sample_dir]/$SN)

EOF

#--- Read command line args and if arg is -h, provide the usage information
[[ -n "$1" ]] && [[ "$1" == '-h' ]] && echo "$USAGE" && exit 0

#--- Read command line args
[[ -n "$1" ]] && raw1="$1"
[[ -n "$2" ]] && raw2="$2"
[[ -n "$3" ]] && refFA="$3"
[[ -n "$4" ]] && sampid="$4"
[[ -n "$5" ]] && outdir="$5"

#--- Check that files are provided and exist
[[ -z ${raw1+x} ]] && echo "FAILED: read1 is not set" && echo "$USAGE" && exit 1
[[ ! -e "$raw1" ]] && echo "[---$SN---] ($(date)) FAILED: file $raw1 does not exist" && exit 1

[[ -z ${raw2+x} ]] && echo "FAILED: read2 is not set" && echo "$USAGE" && exit 1
[[ ! -e "$raw2" ]] && echo "[---$SN---] ($(date)) FAILED: file $raw2 does not exist" && exit 1

[[ -z ${refFA+x} ]] && echo "FAILED: refFA is not set" && echo "$USAGE" && exit 1
[[ ! -e "$refFA" ]] && echo "[---$SN---] ($(date)) FAILED: file $refFA does not exist" && exit 1

[[ -z ${sampid+x} ]] && echo "FAILED: sampid is not set" && echo "$USAGE" && exit 1

#--- Set outdirectory
[[ -z ${outdir+x} ]] && outdir=$(dirname $raw1)/$SN
mkdir -p $outdir

#--- Determine CPUs to use
# First examines NCPU environment variable, then nproc, finally sets to  1
[[ -n "$NCPU" ]] && ncpu=$NCPU
[[ -z $ncpu ]] && ncpu=$(nproc 2> /dev/null)
[[ -z $ncpu ]] && ncpu=1


echo "[---$SN---] ($(date)) read1:             $raw1"
echo "[---$SN---] ($(date)) read2:             $raw2"
echo "[---$SN---] ($(date)) reference_fasta:   $refFA"
echo "[---$SN---] ($(date)) samp_id:           $sampid"
echo "[---$SN---] ($(date)) outdir:            $outdir"
echo "[---$SN---] ($(date)) num CPU:           $ncpu"

#--- Start the timer
t1=$(date +"%s")

###############################################################################
# Step 1: Sample Reads
###############################################################################
stage="sample_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

# if the sampled files are present, skip this stage. Otherwise, call sample_reads
if [[ -e $outdir/sample_1.fastq && -e ${outdir}/sample_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage sample_1.fastq, sample_2.fastq"
else
	# this reads in the sample_reads command and saves it in the variable cmd
    read -r -d '' cmd <<EOF
haphpipe sample_reads\
 --fq1 $raw1\
 --fq2 $raw2\
 --nreads 50000\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 2: Trim Reads
###############################################################################
stage="trim_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/trimmed_1.fastq && -e ${outdir}/trimmed_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage trimmed_1.fastq,trimmed_2.fastq"
else
    read -r -d '' cmd <<EOF
haphpipe trim_reads\
 --ncpu $ncpu\
 --fq1 $raw1\
 --fq2 $raw2\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 3: Error correction using Spades
###############################################################################
stage="ec_reads"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/corrected_1.fastq && -e $outdir/corrected_2.fastq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage corrected_1.fastq,corrected_2.fastq"
else
    read -r -d '' cmd <<EOF
haphpipe ec_reads\
 --ncpu $ncpu\
 --fq1 ${outdir}/trimmed_1.fastq\
 --fq2 ${outdir}/trimmed_2.fastq\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 4: Denovo assembly
###############################################################################
stage="assemble_denovo"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/denovo_contigs.fna ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage denovo_contigs.fna"
else
    read -r -d '' cmd <<EOF
haphpipe assemble_denovo\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --no_error_correction\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 5: Scaffold assembly
###############################################################################
stage="assemble_scaffold"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e $outdir/scaffold_assembly.fa ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage scaffold_assembly.fa"
else
    read -r -d '' cmd <<EOF
haphpipe assemble_scaffold\
 --contigs_fa ${outdir}/denovo_contigs.fna\
 --ref_fa ${refFA}\
 --seqname ${sampid}\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi


###############################################################################
# Step 5: Refine assembly
###############################################################################
stage="refine_assembly"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e ${outdir}/refined.fna ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage refined.fna"
else
    read -r -d '' cmd <<EOF
hp_refine_assembly\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --ref_fa ${outdir}/scaffold_assembly.fa\
 --sample_id ${sampid}\
 --max_step 5\
 --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

###############################################################################
# Step 5: Finalize assembly
###############################################################################
stage="finalize_assembly"
echo -e "\n[---$SN---] ($(date)) Stage: $stage"

if [[ -e ${outdir}/final.fna && -e ${outdir}/final.bam && -e ${outdir}/final.vcf.gz ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage final.fna,final.bam,final.vcf.gz"
else
    read -r -d '' cmd <<EOF
hp_finalize_assembly\
 --ncpu $ncpu\
 --fq1 ${outdir}/corrected_1.fastq\
 --fq2 ${outdir}/corrected_2.fastq\
 --sample_id ${sampid}\
 --ref_fa ${outdir}/refined.fna\
 ${quiet} --logfile ${outdir}/haphpipe.out\
 --outdir ${outdir}
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi


#### Put haphpipe module scripts here

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
```

---

**Step 5 - Executing the script.**

Example to run:

`bash covid_genome_assembly.sh SRR11140750/SRR11140750_1.fastq SRR11140750/SRR11140750_2.fastq SARSCoV2.NC_045512.COVID19.fasta SRR11140750`

Example to run in a loop over all the samples:

```bash
for sra in SRR11140744 SRR11140746 SRR11140748 SRR11140750; do
    bash covid_genome_assembly.sh ${sra}/${sra}_1.fastq ${sra}/${sra}_2.fastq SARSCoV2.NC_045512.COVID19.fasta ${sra} covid_genome_assembly
done
```

Directories should look like such after running this script:


```
.
├── SRR11140744
|   ├── SRR11140744_1.fastq
|   ├── SRR11140744_2.fastq
|   └── covid_genome_assembly
|       ├── corrected_1.fastq
|       ├── corrected_2.fastq
|       ├── corrected_U.fastq
|       ├── denovo_contigs.fna
|       ├── denovo_summary.txt
|       ├── haphpipe.out
|       ├── final.bam
|       ├── final.bam.bai
|       ├── final_bt2.out
|       ├── final.fna
|       ├── final.vcf.gz
|       ├── final.vcf.gz.tbi
|       ├── refined.01.fna
|       ├── refined_bt2.01.out
|       ├── refined.fna
|       ├── refined_bt2.out
|       ├── refined_summary.out
|       ├── sample_1.fastq
|       ├── sample_2.fastq
|       ├── scaffold_aligned.fa
|       ├── scaffold_assembly.fa
|       ├── scaffold_imputed.fa
|       ├── scaffold_padded.out
|       ├── trimmed_1.fastq
|       ├── trimmed_2.fastq
|       ├── trimmed_U.fastq
|       └── trimmomatic_summary.out
├── SRR11140746
|   ├── SRR11140746_1.fastq
|   ├── SRR11140746_2.fastq
|   └── covid_genome_assembly
....
```


