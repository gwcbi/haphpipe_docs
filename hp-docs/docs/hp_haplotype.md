**hp_haplotype** includes haplotype assembly stages. HAPHPIPE implements PredictHaplo ([paper](https://www.ncbi.nlm.nih.gov/pubmed/26355517)), although other haplotype reconstruction programs can be utilized outside of HAPHPIPE using the final output of HAPHPIPE, typically with the final consensus sequence (FASTA) file, reads (raw, trimmed, and/or corrected), and/or final alignment (BAM) file as input.
Use -h after any command for a list of options. 

_Note:_ PredictHaplo is not compatable with read files including read IDs (.1 and .2 appended at the end of read names for read 1 and read 2, respectively). If your file has read IDs, use the following commands to create new files with the read IDs taken out before running PredictHaplo:
```
	cat corrected_1.fastq | sed 's/\.1 / /' > corrected_1_fixed.fastq
	cat corrected_2.fastq | sed 's/\.1 / /' > corrected_2_fixed.fastq
```

### *predict_haplo*
Haplotype identification with PredictHaplo. Input is reads in FASTQ format and and reference sequence in FASTA format. Output is the longest global haplotype file and corresponding HTML file. _Note: PredictHaplo must be installed separately before running this stage._ 

**Usage:**

`haphpipe predict_haplo [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --ref_fa <FASTA> --interval_txt [TXT] [--outdir]`

**(or):**

`hp_predict_haplo [OPTIONS] [SETTINGS] --fq1 <FASTQ> --fq2 <FASTQ> --ref_fa <FASTA> --interval_txt [TXT] [--outdir]`

*Output files:* <br> *best.fa*

*Input/Output Arguments:* 

Option         | Description
---------------|-------------
--fq1          | Fastq file with read 1.
--fq2          | Fastq file with read 2.
--ref_fa       | Reference sequence used to align reads (Fasta).
--interval_txt | File with intervals to perform haplotype reconstruction.
--outdir       | Output directory (default: current directory).

*Options:*

Option           | Description
-----------------|-------------
--min_readlength | Minimum read length passed to PredictHaplo (default: 36).

*Settings:*

Option      | Description
------------|-------------
--keep_tmp  | Keep temporary directory (default: False).
--quiet     | Do not write output to console (silence stdout and stderr) (default: False).
--logfile   | Append console output to this file.
--debug     | Print commands but do not run (default: False).


_Example usage:_

```
haphpipe predict_haplo corrected_1.fastq --fq2 corrected_2.fastq --ref_fa final.fna
```

### *ph_parser*
Returns PredictHaplo output as a correctly formatted FASTA file. Input is the output file from predict_haplo (longest global .fas file). Output is a correctly formatted FASTA file.

**Usage:**

`haphpipe ph_parser [OPTIONS] [SETTINGS] --haplotypes_fa <best.fas> [--outdir]`

**(or):**

`hp_ph_parser [OPTIONS] [SETTINGS] --haplotypes_fa <best.fas> [--outdir]`


*Output files:* <br> *ph_summary.txt* <br> *ph_haplotypes.fna*

*Input/Output Arguments:* 

Option          | Description
----------------|-------------
--haplotypes_fa | Haplotype file created by PredictHaplo.
--outdir        | Output directory (default: current directory).

*Options:*

Option      | Description
------------|-------------
--prefix    | Prefix to add to sequence names.
--keep_gaps | Do not remove gaps from alignment (default: False).

*Settings:*

Option    | Description
----------|-------------
--quiet   | Do not write output to console (silence stdout and stderr) (default: False).
--logfile |Append console output to this file.



_Example usage:_
```
haphpipe ph_parser PH01.best_1_864.fas
```

Example:

By default, PredictHaplo outputs their own unique version of a fasta file. It includes the frequency, some information regarding their unqiue overlapping scores, and their unique confidence scores. This file is always named `PH#.best_#_#.fas`, where the first number is the reconstructed haplotype number and the next numbers are the start and end of the longest haplotype reconstructed by PredictHaplo.

<details>
  <summary>PredictHaplo Output</summary>
```
	$ cat PH01.best_1_1208.fas
	>reconstructed_0
	;Freq:0.190638
	;Overlap quality scores (5,10):1,1
	;Confidence scores
	;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~n7~~~y~~~~~t~~~i~jjk~~zz
	;~~~~{~~~~~N{~{~Sx~~~~~~z~K~F~~~~~~~y~~~~~~~|~~F~wx|~~{~~|~~s<|]~~~~~~
	;~m~kj|{|v~{|_`~~~z~~~~~~~jy{y~~~~~a~__~~|~~~~~{wXZ~}~~~qm~xV~~~~~~~}~
	;Q~}~y||~~~}}~~~z~~~~~~{}A~}b|~~u~~|}}|~}}z~}bx~~n||~~||{~}~~d}~bz~~~}
	;|~}}~~~}~~~{|}}g{~~~}~r~}}~~~~u~~{kx{~}~}~|~}~}{}~}~}||~~~~[~}~}}~~~~
	;~~}~~~U||U}~|}}~~}}~~~~}u~b|}~w~~~~~{}|wv~}Dxzp{|}~@~~P~}~}~V~~z}~|ry
	;q~}|}~}~~}t~o~}~f}~~}{~~~}~{~~~~~~~}~~}~~|~~~M~~}~}x~}c~}^v~~~yzA~}~}
	;y}~z}~~~~~~~~~{z~~}~~~}~{}~~~~~}~~~~~~~|~~v~~}~~|y~]|{~||~~~~~~~~||~~
	;Y||~~|Q~|~~~|~~~~~|~~~z~~z~{{{y~~~~~~~~~~w{{~wz|~Z|~z|~~}p|~~|}}~~x}}
	;z}~}}|a|}}}}{}|~~}}~}}~}~{~|~}}~}{}|}|}}}~|}}}}{}}}|}}|}|}}}|}}}}~}}}
	;||}}}}{}}~}~}}}}}}}}}}}~}}}}}}}}}}}}}}}}}y}~}}}}}}}}}}~|}}}}|}}}}|}}}
	;}}}}}}}}|}}}}}}}}~}}}}}}|}}}}}}}}|}}}}}|}}}}|~|}}}}{}}}}}}}}}}}}}}}}}
	;}}=}}}{}}}}}}}}}}}}}}}}}}}}}}}}~||}}|~}}}}}{}}}}}|}}|}}}}}}}}}}}}}|}}
	;|}~}}}}|}}}}}}}}|}|~J}}}}}}~}}}}}}}}}}}}}}}}}}}}|}}`~}}}}|}}}}}}}}}}~
	;|}}}}}}}}}}}~}}}|}}}}}}}}}}|}}}}}}}}}}}t}}{}`}}}}}}|}~}~}|~}}}|k}}}}}
	;|}|~}}|}}~}}~}|}}z}}}}}}}}}}}}}}~}}~}}}~~}|}}}~}}}}}}}}~}~}}||~~}~z}}
	;~~}~}||~|~}|{}||~|z~~||}~|~}}|}~~}|}}~}}|}z~~~~{~}{}~y~~~~~{|}}~~|y~~
	;~|~~||~~~~~~~|n~~{~~~~~~~~~~~~~~~
	;EndOfComments
	ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAATATCTGGCCTTCCCGCAAGGGAAG
	GCCAGGGAACTTTCCTCAGAGCAGGCTAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGCTTGGGGAA
	GCAACAACAAAGCCCTCTCAGAAACAGGAGCCAATAGACAAGGAAATGTATCCTTTAACTTCCCTCAGAT
	CACTCTTTGGCAACGACCCCTTGTCACAATAAAGATAGGGGAGCAACTAAAGGAAGCCCTGTTAGATACA
	GGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAGAATGATAGGGGGAA
	TCGGGGGTTTTATCAAAGTAAGACAGTATGATCAGATAGCCATAGACATCTGTGGCCATAAAGCTATAGG
	TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGCTTGGTTGCACT
	TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTGAAGCCAGGAATGGATGGCCCAAAGG
	TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAAGA
	AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGGAAAAAAGAC
	AGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAGGACTTCTGGGAAGTCC
	AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGCGGGTGATGC
	ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTATACTGCATTTACCATACCTAGTGTAAACAAT
	GAGACACTAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
	AAGCTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACAT
	GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
	CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
	GTTATGAACTCCATCCTT
	>reconstructed_1
	;Freq:0.294104
	;Overlap quality scores (5,10):1,1
	;Confidence scores
	;~~~~~~~~z~~~~~~~~~|~y~~~~u}|~~~~~}}~~~~~~~ya~|~~}}~~~y~~~j~YXH~~s~
	;~~~~z~~|~~b|}^}xZ}}~}z~t~r~{}~}x}}yb}~~~}~}u~~}~gu|}}{~|}}~ozsw}|}}h~
	;|l}v|^][t}zz}vw}}s}}}}~v}|~|~~}}}~~}}}}}~}}}}}}zo\}}}~}vn}iv}}}}}}}}}
	;v}}~|u}}}}}}}}}|}}}}}}|}^}}[}~}r}}}{|}}|t|}}{z}}pz}}{}}}}}}}v~}y|~}}}
	;}}}}~}}~}}}tz}}`}}}~}}q}}}}}~}y}}|r||}}}}~}}~}}{}}}}}}}}}}}y~}}}}t}}}
	;}}}}}}{|}|}}{r}}}}tv~}}}t}m|~}h}}}}}v}}qt}|y|{|u~}}|}~_}~|~~|}~|}~}|e
	;Y}~}}~}}~}v}t}}}j}}~}}}~}}~|}}}}|}}|}}}~}}}}};}}~}}x}}|}}yu}}}wy`}|~~
	;{}}x}}}}}~}}}~{v}}}}}{}~}}~}~}}|}}}}}}}}}}u}}y~~{g}L}}}}}|~}{~|}|}|}}
	;l}}}}{_}}}}~z~}|}~~|~~x|~}~z~~|~}~~~}~~~}y}~~|k}}Pf}w~~~~T}}}F~~}~z}{
	;{~|}|zZ~s~|}r~}}~~|}}}}}~|}}}}n}}}~}}}}~}}}}}}}}}}}|~}}}}}}~~}}~~~~~}
	;N}}}}}|~}}~}}}}{}}}~}}|~~|}|}}}~|}~}}}}}}D}~}}}}}}}}}}}}}}}}|}}z}|}}}
	;}k}}}}}}}}}}}}}}}}}}}|}}}}}|}}}}}|}}}}}|}}}}|~{~}}}G}~}}}}}~}}}}}}}}}
	;}}h}}~|}|}}}}}}}}}}}}}}}}}}}}}}}}{}}|}}}}|~|}}}}}}}M}}}}}}}~}}}}}}}}}
	;}}~|}}}}}~~~}}}~}}|}g}}}}}}}}~|}}}}}}}}}}}}}}}}}|}~t~}|}}}}}}}||~}}}}
	;|}}}}~}}}}}}~}}}gi}}~}}}}}}}}}}}}|}}~}}g}~}}v}}}}}|{}}}}}}~}}~|e}}}}}
	;}}|}}~}~}~}}~}|~}|}~|||}}}}|~}}~}}}}}}}}|}}~}}~}|~}|y~}~}}}}}|}}~~}~}
	;}}~}~}|~}}}~}}}}~}{}}~}}~}}}||~~}||~}x|}{~|}|~~|~}}}|{~}}~}{}}||~}}}~
	;|~}~}}}~|}}~~|,~~~~~~~~~~~~~~m~~~
	;EndOfComments
	ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCCGCGGCGGAAG
	GCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAG
	GAGACAACAACTCCCTCTCAGAAGCAGGAGCCGAAGGACAAGGAAACATATCCTTTAGCTTCCCTCAAAT
	CACTCTTTGGCAACGACCCCTTGTTACAATAAAGATAGGGGGGCAGCTAAAGGAAGCTCTATTAGATACA
	GGAGCAGATGACACAGTATTAGAAGAAATGGATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAA
	TTGGAGGTTTTATCAAAGTAAAACAGTATGATCAGATACCCATAGAAATTTGTGGACATAAAGTGATAGG
	TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGCTTGGTTGCACT
	TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTGAAGCCAGGAATGGATGGCCCAAAGG
	TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAATAGAAATCTGTGCAGAAATGGAAAAGGA
	AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGAAAAAAAGAC
	AGTACTAAATGGAGAAAATTAGTAGATTTCAAAGAACTTAATAAGAGAACTCAGGACTTCTGGGAAGTCC
	AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAGTCAGTAACAGTACTGGATGTGGGTGATGC
	ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTACACTGCATTTACCATACCTAGTGTAAACAAT
	GAGACACCAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
	AATGTAGCATGACAAAAATCTTAGAACCTTTTAGAAAACAAAATCCAGATATAGTTATCTATCAATACAT
	GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
	CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
	GTTATGAACTCCATCCTT
	>reconstructed_2
	;Freq:0.515259
	;Overlap quality scores (5,10):4,4
	;Confidence scores
	;~~~~~~~w~~zz~~~~z~~w~|wy~|~~{|~~|~~~~~~~~}~{~hD~~~}|~|}}t~~_;~tue}}S}
	;~}~~p}}l}~P}}ZlvL~z}}h~g~A~x}}}s}}{S~}~}}}}l}~|}SOq}}q}w}}~i\pz~w~~z~
	;}7~|vzUu4|smtgm}~W}}~~}r}n`f`}{PZLSI:Vs_V_}}}}_<db}}}}}CQ}YJ}t}}}~}}}
	;T}}~|}}}}}}}~}}k}}}}}}q}L}}]}}}L}}}tt|~}{q}}{M}lAIu}}}s{}}}~y~}st}|}}
	;z}}}}|}}}}}{q|}Z}}}}}}P~}|}|}}E}}cLbi}vu}}|}}|}g}}}}}}|}~}~{}~}|}{}}}
	;}~}|}}}{{|~}pp}~~~{|}}~}[yEx}}\}}}{}N}uoK}}pauSy}}}v~~J}}q}}v~}c}~}\i
	;]|}}~}|}~}r}Y}}}B}}~||~}}}}{}}}|}|~|}~~~~}}~{w}{}~|R|}s|}k{~}}}ra}q}}
	;o}}d}~}}}}}~}}|X~}}~}e}}}}}}}|}y|}}}}~}|~}c~~b~~nL~a}}}}}~}}m}}~{}}}}
	;g}}}}uLm~l}}`~|n}}~}}|d~~y~u}|tz~~~~|}}~~t}}}tx}|Tu~v}}~}^|~~w~|~}m{}
	;u}}}|w;}{}{}r}||}}z|}|~|}z}||}y}}}~}}}|}}~}}}|}y}}}}}}}~}||}}}}}}~}}|
	;v{}|}}{}}}|~|}}|}}}}}}}}}}|}}}}}|}~}}}}}}p}}||}}}}}~}}}{}}}}}}}z}}}}}
	;}|}}~}}}|}}}}}}}}}}}}}}}}}}}}}}|}||}}|}|}}}}|}|}}}}u}}}}}}}}}}}}}}}}}
	;}}g}}}{}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}y}}}}}|}m}}}}}}}}}}}}}}}}}
	;}}}}}}}}}}}}}}}}}}|}i}}}}}}}}}}|}}}}}}}}}}}}}}}}}}}H}}}}}}}}}|}{}}}}}
	;|}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}y}}|}H}}}}}}}}}}}}}}}}}}u}}}}}
	;}}}}}}}}}}}}}}}}}z}}}}}}}}}}}}}}}}|}}}|}}}}}}}}}|}}}}}}}}}}}}|}}}}|}}
	;}}}}}}}}}}}}||}|}}|}}}}}}|}}|}}}~}|}|}}}|}{}}}}|}}}}|{|}}}}x||}}}|{}}
	;}}}}}}}~}}~}~|n}~}}}}}~}}~~|~~~y~,
	;EndOfComments
	ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAAATCTGGCCTTCCCACAAGGGAAG
	GCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAA
	GAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGAAATGTATCCTTTAGCTTCCCTCAAAT
	CACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACA
	GGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAA
	TTGGAGGTTTTATCAAAGTAAAACAGTATGATCAGATACCCATAGAAATCTGTGGACATAAAGCTATAGG
	TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACT
	TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAG
	TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGA
	AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGGAAAAAAGAC
	AGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCCGGGAAGTCC
	AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGC
	ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTATACTGCATTTACCATACCTAGTGTAAACAAT
	GAGACACCAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
	AAGCTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACAT
	GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
	CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
	GTTATGAACTCCATCCCC
```
</details>

*ph_parser* takes this output and creates a proper fasta file with each resontructed haplotype and a text file that has the hpalotype diversity estimate.

<details>
  <summary>PH Parser Output</summary>
```
	$ cat ph_summary.txt
	PH_num_hap 3
	PH_hap_diversity 0.611668153059
	PH_seq_len 1208


	$ cat ph_haplotypes.fna
	>reconstructed_0 Freq=0.190638
	ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAATATCTGGCCTTCCCGCAAGGGAAG
	GCCAGGGAACTTTCCTCAGAGCAGGCTAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGCTTGGGGAA
	GCAACAACAAAGCCCTCTCAGAAACAGGAGCCAATAGACAAGGAAATGTATCCTTTAACTTCCCTCAGAT
	CACTCTTTGGCAACGACCCCTTGTCACAATAAAGATAGGGGAGCAACTAAAGGAAGCCCTGTTAGATACA
	GGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAGAATGATAGGGGGAA
	TCGGGGGTTTTATCAAAGTAAGACAGTATGATCAGATAGCCATAGACATCTGTGGCCATAAAGCTATAGG
	TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGCTTGGTTGCACT
	TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTGAAGCCAGGAATGGATGGCCCAAAGG
	TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAAGA
	AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGGAAAAAAGAC
	AGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAGGACTTCTGGGAAGTCC
	AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGCGGGTGATGC
	ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTATACTGCATTTACCATACCTAGTGTAAACAAT
	GAGACACTAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
	AAGCTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACAT
	GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
	CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
	GTTATGAACTCCATCCTT
	>reconstructed_1 Freq=0.294104
	ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCCGCGGCGGAAG
	GCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAG
	GAGACAACAACTCCCTCTCAGAAGCAGGAGCCGAAGGACAAGGAAACATATCCTTTAGCTTCCCTCAAAT
	CACTCTTTGGCAACGACCCCTTGTTACAATAAAGATAGGGGGGCAGCTAAAGGAAGCTCTATTAGATACA
	GGAGCAGATGACACAGTATTAGAAGAAATGGATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAA
	TTGGAGGTTTTATCAAAGTAAAACAGTATGATCAGATACCCATAGAAATTTGTGGACATAAAGTGATAGG
	TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGCTTGGTTGCACT
	TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTGAAGCCAGGAATGGATGGCCCAAAGG
	TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAATAGAAATCTGTGCAGAAATGGAAAAGGA
	AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGAAAAAAAGAC
	AGTACTAAATGGAGAAAATTAGTAGATTTCAAAGAACTTAATAAGAGAACTCAGGACTTCTGGGAAGTCC
	AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAGTCAGTAACAGTACTGGATGTGGGTGATGC
	ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTACACTGCATTTACCATACCTAGTGTAAACAAT
	GAGACACCAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
	AATGTAGCATGACAAAAATCTTAGAACCTTTTAGAAAACAAAATCCAGATATAGTTATCTATCAATACAT
	GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
	CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
	GTTATGAACTCCATCCTT
	>reconstructed_2 Freq=0.515259
	ACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAAATCTGGCCTTCCCACAAGGGAAG
	GCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAA
	GAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGAAATGTATCCTTTAGCTTCCCTCAAAT
	CACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACA
	GGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAA
	TTGGAGGTTTTATCAAAGTAAAACAGTATGATCAGATACCCATAGAAATCTGTGGACATAAAGCTATAGG
	TACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACT
	TTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAG
	TTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGA
	AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAGGAAAAAAGAC
	AGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCCGGGAAGTCC
	AATTAGGAATACCACATCCCTCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGC
	ATATTTTTCAGTTCCCTTAGATGAAGACTTCAGAAAGTATACTGCATTTACCATACCTAGTGTAAACAAT
	GAGACACCAGGGATTAGGTATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCAGCAATATTCC
	AAGCTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACAT
	GGATGATCTGTATGTAGGATCTGACCTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGAGAA
	CATCTGTTGAGGTGGGGGTTTTGCACACCAGACAAGAAACATCAGAAGGAACCTCCATTCCTTTGGATGG
	GTTATGAACTCCATCCCC
```
</details>

			
