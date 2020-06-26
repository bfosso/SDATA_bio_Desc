# Passaro_et_al 2020

* [Rationale](#rationale)
* [Requirements](#requirements)  
    - [Bioinformatic tools and packages](#bioinformatic-tools-and-packages)
    - [Required data files](#required-data-files)  
* [Bioinformatic Workflow](#bioinformatic-workflow)  
    1. [Taxonomic assignment of Illumina PE reads](#1-taxonomic-assignment-of-illumina-pe-reads)  
        * [Raw data retrieval](#raw-data-retrieval)
        * [Metashot application](#metashot-application)
        * [Unassigned PE reads extraction](#unassigned-pe-reads-extraction)
    2. [Meta-assembly of unassigned reads](#2-meta-assembly-of-unassigned-reads)
        * [Unassigned data retrival](#unassigned-pe-reads-extraction)
        * [Human reads removal](#human-reads-removal)  
        * [Metagenome Assembly](#metagenome-assembly)  
    3. [Taxonomic assignments of the obtained contigs/scaffolds](#3-taxonomic-assignments-of-the-obtained-contigsscaffolds)
        * [Scaffold data retrieval](#scaffolds-data-retrieval)
-------------
##Rationale
In this repository are described all the bioinformatic steps discussed in [Passaro et al 2019](https://www.nature.com/articles/s41598-019-56240-1). 
The study was focused on the metagenomic investigation of tumors samples in order to identify putative oncoviruses in immunosuppressed patients. Consistently with the major findings of several recent papers no new human tumorigenic viruses were identified. 
The raw data mentioned in the paper are available in the *SRA* repository under the [**Bioproject PRJNA544407**](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA544407).
In particular, the 13 biological samples used in this study were tumors ablated for therapeutic purposes from 12 patients (Table 1).  
DNA or RNA were extracted and, according to their quality, sequenced by using the **Illumina NestSeq 500** platform. Only for patients T7 we were able to obtain both high quality DNA and RNA.  

*Table1: Sample metadata and SRA data references.*  

|Code|Tumor type|Nucleic acid sequenced|Immunosuppressive condition (IC)|Years from onset of IC|BioSample ID|SRA ID|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|T1|Skin squamous cell carcinoma|RNA|Renal transplantation, immunosuppressive therapy|20|SAMN11835442|SRR10202451|
|T5|Native kidney (oncocytoma)|RNA|Renal transplantation, immunosuppressive therapy|19|SAMN11835462|SRR10202450|
|T7|Transplanted kidney (clear cell carcinoma)|DNA|Renal transplantation, immunosuppressive therapy|3|SAMN11835463|SRS5453517|
|T7_RNA|Transplanted kidney (clear cell carcinoma)|RNA|Renal transplantation, immunosuppressive therapy|3|SAMN12838684|SRR10202444|
|T8|Native kidney (oncocytoma)|DNA|Renal transplantation, immunosuppressive therapy|20|SAMN11835464|SRR10202443|
|T9|Non-Hodgkin Lymphoma|DNA|Renal transplantation, immunosuppressive therapy|12|SAMN11835496|SRR10202442|
|T10|Colon adenocarcinoma|DNA|Renal transplantation, immunosuppressive therapy|5|SAMN11835497|SRR10202441|
|T11|Native kidney (clear cell carcinoma)|RNA|Renal transplantation, immunosuppressive therapy|7|SAMN11835498|SRR10202440|
|T12|Skin carcinomas|RNA|Renal transplantation, immunosuppressive therapy|12|SAMN11835499|SRR10202439|
|T13|Skin carcinomas|RNA|Renal transplantation, immunosuppressive therapy|12|SAMN11835501|SRS5453523|
|T14|Skin squamous cell carcinoma|RNA|Renal transplantation, immunosuppressive therapy|8|SAMN11835502|SRR10202438|
|N4|Carcinoma of the tongue and oropharynx|RNA|Non-Hodgkin lymphoma|15|SAMN11835504|SRR10202448|
|N6|Lip squamous cell carcinoma (HPV-neg.)|RNA|Acute lymphocytic leukemia|11|SAMN11835505|SRR10202447|

The obtained sequencing data were analysed by applying a bioinformatic pipeline relying on 3 main steps:  
    1. Taxonomic assignment of Illumina PE reads by exploiting **MetaShot**;  
    2. Meta-assembly of unassigned reads;  
    3. Taxonomic assignments of the obtained contigs/scaffolds.  
The described procedure allows users to replicate the whole procedure or just reproduce on of the steps and the intermediate data are available as a [**Zenodo**]() repository.  

## Requirements
### Bioinformatic tools and packages
All the steps described below rely on several tools and packages whose installation and configuration is required to properly reproduce all the listed steps.  
Following the list of required tools:  
  * [**MetaShot (Metagenomics Shotgun)**](https://github.com/bfosso/MetaShot) \[[PMID: 28130230](https://pubmed.ncbi.nlm.nih.gov/28130230/)\] is a pipeline designed for the complete taxonomic assessment of the human microbiota. 
  In MetaShot, third party tools and *ad hoc* developed `Python` and `Bash` scripts are integrated to analyse *paired-end (PE) Illumina reads*, offering an automated procedure covering all the steps from raw data management to taxonomic profiling. 
  It is designed to analyse both DNA-Seq and RNA-Seq data.   
  * [**metaSPAdes**](http://cab.spbu.ru/software/meta-spades/) \[[PMID: 28298430](https://pubmed.ncbi.nlm.nih.gov/28298430/)\] is an assembler designed to obtain high quality metegenomes assemblies.  
  * [**WindowMasker**](https://www.ncbi.nlm.nih.gov/toolkit) \[[PMID: 16287941](https://pubmed.ncbi.nlm.nih.gov/16287941/)\]  identifies and masks highly repetitive DNA sequences in a contigs/scaffolds. It is included in the NCBI C++ toolkit.  
  * [**RepeatMasker**](http://www.repeatmasker.org/) \[[PMID: 19274634](https://pubmed.ncbi.nlm.nih.gov/19274634/)\] allows to identify, classify, and mask repetitive elements, including low-complexity sequences and interspersed repeats.  
  * [**BLAST+**](https://www.ncbi.nlm.nih.gov/books/NBK279690/) \[[PMID: 20003500](https://pubmed.ncbi.nlm.nih.gov/20003500/)\] inds regions of local similarity between sequences.  
  * [**Custom Perl scripts**](https://github.com/matteo14c/Passaro_et_al) developed by [Dr. Matteo Chiara](mailto:matteo.chiara@unimi.it).  
  * [**bowtie2**](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) \[[PMID: 22388286](https://pubmed.ncbi.nlm.nih.gov/22388286/)\] aligns Illumina PE reads to references genomes.  
  * [**SRA toolkit**](https://github.com/ncbi/sra-tools) is a suite of tools allowing to access the *INSDC* content.  
  * [**Zenodo_get**](https://gitlab.com/dvolgyes/zenodo_get) is a dowloader for Zenodo records.  

### Required data files
1. The link to the MetaShot reference data is available in its *github* repository, where it is also available a description of the required configurations.  
    ***Please note that the download and configuration of MetaShot reference data is required only in case you want to reproduce the whole pipeline.***   
2. The Meta-assembly steps begins with the mapping on the human genome, in order to remove human sequences MetaShot was unable to identify:  
        * If you have downloaded the MetaShot reference data it is also contains the **human genome assembly hg19** in `MetaShot_reference_data/Homo_sapiens/`;  
        * Otherwise you may download it from [**UCSC**](http://hgdownload.soe.ucsc.edu/downloads.html#human) or [**ENSEMBLE**](http://grch37.ensembl.org/Homo_sapiens/Info/Index).      
    In order to build the *hg19* **bowtie2** indexes, type the following line:
    ```
    bowtie2-build -f hg19.fa hg19
    ```
    `bowtie2-build` is part of the **bowtie2** package and builds a Bowtie index from a set of DNA sequences. The `-f` option indicates the input sequences are in *FASTA* format.    
5. The *FASTA* file containing all the human **RefSeq** transcript can be downloaded by using the [**UCSC TABLE browser**](https://genome.ucsc.edu/cgi-bin/hgTables). The relative bowtie2 indexes were obtained by typing:
    ```
    bowtie2-build -f refseq.fa refseq
    ```
4. The human micro-satellites sequences for the hg19 assembly were obtained from the [**UCSC Genome Browser**](http://genome.ucsc.edu/) by using the table browser tool.   
5. The list of human repeats was retrieved from the [**GIRI (Genetic Information Research Institute) Repbase**](https://www.girinst.org).    

## Bioinformatic Workflow
For reproducibility purposes, sequencing data were deposited as raw reads. 
Nonetheless, considering that the most intense and computational expensive steps were the taxonomic classification and the metagenomes assemblies performed by MetaShot and metaSPAdes, respectively, and in order to facilitate the analysis reproducibility by avoiding to repeat one or both these steps, the **unassigned PE reads** and the **assembled scaffolds** are available in a [*Zenodo*]() repository:  
* [MetaShot unassigned read]());  
* [Scaffolds]()).  


### 1. Taxonomic assignment of Illumina PE reads
#### Raw data retrieval
To begin the analysis from this step, raw data download is required. To retrieve the *FASTQ* files from SRA, just use the `fastq-dump` tool from the **SRA toolkit** suite.  
To download the sample *N6* data type the following line:
```
fastq-dump --split-files -O SRR10202447 `SRR10202447` 
```
This will generate a folder called `SRR10202447` containing 2 fastq files: `SRR10202447_1.fastq` and `SRR10202447_2.fastq`.  

#### Metashot application
The whole MetaShot workflow was performed by typing the following command:  
```
MetaShot_Master_script.py -m read_list.tsv -p parameters_file
```
In particular:
- `-m` refers to a *tsv* file containing a list of PE reads files. Have a look at the [guide](https://github.com/bfosso/MetaShot#Usage);   
- `-p` refers to a structured file containing all the info MetaShot needs to perform the analysis. For more info, have a look at the *MetaShot setting up* section in the [README](https://github.com/bfosso/MetaShot#metashot-setting-up) file.   

It produces several files and folders but the most important are:  
* *ambiguos_pe_read.lst*: textual fine containing all the ambiguous PE reads (reads mapping on more than one reference division);  
* *bacteria_CSV_result.csv*: tabular file summarising the taxonomic assignments for prokaryotes;  
* *fungi_CSV_result.csv*: tabular file summarising the taxonomic assignments for fungi;  
* *protist_CSV_result.csv*: tabular file summarising the taxonomic assignments for protist;  
* *virus_CSV_result.csv*: tabular file summarising the taxonomic assignments for viruses.  

A more extensive description about MetaShot results is available [here](https://github.com/bfosso/MetaShot#result-files-interpretation).  
#### Unassigned PE reads extraction
Following, the unassigned reads were extracted by using the `PE_extraction.py` script:  
```
PE_extraction.py -u
```
This script generates a folder containing all the unassigned reads.  
The folder name is automatically generated by using the following format `taxid_PE_file_DAY_MONTH_YEAET_HOUR_MINUTES_SECONDS` in order to avoid data overwriting.


### 2. Meta-assembly of unassigned reads
#### Unassigned data retrieval
If you want to skip the MetaShot step, just download the unassigned reads available in the **Zenodo** repository by using *zenodo_get*.  
There are two ways to download the files:
1. retrieve the whole repository:  
```
    mkdir zenodo_<RECORD_ID> && cd zenodo_<RECORD_ID>
    zenodo_get 3893846
```  
   This will automaticaly download the whole repository content and create a `md5sums.txt`file to crosscheck if the files were properly retrieved.
2. retrieve selected files. Suppose we are interested in retrieve just the *N6* unassigned reads:
```
    mkdir N6 && cd N6
    zenodo_get -w `files_link_list` zenodo_<RECORD_ID>
```
   This will automatically create a file named `files_link_list` containing the **https** link to each file in the repository.  
   Following by using `wget` you can retrieved the desired data.  
```
    wget <LINK_https>
```   
#### Human reads removal
Before to perform the reads assembly, we have remove *human* reads MetaShot was unable to identify.  
```
    bowtie2 -1 <R1files> -2 <R2files> -x /path/to/hg19_bowtie_index/hg19 --very-sensitive-local -p 12 -S <name.sam> --un-conc <name\_nonhumanPE>
```
In particular:
    * **-1**: file(s) containing the R1 reads;  
    * **-2**: file(s) containing the R2 reads;  
    * **-x**: bowtie2 indexes;  
    * **--very-sensitive-local**: mapping preset giving prority to sensitivity;  
    * **-p**: number of available processors;  
    * **-S**: SAM file output name;   
    * **--un-conc**: by using this option 2 files containing unmapped R1 and R2 reads were generated and named *name\_nonhumanPE*.  

In case of RNA-Seq data a mapping against *RefSeq* human transcripts was also performed:  
```
    bowtie2 -1 <name\_nonhumanPE.1> -2 <name\_nonhumanPE.2> -x /home/mchiara/refseq --very-sensitive-local -p 12 -S <name.sam> --un-conc <name\_RNA_nonhumanPE>
```
#### Metagenome Assembly
The metagenome assembly was perfomed by using **metaSPAdes**.  
Please note that for DNA-Seq and RNA-SEQ data were assembled the obtained *name\_nonhumanPE* and *name\_RNA_nonhumanPE* files, respectively.  
```
metaspades.py -1 <name\_nonhumanPE.1>  -2 <name\_nonhumanPE.2> -t 12  -k 21,33,55,77,99 -o <name_meta>
```
In particular:
    * **-1**: file(s) containing the R1 reads;  
    * **-2**: file(s) containing the R2 reads;  
    * **-t**: number of available threads;  
    * **-k**: k-mer dimensions;  
    * **-o**: output folder name.
    
Following in order to mask human repeat (we used the `species human` ) we applied **RepeatMasker** on the obtained contigs: 
```
RepeatMasker -species human <name_meta>
```
At the end of the analysis a file with the `.masked` suffix was generated.  
**WindowMasker** allowed to indentify and mask low complexity and highly repetitive sequences.  
Two steps were required:
1.  Words occurrence inference: it produced the `name_meta.MK` file containing the word occurrence.
```
    windowmasker -in <name_meta.masked> -mk_counts > <`name_meta.MK`>
```
2. Masking: It masked the repetitive words and low-complexity regions by using the word count file generated in the previous step and the dust algorithm, respectively.  
```
    windowmasker -in <name_meta.masked> -ustat <name_meta.MK> -dust T -outfmt fasta > <name_meta.double-masked.fasta>
```
Finally, we removed contings containing more than 15% of *N*: 
```
perl filter.pl <name_meta.double-masked.fasta> > <name_meta.BLAST.fasta>
```

## 3. Taxonomic assignments of the obtained contigs/scaffolds  
#### Scaffolds data retrieval
If you want to skip the MetaShot step, just download the unassigned reads available in the **Zenodo** repository by using *zenodo_get*.  
There are two ways to download the files:
1. retrieve the whole repository:  
```
    mkdir zenodo_<RECORD_ID> && cd zenodo_<RECORD_ID>
    zenodo_get 3893846
```  
   This will automaticaly download the whole repository content and create a `md5sums.txt`file to crosscheck if the files were properly retrieved.
2. retrieve selected files. Suppose we are interested in retrieve just the *N6* unassigned reads:
```
    mkdir N6 && cd N6
    zenodo_get -w `files_link_list` zenodo_<RECORD_ID>
```
   This will automatically create a file named `files_link_list` containing the **https** link to each file in the repository.  
   Following by using `wget` you can retrieved the desired data.  
```
    wget <LINK_https>
```   
The retained contigs were taxonomically classified by using the **blastn**. The ``-remote`` options allows to query remote blast db available on the NCBI servers.  
```
blastn -remote -query  <name_meta.BLAST.fasta> -db nr > name_meta_BLAST.res
```

Scaffolds were assigned to the blastn best match if it covered at least the 30% of the query sequence with a similarity percentage equal or higher than 70%, by using the following command:
```
perl simple.parse.blast.pl G <name_meta_BLAST.res>
```
