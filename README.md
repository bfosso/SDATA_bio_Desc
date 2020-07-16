# Passaro_et_al 2020
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3948399.svg)](https://doi.org/10.5281/zenodo.3948399)  

* [Rationale](#rationale)
* [Requirements](#requirements)  
    - [Tools and packages](#tools-and-packages)
    - [Required data files](#required-data-files)  
* [Bioinformatic Workflow](#bioinformatics-workflow)  
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
        * [Masking of low complexity and repeated regions in the scaffolds](#masking-of-low-complexity-and-repeated-regions-in-the-scaffolds)
        * [Scaffolds taxonomic classification](#scaffolds-taxonomic-classification) 
-------------
## Rationale
This repository provides all the bioinformatics tools and workflows used to analyse the data in [Passaro et al 2019](https://www.nature.com/articles/s41598-019-56240-1).  
The aim of the study was to perform investigation of tumors samples by metagenomics analyses, in order to identify putative oncoviruses in immunosuppressed patients. Consistently with the major findings of several recent papers no novel human tumorigenic viruses could be identified. 
The 13 biological samples used in this study were tumors ablated for therapeutic purposes from 12 patients (*Table 1*).
DNA or RNA were extracted and, according to their quality, sequenced by using the **Illumina NestSeq 500** platform and a **Paired-Ends (PE)** layout. Only for patients **T7** we were able to obtain both high quality DNA and RNA.
The raw data are available in the *SRA* repository under the [**Bioproject PRJNA544407**](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA544407).

*Table 1: Sample metadata and SRA data references.*  

|Code|Tumor type|Nucleic acid sequenced|Immunosuppressive condition (IC)|Years from onset of IC|BioSample ID|SRA ID|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|T1|Skin squamous cell carcinoma|RNA|Renal transplantation, immunosuppressive therapy|20|SAMN11835442|SRR10202451|
|T2|Skin squamous cell carcinoma|RNA|Renal transplantation, immunosuppressive therapy|9|SAMN11835461|SRR10202450|
|T5|Native kidney (oncocytoma)|RNA|Renal transplantation, immunosuppressive therapy|19|SAMN11835462|SRR10202445|
|T7|Transplanted kidney (clear cell carcinoma)|DNA|Renal transplantation, immunosuppressive therapy|3|SAMN11835463|SRR10202444|
|T7_RNA|Transplanted kidney (clear cell carcinoma)|RNA|Renal transplantation, immunosuppressive therapy|3|SAMN12838684|SRR10202446|
|T8|Native kidney (oncocytoma)|DNA|Renal transplantation, immunosuppressive therapy|20|SAMN11835464|SRR10202443|
|T9|Non-Hodgkin Lymphoma|DNA|Renal transplantation, immunosuppressive therapy|12|SAMN11835496|SRR10202442|
|T10|Colon adenocarcinoma|DNA|Renal transplantation, immunosuppressive therapy|5|SAMN11835497|SRR10202441|
|T11|Native kidney (clear cell carcinoma)|RNA|Renal transplantation, immunosuppressive therapy|7|SAMN11835498|SRR10202440|
|T12|Skin carcinomas|RNA|Renal transplantation, immunosuppressive therapy|12|SAMN11835499|SRR10202439|
|T13|Skin carcinomas|RNA|Renal transplantation, immunosuppressive therapy|12|SAMN11835501|SRR10202438|
|T14|Skin squamous cell carcinoma|RNA|Renal transplantation, immunosuppressive therapy|8|SAMN11835502|SRR10202449|
|N4|Carcinoma of the tongue and oropharynx|RNA|Non-Hodgkin lymphoma|15|SAMN11835504|SRR10202448|
|N6|Lip squamous cell carcinoma (HPV-neg.)|RNA|Acute lymphocytic leukemia|11|SAMN11835505|SRR10202447|

The obtained sequencing data were analysed by applying a bioinformatics pipeline relying on 3 main steps:  
    1. Taxonomic assignment of Illumina PE reads by exploiting **MetaShot**;  
    2. Meta-assembly of unassigned reads;  
    3. Taxonomic assignments of the obtained scaffolds.  
The described method allows users to replicate the whole procedure or just reproduce one specific step. The intermediate data are available as a [**Zenodo**]() repository.  

## Requirements
### Tools and packages
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
  * [**Zenodo_get**](https://gitlab.com/dvolgyes/zenodo_get) is a downloader for Zenodo records.  
  * [**samtools**](https://github.com/samtools/samtools) \[[PMID: 19505943](https://pubmed.ncbi.nlm.nih.gov//19505943/)\] is a suite of tools for the manipulation of files in Sequence Alignment/Map format (sam). 

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
    `bowtie2-build` is part of the **bowtie2** package and builds a bowtie2 index from a set of DNA sequences. The `-f` option indicates the input sequences are in *FASTA* format.    
5. The *FASTA* file containing all the human **RefSeq** transcript can be downloaded by using the [**UCSC TABLE browser**](https://genome.ucsc.edu/cgi-bin/hgTables). To create the relative bowtie2 indexes type as following:
    ```
    bowtie2-build -f refseq.fa refseq
    ```
4. The human micro-satellites sequences for the hg19 assembly were obtained from the [**UCSC Genome Browser**](http://genome.ucsc.edu/) by using the table browser tool.   
5. The list of human repeats was retrieved from the [**GIRI (Genetic Information Research Institute) Repbase**](https://www.girinst.org). **Due to the Repbase data sharing policy you need to request the data on your-own*.  

## Bioinformatics Workflow
For reproducibility purposes, sequencing data were deposited as raw reads. 
Nonetheless, considering that the most intense and computational expensive steps were the taxonomic classification and the metagenomes assemblies performed by MetaShot and metaSPAdes, respectively, and in order to facilitate the analysis reproducibility by avoiding to repeat one or both these steps, the **unassigned PE reads** and the **assembled scaffolds** are available in a [*Zenodo*]() repository.  


### 1. Taxonomic assignment of Illumina PE reads
*Following it is described how to use MetaShot. If you want to skip this step, just jump to the [Meta-assembly of unassigned reads](#2-meta-assembly-of-unassigned-reads).*

#### Raw data retrieval
To begin the analysis from this step, raw data download is required. To retrieve the *FASTQ* files from SRA, just use the *fastq-dump* tool from the **SRA toolkit** suite.  
For instance, to download the sample *N6* data type the following line:
```
fastq-dump --split-files -O SRR10202447 `SRR10202447` 
```
This will generate a folder called `SRR10202447` containing 2 fastq files: `SRR10202447_1.fastq` and `SRR10202447_2.fastq`.  

#### Metashot application
The whole MetaShot workflow can be performed by typing the following command:  
```
MetaShot_Master_script.py -m read_list.tsv -p parameters_file
```
In particular:
- `-m` refers to a *tsv* file containing a list of PE reads files. Have a look at the [guide](https://github.com/bfosso/MetaShot#Usage);   
- `-p` refers to a structured file containing all the info MetaShot needs to perform the analysis. For more info, have a look to [MetaShot setting up](https://github.com/bfosso/MetaShot#metashot-setting-up) guide.   
Metashot will performs several steps:
1. removal of Phix reads;  
2. trimming of low quality and low-complexity reads;  
3. Mapping on the host-genome;  
4. Comparison with prokaryotic, fungal, viral and protist reference collections;
5. removal of ambiguous reads (i.e. reads mapping on more than one reference collection);  
6. Taxonomic classification of unambiguous reads;  
7. Report preparation.

It produces several files and folders but the most important are:  
* *ambiguos_pe_read.lst*: textual fine containing all the ambiguous PE reads (reads mapping on more than one reference division);  
* *bacteria_CSV_result.csv*: tabular file summarising the taxonomic assignments for prokaryotes;  
* *fungi_CSV_result.csv*: tabular file summarising the taxonomic assignments for fungi;  
* *protist_CSV_result.csv*: tabular file summarising the taxonomic assignments for protist;  
* *virus_CSV_result.csv*: tabular file summarising the taxonomic assignments for viruses.  

A more extensive description about MetaShot results is available [here](https://github.com/bfosso/MetaShot#result-files-interpretation).  
#### Unassigned PE reads extraction
Following, the unassigned reads are extracted by using the `PE_extraction.py` script:  
```
PE_extraction.py -u
```
This script generates a folder containing all the unassigned reads.  
The folder name is automatically generated by using the following format `taxid_PE_file_DAY_MONTH_YEAET_HOUR_MINUTES_SECONDS` in order to avoid data overwriting.


### 2. Meta-assembly of unassigned reads
#### Unassigned data retrieval
If you have skipped the MetaShot step, just download the unassigned reads available in the **Zenodo** repository by using *zenodo_get*.  
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
Before to perform the reads assembly, we need to remove *human* reads MetaShot was unable to identify.  
```
    bowtie2 -1 <R1files> -2 <R2files> -x /path/to/hg19_bowtie_index/hg19 --very-sensitive-local -p 12 -S <name>.sam --un-conc <name>_nonhumanPE
```
In particular:
    * **-1**: file(s) containing the forward reads;  
    * **-2**: file(s) containing the reverse reads;  
    * **-x**: bowtie2 indexes;  
    * **--very-sensitive-local**: mapping preset giving priority to sensitivity and performing the alignment in local mode;  
    * **-p**: number of available processors;  
    * **-S**: SAM file output name;   
    * **--un-conc**: by using this option 2 files containing unmapped R1 and R2 reads are generated and named *name_nonhumanPE.1* and *name_nonhumanPE.2*. Just replace `<name>` with the sample name.   

In case of RNA-Seq data a mapping against *RefSeq* human transcripts needs also to be performed:  
```
    bowtie2 -1 <name>_nonhumanPE.1 -2 <name>_nonhumanPE.2 -x /home/mchiara/refseq --very-sensitive-local -p 12 -S <name>.sam --un-conc <name>_RNA_nonhumanPE
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
    
## 3. Taxonomic assignments of the obtained contigs/scaffolds  
#### Scaffolds data retrieval
If you want to skip the metaSPAdes assembly step, just download the scaffolds available in the **Zenodo** repository by using *zenodo_get*.  
As mentioned [before](#unassigned-pe-reads-extraction) you can retrieve the whole repository or just some selected files.  
    
Following in order to mask human repeat (we used the `species human` ) we applied **RepeatMasker** on the obtained contigs: 

#### Masking of low complexity and repeated regions in the scaffolds
In order to save time we mask repeated and low complexity regions in the scaffolds by using *RepeatMasker* and *WindowMasker*, respectivaly.  
By using *RepeatMasker* we can mask known human repetitive elements:  
```
RepeatMasker -species human <name_meta>
```
At the end of the analysis a file with the `.masked` suffix will be produced.  
*WindowMasker* allows to identify and mask low complexity and highly repetitive sequences.  
Two steps are required:
1.  *Words occurrence inference*: it stores the the words occurrence in the `name_meta.MK` file.
```
    windowmasker -in <name_meta.masked> -mk_counts > <`name_meta.MK`>
```
2. *Masking*: It masks the repetitive words and low-complexity regions by using the words count file generated in the previous step and the dust algorithm, respectively.  
```
    windowmasker -in <name_meta.masked> -ustat <name_meta.MK> -dust T -outfmt fasta > <name_meta.double-masked.fasta>
```
Finally, we removed contings containing more than 15% of *N* (corresponding to masked nucleotides): 
```
perl filter.pl <name_meta.double-masked.fasta> > <name_meta.BLAST.fasta>
```

#### Scaffolds taxonomic classification
The retained contigs were taxonomically classified by using the *blastn*. The -remote options allows to query remote blast db available on the NCBI servers.
```
blastn -remote -query  <name_meta.BLAST.fasta> -db nr > name_meta_BLAST.res
```
Scaffolds were assigned to the blastn best match if it covered at least the 30% of the query sequence with a similarity percentage equal or higher than 70%, by using the following command:  
```
perl simple.parse.blast.pl G <name_meta_BLAST.res>
```
The final output consists in a simple table, where for every species to which one or more contigs were assigned, the total number of contigs assigned to that species, and their total size is reported.  
An example is enclosed in *Table 2*:  

*Table 2: result example.*  

|Species Name|# of assigned contigs|total size|
|:---:|:---:|:---:|
|human|513|301243|
|chlorocebus|1|371|
|Nohit|16|5125|
|monkey|90|41487|
|onchocerca_flexuosa|1|127|
|rhesus_macaque|4|3613|

Nohit is used to indicate sequence that show no significant similarity/were not assigned to any species.

Only for scaffolds derived from the assembly of metatranscriptomic data (RNA samples, namely T1,T5,T7_RNA,T11,T12,T13,T14,N4,N6), an additional round of taxonomic assignment was performed by a sequence similarity search against the Gencode V31 annotation of the human genome. A copy of the fasta file `Gencode_V31_transcripts.fa` is available in this repository.
The following commands are required to execute this analysis.
First a blast nucleotide database needs to be created from the fasta file. This can be done as follows
```
gzip -d Gencode_V31_transcripts.fa.gz 

makeblastdb -in Gencode_V31_transcripts.fa -dbtype nucl -out Gencode_V31
```
The makeblastdb command should be included in any standard blast+ package installation.   
Sequence similarity searches, were performed as outlined previously by means of the blastn command
```
blastn -query <name_meta.BLAST.fasta> -db Gencode_V31 > name_meta_BLAST_transcriptome.res
```
Please notice that the name of the "-db" argument of the **blastn** command needs to match exactly the name of the blast database that you should have created with **makeblastdb** (that is the name of the -out argument).  

Finally assignment of metatranscriptomic contigs, is performed by parsing the output file of blastn by means of *simple.parse.blast.pl*
```
perl simple.parse.blast.pl G <name_meta_BLAST_transcriptome.res>
```
#### Identification of possible integration sites of viral genomes

For the 3 samples (namely T1, T8 and N6) in which viral specimens were detected by our analyses, we applied an *ad-hoc* method to identify possible sites of integration of such viruses in the genome of the host.
For this analysis, unassigned metagenomic reads were mapped by using bowtie2 on a custom sequence database (called index) containing the metagenomic assemblies of the viral isolates identified in our the paper. The corresponding file, in fasta format, can be obtained github repository with our custom Perl scripts: [this repository ](https://github.com/matteo14c/Passaro_et_al) .  The file is called *viral_seq_Passaro_et_al.fa*
To obtain a bowtie2 index file from the fasta, the following command is required
```
bowtie2-build <viral_seq_Passaro_et_al.fa> viraldb
```
To map unassigned metagenomic reads to these viral scaffolds, the following command was used with bowtie2:
```
bowtie2 -1 <R1files> -2 <R2files> -x viraldb --very-sensitive-local -p 12 -S <name>.sam 
```
Please notice that *viraldb* here denotes the database of viral metagenomic scaffold as obtained by the *bowtie2-build* command.

After mapping the reads, a custom Perl script called *parse_Vir_map.pl* ,which is available from [this repository ](https://github.com/matteo14c/Passaro_et_al), was applied to the output file *\<name>.sam* to retrieve reads with partial similarity to a viral genome assembly (i.e incomplete mapping) or pairs of reads for which only one mate of the pair could be confidently mapped to to a viral scaffolds. By further mapping these reads to the reference hg19 human genome assembly, we searched for possible sites of integration of the viruses in the genome of the host.  The presence of single reads or pairs of reads with partial similarity to both genomes (the virus and the host genome) indeed is considered to be indicative of integration of the virus in the host genome. 
The command for *parse_Vir_map.pl* is as follows:
```
perl parse_Vir_map.pl name.sam
```
The program produces 2 output files:
    * **-1**: one, with the suffix "singleton.fq" contains "orphan" reads, that is the reads in a pair where only the mate of the read, but not the read itself could be assigned confidently to a viral scaffold   
    * **-2**: a second file,with the suffix: "partial.fq" contains reads that have only a partial similarity to a viral genome scaffold (greater than 25% but smaller than 75% of the read size)
   
   At this point Sites/events of possible integration of viral genomes in the human host genome can be inferred by mapping back *\_singleton.fq* and *\_partial.fq* to the hg19 human genome reference assembly using bowtie2. 
   For example with these commands
```
bowtie2 --U name.sam_singleton.fq -x path/to/hg19_bowtie_index/hg19  --very-sensitive-local -S singleton.sam 

bowtie2 --U name.sam_partial.fq -x path/to/hg19_bowtie_index/hg19  --very-sensitive-local -S partial.sam 
```      
In the example, the results of the mapping, will be stored in 2 files: *singleton.sam* and *partial.sam*, both in sam format.
In Passaro et al, no evidence of integration was observed, and no read or pairs or reads showed hints of possible cross mapping on the reference and host genome.
To check the number of reads mapped to the reference hg19 assembly, you can simply use *samtools flagstat*:
```
    samtools flagstat singleton.sam 
```
```
    samtools flagstat partial.sam 
```
Should you find any reads mapped to hg19, if you want to extract them from the corresponding sam file, to have an indication of the possible loci of integration,  you can again use the samtools to extract these reads. The command would be something like:
```
    samtools view -F 4 partial.sam 
```
```
    samtools view -F 4 singleton.sam 
```

If you read up to this point, this means that either you have successfully completed the workflow, or that you just skipped to the last line. In the first case **congrats!** in the latter, don't worry I usually do that to!
