# Passaro_et_al 2020
This repository describes all the bioinformatic steps performed in [Passaro et al 2019](https://www.nature.com/articles/s41598-019-56240-1) and **NNNNNNNN**.  
The raw data mentioned in the paper are available in the *SRA* repository under the Bioproject [PRJNA544407](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA544407).  
The bioinformatic pipeline relies on 3 main steps:  
1. Taxonomic assignment of Illumina PE reads by exploiting **MetaShot**;  
2. Meta-assembly of unassigned reads;  
3. Taxonomic assignments of the obtained contigs/scaffolds.  
The described procedure allows users to replicate the whole procedure or just reproduce on of the steps and the intermediate data are available as a **Zenodo** repository.  

## Requirements
### Bioinformatic tools and packages
All the steps described below rely on several bioinformatic tools and packages whose installation and configuration is required to properly reproduce all the listed steps.  
Following the list of required tools:  
  * [**MetaShot (Metagenomics Shotgun)**](https://github.com/bfosso/MetaShot) \[[PMID: 28130230](https://pubmed.ncbi.nlm.nih.gov/28130230/)\] is a pipeline designed for the complete taxonomic assessment of the human microbiota. In MetaShot, third party tools and new developed Python and Bash scripts are integrated to analyse *paired-end (PE) Illumina reads*, offering an automated procedure covering all the analysis steps from raw data management to taxonomic profiling. It is designed to analyse both DNA-Seq and RNA-Seq data. Among the tools required in MetaShot, **bowtie2** is also listed, so it can be installed only once.   
  * [**metaSPAdes**](http://cab.spbu.ru/software/meta-spades/) \[[PMID: 28298430](https://pubmed.ncbi.nlm.nih.gov/28298430/)\] is an assembler designed to obtain high quality assembly of metegenomes.  
  * [**WindowMasker**](https://www.ncbi.nlm.nih.gov/toolkit) \[[PMID: 16287941](https://pubmed.ncbi.nlm.nih.gov/16287941/)\]  identifies and masks highly repetitive DNA sequences in a contigs/scaffolds. It is included in the NCBI C++ toolkit.  
  * [**RepeatMasker**](http://www.repeatmasker.org/) \[[PMID: 19274634](https://pubmed.ncbi.nlm.nih.gov/19274634/)\] allows to identify, classify, and mask repetitive elements, including low-complexity sequences and interspersed repeats.  
  * [**Custom Perl scripts**](https://github.com/matteo14c/Passaro_et_al) developed by [Dr. Matteo Chiara](mailto:matteo.chiara@unimi.it)  
  * [**bowtie2**](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) \[[PMID: 22388286](https://pubmed.ncbi.nlm.nih.gov/22388286/)\] aligns Illumina PE reads to references genomes. This in the only tool requiring installation if the MetaShot execution is skipped.
  * [**SRA toolkit**](https://github.com/ncbi/sra-tools) is a suite of tools allowing to access the INSDC content.  

### Required data files
The link to the MetaShot reference is available its *github* repository, where it is also available a description of the required configurations.  
*Please note that the download and configuration of MetaShot reference data is required only in case you want to reproduce the whole pipeline.*   
It also contains the **human genome assembly hg19**.  
In order to build the *hg19* **bowtie2** indexes, type the following line:
```
bowtie2-build /path/to/MetaShot_reference_data/Homo_sapiens/hg19.fa.gz hg19
```

Please note that `/path/to/MetaShot_reference_data` is the path were the MetaShot reference data were downloaded.  

The human micro-satellites sequences for the hg19 assembly were obtained from the [**UCSC Genome Browser**](http://genome.ucsc.edu/) by using the table browser tool. The list of human repeats was retrieved from the [**GIRI (Genetic Information Research Institute) Repbase**](https://www.girinst.org).
In order to improve the analysis reproducibility both the files are available **HERE**.  

## Bioinformatic Workflow


For reproducibility purposes, sequencing data were deposited as raw reads. Nonetheless, considering that the most intense and computational expensive steps of the bioinformatic workflow were the taxonomic assignment and the metagenome assembly performed by MetaShot and metaSPAdes, respectively, and in order to facilitate the analysis reproducibility by avoiding to repeat one or both these steps, the **unassigned PE reads** and the **assembled scaffolds** are available in a [*Zenodo*]() repository:  
* [MetaShot unassigned read]());  
* [Scaffolds]()).  


### 1. Taxonomic assignment of Illumina PE reads by exploiting MetaShot
The whole MetaShot workflow was performed by typing the following command:  
```
MetaShot_Master_script.py -m read_list.tsv -p parameters_file
```
In particular:
- `-m` refers to a *tsv* file containing a list of PE reads files;  
- `-p` refers to a structured file containing all the info MetaShot needs to perform the analysis. For more info, have a look at the *MetaShot setting up* section in the [README](https://github.com/bfosso/MetaShot#metashot-setting-up) file.   

It produces several files and folders but the most important are:  
* *ambiguos_pe_read.lst*: textual fine containing all the ambiguous PE reads (reads mapping on more than one reference division);  
* *bacteria_CSV_result.csv*: tabular file summarising the taxonomic assignments for prokaryotes;  
* *fungi_CSV_result.csv*: tabular file summarising the taxonomic assignments for fungi;  
* *protist_CSV_result.csv*: tabular file summarising the taxonomic assignments for protist;  
* *virus_CSV_result.csv*: tabular file summarising the taxonomic assignments for viruses.  

A more extensive description about MetaShot results is available [here](https://github.com/bfosso/MetaShot#result-files-interpretation).  

Following, the unassigned reads were extracted by using the `PE_extraction.py` script:

```
PE_extraction.py -u
```
This script generates a folder containing all the unassigned reads. The folder name is automatically generated by using the following format `taxid_PE_file_DAY_MONTH_YEAET_HOUR_MINUTES_SECONDS` in order to avoid data overwriting.


### 2. Meta-assembly of unassigned reads
mappatura sul genoma umano

```
bowtie2 -1 <R1files> -2 <R2files> -x /home/mchiara/hg19 --very-sensitive-local -p 12 -S <name.sam> --un-conc <name\_nonhumanPE>
```
--un-conc scrive 2 file fastq con i read non mappati

per le rnaseq, con comando analogo, ho mappato i file ottenuti al punto 1 sul trascrittoma
```
bowtie2 -1 <name\_nonhumanPE.1> -2 <name\_nonhumanPE.2> -x /home/mchiara/refseq --very-sensitive-local -p 12 -S <name.sam> --un-conc <name\_RNA_nonhumanPE>
```
i metacontig sono stati assemblati con il seguente comando (per i dati di RNAseq i file di input sono quelli generati al punto 2, per i campioni metagenomici, quelli al punto 1):
```
SPAdes-3.12.0-Linux/bin/metaspades.py       -1 <name\_nonhumanPE.1>  -2 <name\_nonhumanPE.2>   -t      12      -k      21,33,55,77,99  -o   </home/mchiara/metassemblies/name_meta>
```
gli assemblaggi sono stati mascherati prima con repeatMasker:
```
RepeatMasker -species human <name_meta> che genera un file con suffisso .masked
```
poi con windowmasker (2 comandi, 1 per calcolare le occorenze delle parole, 1 per mascherare)
```
windowmasker -in <name_meta.masked> -mk_counts > <name_meta.MK>

windowmasker -in <name_meta.masked> -ustat <name_meta.MK> -dust T -outfmt fasta > <name_meta.double-masked.fasta>
```
infine ho tolto le sequenze con più di tot N con uno script molto cretino in perl, che posso caricare su github se necessario
```
perl filter.pl <name_meta.double-masked.fasta> > <name_meta.BLAST.fasta>
```
e per il blast ho usato blastn dal pacchetto blast+ con l'opzione "remote"
```
blastn -remote -query  <name_meta.BLAST.fasta> -db nr > name_meta_BLAST.res
```
ultimo passaggio, per assegnare i metacontigs a una specie, dall'output di blast, il comando è come segue:
```
perl simple.parse.blast.pl G <name_meta_BLAST.res>
```
trovi i miei 2 script Perl in questo repo su github https://github.com/matteo14c/Passaro_et_al
