# SigMA
1. [ Description ](#desc)
2. [ Manual ](#man)
3. [ Functionalities ](#func)
4. [ Quick start ](#quick)
5. [ Quick start shiny ](#shiny)
6. [ Input file format ](#input)
6. [ Tumor type tags ](#tissue)

<a name="desc"></a>

### 1. Description

SigMA is a signature analysis tool optimized to detect the mutational signature associated to HR defect, Signature 3, from hybrid capture panels, exomes and whole genome sequencing. For panels with low SNV counts, conventional signature analysis tools do not perform well while the novel approach of SigMA allows it to detect Signature 3-positive tumors with 74% sensitivity at 10% false positive rate. One novelty of SigMA is a likelihood based matching: We associate a new patient's mutational spectrum to subtypes of tumors according to their signature composition. The subtypes of tumors are defined using the WGS data from ICGC and TCGA consortia, by a clustering of signature fractions with hierarchical clustering. The likelihood of the sample to belong to each tumor subtype is calculated, and the likelihood of Signature 3 is the sum of the likelihoods of all Signature 3-positive tumor subtypes. The second novel step is the multivariate analysis with gradient boosting machines, which allows us to obtain a final score for presence of Signature-3 combining likelihood with cosine similarity and exposure of Signature 3 obtained with non-negativel-least-squares (NNLS) algorithm. The multivariate analysis allows us to automatically handle different sequencing platforms. For different platforms different methods for signature analysis become more efficient, e.g. for WGS data it is not necessary to associate the tumor to a subtype of tumors, because it is possible to determine Signature 3 with NNLS acurately. We have a new feature also for these cases and we calculate the likelihood of NNLS decomposition to be unique. This likelihood value was found to be the most influential feature in the multivariate analysis.

<a name="man"></a>

### 2. Manual

[ manual_sigma.pdf ](manual_sigma.pdf)

<a name="func"></a>

### 3. Functionalities:

* Optional input format for mutation data vcf, maf, or directly a 96-dimensional matrix
* If vcf or maf files are provided it constructs the 96-dimensional SNV matrix with trinucleotide context
* Contains a database of signature clusters, defining tumor subtypes, for multiple tumor type
* For each input sample calculates:
   * Likelihood values of the tumor subtypes
   * Cosine similarity for all signatures
   * Exposures of signatures which are found in the matching tumor type
* Output
  * Tumor type specific final decision for existence of Signature 3
  * Plot of the mutational spectra and other relevant variables obtained from this such as likelihoods, cosine similarity and exposures, as well as the information on the signatures of the matched tumor subtype

<a name="quick"></a>

### 4. Quick start 
`git clone https://github.com/parklab/SigMA`

`cd SigMA`

`Rscript install.R`

`Rscript test.R`

More details on the optional settings can be found in the User's Manual

<a name="shiny"></a>

### 5. Quick start shiny 
`git clone https://github.com/parklab/SigMA`

`cd SigMA`

`git checkout shiny`

`Rscript install.R`

`Rscript runShiny.R`

<a name="input"></a>

### 6. Input file format

MAF file example: For MAF file format, both a single file with multiple samples, or a directory which contains multiple files with a single sample in each, is acceptable
You can download the test_mutations_50sample.maf file on this repository and use this for testing. 
VCF file example: Direct SigMA to a directory with multiple VCF files. Two test samples can be found in:
`cd inst/extdata/`

<a name="tissue"></a>

### 7. Tumor type tags

Tags and descriptions of performance of SigMA for signature 3 detection

| Tumor type                           | Tag           | Details     | Platform        |
| ------------------------------------ | ------------- | ----------- | --------------- |
| Ewing Sarcoma                        | bone_other    | Tested      | Exome/WGS       |
| Urothelial Bladder Cancer            | bladder       | Exploratory | Panel/Exome/WGS |
| Breast Cancer                        | breast        | Tested      | Panel/Exome/WGS |
| Colorectal Adenocarcinoma            | crc           | Exploratory | Panel/Exome/WGS |
| Oesophageal Carcinoma                | eso           | Tested      | Panel/Exome/WGS |
| Glioblastoma                         | gbm           | Exploratory | Panel/Exome/WGS |
| Lung Adenocarcinoma                  | lung          | Exploratory | Panel/Exome/WGS |
| Lymphoma                             | lymph         | Exploratory | Panel/Exome/WGS |
| Medulloblastoma                      | medullo       | Tested      | Exome/WGS       |
| Osteosarcoma                         | osteo         | Tested      | Panel/Exome/WGS |
| Ovarian Cancer                       | ovary         | Tested      | Panel/Exome/WGS |
| Pancreas Adenocarcinoma              | panc_ad       | Tested      | Panel/Exome/WGS |
| Pancreas Neuroendocrine              | panc_en       | Tested      | Panel/Exome/WGS |
| Prostate Adenocarcinoma              | prost         | Tested      | Panel/Exome/WGS |
| Stomach Adenocarcinoma               | stomach       | Tested      | Panel/Exome/WGS |
| Thyroid Cancer                       | thy           | Exploratory | Panel/Exome/WGS |
| Uterine corpus endometrial carcinoma | uterus        | Tested      | Panel/Exome/WGS |

