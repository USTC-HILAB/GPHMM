GPHMM:Parameter Hidden Markov Model
==================================
##Introduction
GPHMM is a novel statistical method dedicated to identify copy number alteration and loss of heterozygosity (LOH) in tumor samples using whole genome SNP arrays. In contrast to other HMM methods, a distinct feature of GPHMM is that the issues reported to tangle genotyping signals, such as baseline shift of LRR signal due to aneuploidy, normal cell contamination and genomic waves, are quantitatively modeled by global parameters and integrated in the whole statistical framework. Therefore GPHMM provides not only optimal solutions to these issues but also automatic and accurate identification of copy number and LOH status of each SNP in the array.<br>

---------------
##Recent Update
2015/05/01: GPHMM 1.4 is released. The source codes are available now (details see Download section). Also, the visualization is improved for better illustration.
- 2012/10/30: An important statement is issued for the controversial performance evaluation in a published paper.
- 2012/10/05: GPHMM version 1.3 introduced an advanced strategy to determine optimal initial global parameters. This new version can accurately estimate aberrations in tetraploid samples (with average copy number (ACN) >=4.0), which might be occasionally identified as near-diploid by previous versions. Also, GPHMMv1.3 will generate a log file (LOG.txt in present directory), which gathers all important parameters of samples.
-2012/06/21: GPHMM version 1.2.2 is released with minor bugs (error in loading predefined tumor cell proportion ranges from config file) fixed.
- 2012/03/01: GPHMM version 1.2.1 is released with minor bugs fixed.
- 2011/11/09: GPHMM version 1.2 is released with maximum copy number up to 7. Also, an improved approach based on Newton-Raphson Method is adopted in GPHMMv1.2 for better estimation of cancer cell proportion. Besides, this new version can handle a configuration list file(provided by user) which includes predefined ranges of cancer cells admixed in tumor samples.
- 2011/11/01: User contributed program(provided by Geert Vandeweyer at University of Antwerp): compiled standalone binaries for ubuntu(10.04 LTS, 64bit) users. Usage and Package(Original / Mirror).
- 2011/07/28: GPHMM version 1.1 is released with minor bugs fixed and copy number/allelic imbalance information added in result files. Also Mex-functions for different operating systems (Windows, Linux, Mac) are now provided in this new version. Download and usage page of GPHMMv1.1 from here.

-------------------
##Reference
Li, A., Liu, Z., Lezon-Geyda, K., Sarkar, S., Lannin, D., Schulz, V., Krop, I., Winer E., Harris, L., Tuck, D. GPHMM: A novel method for identification of copy number alteration and loss of heterozygosity in complex tumor samples using whole genome SNP arrays. Nucleic Acids Research, 2011 39(12):4928-41.Link
The supplementary file of this manuscript can be downloaded from here.

-------------------
##Download and Installation
The GPHMM tool is implemented with MATLAB/C and is freely available.<br>
\#NEW: the source codes of GPHMM1.4 are available now, user can download from here. Details see the "example.m" after unzipping the software package.<br>
Users can also download the GPHMM package(v1.3) from here, which requires pre-installed MATLAB environment. For installation, extract the ZIP file to a directory and make sure the path of the extracted files is saved in MATLAB.<br>
We also provide an alternative solution by which users can deploy this tool to a computer without MATLAB (currently it only works for Windows platform). To do so, first download MATLAB Component Runtime library utility from here(extracted from MATLAB2010b 64bit) and install it on the deployment machine (Detailed information is available at the MathWorks website). Second, download the package(GPHMMv1.3) from here(built by MATLAB2010b 64bit, in console application mode) and extract the ZIP file to a directory, then use the command line in Windows to run GPHMM as a stand-alone program.<br>
Test sample: a sample of SNP array data for testing is also available by clicking here.
System requirement: 2G or more memories are required to properly run GPHMM and in case out-of-memory issue is reported by MATLAB, try to start it without loading the Java virtual machine (JVM). More information is available at the MathWorks website.
Additional information: population frequency of B allele (PFB) and GC content of each SNP in the array are required to run this utility. The two data files (‘hhall.hg18_m.pfb’ and ‘hhall.hg18.gcmodel’) in GPHMM that include such information are obtained from PennCNV package (Courtesy of Dr. Kai Wang).

-----------------------
##Quick Examples
####(1) GPHMM (file "GPHMM.p" or "GPHMM.exe")
Description:the original function implementing the proposed algorithm
Function:GPHMM (SNP_Array_Path, Pfb_File_Path, GC_File_Path, Result_Path, Genotype_Path, Sameformat)
Arguments:
SNP_Array_Path:   Path of the directory containing SNP array data files.
Pfb_File_Path:       Path of the PFB file.
GC_File_Path:       Path of the GC file.
Result_Path:          Path of the directory containing GPHMM results.
Genotype_Path:   Path of the directory containing predicted genotypes of all SNPs. If it is an empty string (default), no genotype prediction will be made.
Sameformat:           If true (1), for all SNP data files GPHMM will parse PFB and GC files and match SNP ids only once to save time. Default value is false (0).
configlist:             (New in GPHMMv1.2) Path of the configure list file. It contains the names of the data files which will be analyzed by GPHMM, and the tumor cell proportion range. if provided, GPHMM will only analyze the data files from list with preset tumor proportion range; if not, GPHMM will analyze all the files in the SNPdata directory(SNP_Array_Path) with default tumor proportion range(0.1 to 1.0).

--------------------------
##Examples of usage:
Suppose GPHMM is extracted to ‘C:\GPHMM’ and SNP array data is located in 'C:\SNPdata'. To use GPHMM, first create a directory in which the result files will be saved, for example ‘C:\results’. If genotype prediction is also required, create another directory for the genotype files, for example ‘C:\genotypes’. Then use following commands in MATLAB to run GPHMM. (GPHMM will not automatically create these directories for you, so make sure you have successfully created them and all the paths you inputted are correct)
- Case 1. Different SNP array data formats, no genotype prediction
GPHMM( 'C:\SNPdata', 'C:\GPHMM\hhall.hg18_m.pfb', 'C:\GPHMM\hhall.hg18.gcmodel', 'C:\results', '', 0,'config.txt'), or
GPHMM( 'C:\SNPdata', 'C:\GPHMM\hhall.hg18_m.pfb', 'C:\GPHMM\hhall.hg18.gcmodel', 'C:\results', '', 0)
- Case 2. Different SNP array data formats, genotype prediction needed
GPHMM( 'C:\SNPdata', 'C:\GPHMM\hhall.hg18_m.pfb', 'C:\GPHMM\hhall.hg18.gcmodel', 'C:\results', 'C:\genotypes', 0 ,'C:\results', '', 0 ,'config.txt'), or
GPHMM( 'C:\SNPdata', 'C:\GPHMM\hhall.hg18_m.pfb', 'C:\GPHMM\hhall.hg18.gcmodel', 'C:\results', 'C:\genotypes', 0)
- Case 3. Same SNP array data format and genotype prediction needed
GPHMM( 'C:\SNPdata', 'C:\GPHMM\hhall.hg18_m.pfb', 'C:\GPHMM\hhall.hg18.gcmodel', 'C:\results', 'C:\genotypes', 1 ,'config.txt'), or
GPHMM( 'C:\SNPdata', 'C:\GPHMM\hhall.hg18_m.pfb', 'C:\GPHMM\hhall.hg18.gcmodel', 'C:\results', 'C:\genotypes', 1)
- Case 4. To run GPHMM as a stand-alone program for Case 1. In the command line, type
GPHMM.exe C:\SNPdata C:\hhall.hg18_m.pfb C:\hhall.hg18.gcmodel C:\results NA 0 config.txt, or
GPHMM.exe C:\SNPdata C:\hhall.hg18_m.pfb C:\hhall.hg18.gcmodel C:\results NA 0
- Case 5. To run GPHMM as a stand-alone program for Case 2. In the command line, type
GPHMM.exe C:\SNPdata C:\hhall.hg18_m.pfb C:\hhall.hg18.gcmodel C:\results C:\genotypes 0 config.txt or
GPHMM.exe C:\SNPdata C:\hhall.hg18_m.pfb C:\hhall.hg18.gcmodel C:\results C:\genotypes 0

-------------------------
##Input Files
Requirement of data format: (Important: GPHMM will NOT work properly if the formats of input files are not correct.)
####I. For SNP array data:
Two formats can be used in a SNP array data file
a) Data format 1. A SNP data file with data format 1 has 5 columns with one headline specifying SNP identifier, chromosome identifier, chromosomal position, BAF value, LRR value. An example is shown below:
```
Name   Chr   Position   BAF   LRR
rs12354060   1   10004   1.00   0.00
rs2691310     1   46844   0.00   0.40
...
```
b) Data format 2 (output format of tQN). Normalized SNP array data files created by tQN* can be used directly as the input files of GPHMM, which have 7 columns with one headline specifying SNP identifier, chromosome identifier, chromosomal position, BAF value, LRR value, X and Y values after normalization. An example is shown below:
```
Name   Chr   Position   BAF   LRR   tQN X   tQN Y
rs12354060   1   10004   1.00   0.00   0.08   1.85
rs2691310     1   46844   0.00   0.40   1.04   0.78
...
```
Also, each row in either format should be sorted by the position (physical coordinate) on the corresponding chromosome in ascending order. Rows of X, Y and mitochondrial chromosomes should always come after rows of autosome (Chr1-22). Chromosome identifiers for X, Y and mitochondrial chromosomes can be either numbers or one character such as ‘X’, ‘Y’ and ‘M’.
* Please Note : It has been reported that bias between the two dyes used in the Illumina assay can lead to unsymmetric BAF bands and noisy data. Therefore, an efficient tool (tQN) for removing such effect is highly recommended for data pre-processing before using GPHMM.
####II. For PFB data:
Use PFB file ‘hhall.hg18_m.pfb’ provided in the package.
####III. For GC data:
Use GC file ‘hhall.hg18.gcmodel’ provided in the package.
####IIII. For configure list file (New in GPHMMv1.2) :
The list file contains 3 columns, denoting the data file name, the minimum tumor cell proportion(>=0.1) and maximum tumor cell proportion(<=1.0) respectively. An example is shown below:
```
sample1.txt	0.1	1.0
sample2.txt	0.3	0.5
sample3.txt	0.5	0.9
...
```

-----------------------------
Output Files
####I. Result file:
a) The first section in a GPHMM output file is a summary of the global parameters estimated from SNP array data:
```
             Summary of GPHMM results (version 1.2)          
General information of this cancer sample:                      
   LRR correction factor: -0.1614
   GC coefficient: -0.0140
   Proportion of abnormal cells in the sample: 0.5722
   Standard deviation of LRR signal: 0.1957
   Standard deviation of BAF signal: 0.0679
   Proportion of all abnormal chromosomal regions: 0.8125
   Estimated average cancer DNA index: 1.4225
```
b) In the second section, segmented chromosomal regions are listed with corresponding hidden states. The header is formatted as below:
```
Chr: Chromosome number
StartPos:      Starting chromosomal coordinate
EndPos:       Ending chromosomal coordinate
State:           Predicted state as defined in GPHMM
Score:         Averaged posterior probabilities of the SNPs in current segment
CN:              Copy number of current segment (0-5)
AI:                Allelic imbalance status of current segment (1: deletion; 2: LOH; 3: heterozygosity)
Length:        Length of current segment
StartSNPid:  The ID of the first SNP in current segment
EndSNPid:   The ID of the last SNP in current segment
StartIndx:     The index of the first SNP in current segment
EndIndx:       The index of the last SNP in current segment
```
An example of the second section is shown below:
```
Chr   StartPos   EndPos   State   Score   CN   AI   Length   StartSNPid   EndSNPid   StartIndx   EndIndx
1   995669     19052411   8   0.9937   4   3   2117   rs3934834   rs4920566   1      2117
1   19088005   19929346   7   0.8757   4   3   98     rs4912075   rs4326647   2118   2215
...
```
####II. Genotype file:
The format of the genotype file is as follow:
```
Name   Chr   Position   CN   B Allele CN   Genotype
rs12354060   1   10004   2   0   AA
rs2691310     1   46844   2   1   AB
rs2531266     1   59415   2   2   BB
...
```

-------------------
##Plotting Results
Function:
GPHMM_plot (SNP_Array_Path, Result_Path, Plot_Path)
Arguments:
SNP_Array_Path:    Path of the directory containing SNP array data files.
Result_Path:           Path of the directory containing GPHMM results.
Plot_Path:               Path of the directory containing the plots of every tumor sample
configlist:                (New in GPHMMv1.2) Path of the configure list file.
Examples of usage:
Suppose GPHMM has successfully analyzed the SNP array data located in ‘C:\SNPdata’(or by the sequence in 'config.txt'), and the result files are stored in ‘C:\results’. To plot the results, first create a directory for image files, for example ‘C:\plots’. Then use the following command in MATLAB to run GPHMM_Plot. (make sure you have successfully created it and the paths you inputted are correct)
GPHMM_Plot ( 'C:\SNPdata', 'C:\results', 'C:\plots', 'config.txt' ) or
GPHMM_Plot ( 'C:\SNPdata', 'C:\results', 'C:\plots' )

-----------------------------------
##Contact
Ao Li, Associate Professor, Department of Electronic Science and Technology, USTC (Aoli at ustc dot edu dot cn )
Yuanning Liu, Graduate Student, Department of Electronic Science and Technology, USTC (lynn100 at mail dot ustc dot edu dot cn )
