## Introduction

DisCVR is a viral detection tool which allows the identification of known human viruses in clinical samples from high-throughput sequencing (HTS) data. It uses the k-mers approach in which the sample reads are decomposed into k-mers and then matched against a virus k-mers database. The built-in database is a list consists of all k-mers (k=22) that are not low-complexity and only found in the viral genomes but not in the human genome. Each k-mer is assigned the taxonomic label of all viral genomes that contain that k-mer in the NCBI taxonomy tree. These assignments are made at the species and strains taxonomic level.

DisCVR has a user-friendly Graphical User Interface (GUI) which runs the analysis on the sample and shows the results interactively. It enables the visualisation of the coverage of the virus genomes found in the sample in order to validate the significance of the results. In addition, DisCVR is a generic tool which can be used with other non-human viruses by facilitating the build and use of customised k- mers database.

DisCVR is designed to run on machines with low processing capacity and small memory.

## Getting started

[System Requirements](#system-requirements)  
[Installation](#installation)  
[DisCVR Built-in Databases](#discvr-builtin-db)  
[DisCVR Classification](#discvr-classification)  
[Classification Output](#classification-output)  


## [System Requirements](#system-requirements)

1. Disk Space: DisCVR.jar requires ~700MB space for installation with the built-in databases. It is recommended to have space for 2x sample size when using DisCVR for classification as the process involves writing temporary files to disk. When building a customised database, the amount of space depends on the size of the viral sequences and the k size. For example, extracting k-mers of size 32 from the human genomes generates a file that is 80GB in size. If the viral data sequences are 3GB in size, then the minimum disk space needed to build a customised database is 200GB.

2. Memory: DisCVR runs efficiently on a machine with 4GB RAM, which is the current standard for PCs. It is much faster on machines with higher RAM such as 8GB RAM. However, the amount of RAM depends on the number of the sequences used and the size of the k-mer. The larger the dataset and the k-mer, the larger the amount of RAM needed. Therefore, in the case of "out of memory" errors, the Java heap space should be increased.

## [Installation](#installation)

1. Operating System: DisCVR runs on both Windows and Linux platforms. To use DisCVR, the users need first to download the appropriate folder for their operating system.

2. Java: Java (1.8 or above) must be installed and the full path to the jre\bin folder should be included in the system variables. Java can be downloaded from: http://www.oracle.com/technetwork/java/javase/downloads/jre8- downloads-2133155.html 

3. DisCVR.jar: After downloading the DisCVR zipped folder, it is recommended to use a tool, such as 7-zip, to unzip the Windows OS version and extract all files to a local directory. For Linux and Mac version, open a command prompt and move to the location of the zipped folder. Type the following commands to unzip the folder:
   `tar -xzvf DisCVR_Linux.tar.gz`
This creates a folder, called DisCVR. The contents of DisCVR consists of one jar: DisCVR.jar and a lib folder which are used to run the classification. The script file: downloadDataAndRefSeq.sh and the folders: bin, customisedDB, and TestData which are needed to build a customised database.  
**IMPORTANT:** The full path to DisCVR directory must NOT contain space nor the dot "." to avoid conflict with the files naming during the classification process.

4. Dependencies: DisCVR uses external libraries such as KAnalyze, for k-mers counting, and JFreechart packages, for graphs plotting. It makes use of Tanoti, a Blast-based tool for reference assembly. These are 10 files in total and they are in the lib folder. It is important not to alter the lib folder or its contents and to ensure that it is in the same path as the jar file.

5. If you want to build a customised database, the following NCBI tools and files must be downloaded and installed:

   The NCBI eutilities tools are used to download data. The tools can be found at: ([ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/](ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/)). The full path to the edirect folder should be added to the system variables  

   The NCBI taxdump files are used for taxonomy information retrieval when building a customised database. The file can be downloaded from ([ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)). The file taxdump.tar.gz should be downloaded and unzipped. The two files: names.dmp and nodes.dmp MUST be copied to the customised database folder: customisedDB which is in the same path as DisCVR.jar.

6. To test if the tools are installed properly, open a command prompt and type the following:

   To know what Java version is installed: `java –version`  
   This should state `java version 1.8.0_<some number>`

   To see if Java is added to the path: `java`  
   If the `jre\bin` is not added to the path, you will see the following message: `java is not recognized as an internal or external command, operable program or batch file`.

   To see if eutilities tools is added to the path: `esearch`
   This should state `"Must supply -db database on command line`"

## [DisCVR Built-in Databases](#discvr-builtin-db)

DisCVR tool has three virus k-mers (k=22) databases which are built in DisCVR.jar:

1. Human Hemorrhagic virus dataset (hemorrhagic dataset)

2. Human Respiratory virus dataset (respiratory dataset)

3. Human Pathogenic viruses dataset (pathogenic dataset)

The hemorrhagic and respiratory datasets consist of a list of selected circulating viruses. 
The pathogenic datatset consist of a list of all human pathogenic viruses identified as biological agents from by the Health and Safety Executive (HSE) in the UK available at: [http://www.hse.gov.uk/pubns/misc208.pdf](http://www.hse.gov.uk/pubns/misc208.pdf)

Both the hemorrhagic and respiratory datasets overlap with the HSE database but are smaller in size and therefore are faster to use in classification.
Moreover, each dataset has a list of viral reference genomes that is associated with the viruses in the database. The reference genomes are used in the validation stage of DisCVR.
Information about the viruses in the databases can be obtained from the Database icon on the menu bar from DisCVR GUI. The information includes the virus rank on the taxonomy tree, 
the number of sequences used in the database build to represent it, and the accession number of its reference genome. Segmented viruses have their reference genome accession number listed 
in decreasing order of the segment size. The accession numbers of the reference genomes are linked to their webpages on the NCBI which can be viewed by clicking on the reference genome accession.

## [DisCVR Classification](#discvr-classification)

The DisCVR GUI can be used to carry out a single sample classification. To launch the DisCVR GUI, either double click on DisCVR.jar, or open a command prompt and type the following commands:

```cd full/path/to/DisCVR folder
java –jar full/path/to/DisCVR.jar```


Once DisCVR is launched, there are 4 panels in the graphical interface:

1. Progress Bar: shows the progression of the classification process.

2. Input Sources: allows the user to choose the sample file and a k-mers database for classification.

3. Classification Results: shows the classification output in three sub-panels.

4. Progress Information: updates the user about classification output.

To start a classification process, the user needs to first select the sample file they wish to 
investigate. Currently, DisCVR supports .fasta, .fa, .fastq, and .fq formats. It allows for the selection 
of one file hence, for paired-end files, it is recommended to concatenate both files. The second step is to 
choose the database for the classification. The default setting is to select one of the three built-in databases. 
This sets the k-mer size to 22 and removes all sample k-mers with entropy ≤ 2.5. On the other hand, if the user 
chooses to use a customised database for the classification, then the box next to the **Customised Database** label 
must be checked. The fields for this selection are then activated for the user to upload the customised database 
file which MUST be in the customisedDB folder that is in the same path as DisCVR.jar. The k-mer size for the 
classification is extracted from the database name and the threshold to filter out low-entropy k-mers from the 
sample should be entered in its corresponding field. The entropy threshold should be the same value used in the 
build of the customised database which is recommended to be in the range of [0, 3.0]. If the users choose not to 
specify a value for the entropy threshold then the default value, 2.5, is displayed in the field and used in the 
classification.

Once, the database library is selected, the user clicks the classify button to start the classification process.

There are 4 stages in the classification process and the progress bar is updated after the completion of each 
stage. Messages in relation to the output of each stage are displayed in the progress information panel. 
This includes the number of reads in the sample, the number of distinct k-mers and their total counts in 
the database and the sample file, the number of sample k-mers after removing single copies and low entropy 
k-mers, and the number of classified k-mers and their total counts. The final message shows the number of 
viruses with classified k-mers found in the sample and the time taken to finish the classification process. 
Figure 1 shows screenshots of the sample classification process using the DisCVR GUI.

![Figure 1a](https://josephhughes.github.io/DisCVR/img/img1.png "sample file is uploaded")
*Figure1a: Screenshots of a sample classification process using DisCVR GUI. First the sample file is uploaded*

![Figure 1b](https://josephhughes.github.io/DisCVR/img/img2.png "the k-mers database is then selected from the DisCVR database library")
*Figure1b: Screenshots of a sample classification process using DisCVR GUI. The k-mers database is then selected from the DisCVR database library*

![Figure 1c](https://josephhughes.github.io/DisCVR/img/img3.png "the k-mers database is selected from the user customised database")
*Figure1b: Screenshots of a sample classification process using DisCVR GUI. The k-mers database is then selected from the user customised database*

![Figure 1d](https://josephhughes.github.io/DisCVR/img/img4.png "the classification results are displayed")
*Figure1b: Screenshots of a sample classification process using DisCVR GUI. The classification results are displayed*


## [Classification Output](#classification-output)

At the end of the classification process, the progress information panel states the number of viruses with matched k-mers to the database. 
Detailed information about the classification results are displayed on the centre panel.

The classification results panel consists of three sub-panels: Scoring, Summary, and Full Analysis. Once the classification is completed, a bar 
chart showing up to three viruses is displayed on the scoring panel. These are the viruses with the most number of classified distinct k-mers 
found in the sample. The chart shows for each virus, the number of distinct k-mers that are specific to the virus as well as the number of shared 
k-mers with other viruses in the sample, which are referred to as non-specific k-mers. The Summary panel gives a list of the viruses found in the 
sample along with their number of specific and non-specific distinct k-mers. The Full Analysis panel shows a table with taxonomic and detailed 
information about the viruses with classified k-mers from the sample. The results are uploaded on the table such that the first row shows the virus 
with the highest number of k-mers found in the sample. However, the users can click on any column heading to sort out the results in the table 
according to the information in the column. The full analysis table consists of 8 columns:

1. Virus Name: the scientific name for the virus taking from the NCBI names.dmp file.

2. Taxa ID: the taxonomy identification for the virus in the NCBI taxonomy tree.

3. Virus Rank: the rank of the virus according to the NCBI nodes.dmp file.

4. Total counts of k-mers in DB: the total counts of k-mers that represent the virus in DisCVR’s
database.

5. No. of distinct Classified k-mers: the number of distinct k-mers that represent the virus in the
sample after removing single copies and low entropy k-mers and matched with the k-mers
database.

6. (%) of distinct Classified k-mers: the percentage of distinct classified k-mers that represent the
virus in the sample.

7. Total counts of Classified k-mers: the total number of k-mers that represent the virus in the
sample; some distinct k-mers can occur in the reads multiple times.

8. (%) of total classified k-mers: the percentage of the total number of k-mers that represent the
virus in the sample’s total number of k-mers.

The table can be saved as .csv file from the File icon on the GUI tool bar. Figure 2 shows screenshots
of the classification results sub-panels.

