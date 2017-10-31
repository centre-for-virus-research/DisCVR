# Overview
DisCVR is a computer program which allows diagnosticians to detect known human viruses in clinical samples from High-Throughput Sequencing data. It works by creating a viral database of _k_-mers (22 nucleotide sequences) from a set of known sequence is created. These _k_-mers are taxonomically lavelled according to the viruses from which they originated. Each read in the sample is then screened for the presence of _k_-mers in the viral database, and a list of all viruses with _k_-mers found in the sample is shown. 

DisCVR is a fast, accurate viral detection tool designed to analyse HTS data and validate the results interactively on comput-ers that have limited resources. At present, DisCVR is a human viral diagnostic tool, but it could be readily extended to include non-viral human pathogens and pathogens of other hosts. 
# Framework
