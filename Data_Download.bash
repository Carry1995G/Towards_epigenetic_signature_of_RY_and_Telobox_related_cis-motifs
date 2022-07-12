#!/bin/bash

###Author: Carina Gude
##Download genomic and chromatin data - 


#new folders and path links
mkdir -p Downloads
DOWNLOADS=/Downloads
cd ${DOWNLOADS}

####Obtain genomic data####

wget LINK

#if files are stored as .gz, unzip with:
gzip *.gz


#if files are stored as .zip files, unzip with 
unzip *.zip


###Obtain raw SRA data### 
