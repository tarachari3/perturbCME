#!/bin/bash
#Get all fastqs
#See https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit for updated sra tools

fasterq-dump --include-technical --split-files SRR9833038
