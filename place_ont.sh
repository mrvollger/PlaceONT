#!/usr/bin/env bash

#
# snakemake paramenters
#
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
snakefile=$DIR/snake.py
jobNum=20
waitTime=60 # this really needs to be 60 on our cluster :(
retry=1 # numer of times to retry the pipeline if it failes
# I allow a retry becuase sometimes even the really long waittime is not enough,
# and the files are actaully there

#
# QSUB parameters, these are only the defualts, they can be changed with params.sge_opts
# Allow snakemake to make directories, I think it slows things down when I done with "waitTime"
#
logDir=logs
mkdir -p $logDir
E=$logDir'/snakejob_{rule}_{wildcards}_e'
O=$logDir'/snakejob_{rule}_{wildcards}_o'
ram=4G
defaultCores=1

#
# run snakemake
#
snakemake -p \
	-s $snakefile \
	-j $(nproc)

