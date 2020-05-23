#!/usr/bin/env bash

set -e

source activate genecatalog

snakemake --configfile ../config/default_config.yaml -s ../Snakefile $@  --config database_dir='databases' input_faagz=input_catalog.faa.gz 
