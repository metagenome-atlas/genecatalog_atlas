#!/usr/bin/env bash

set -e

source activate genecatalog

snakemake -s ../Snakefile $@  --config database_dir='databases' input_faa=input_catalog.faa
#--configfile ../config/default_config.yaml 
