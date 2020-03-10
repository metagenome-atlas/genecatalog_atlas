#!/usr/bin/env bash

set -e

source activate genecatalog

snakemake --configfile ../config/default_config.yaml -s ../Snakefile $@  --config database_dir='databases' input_faa=input_catalog.faa 
