#!/usr/bin/env bash

set -e

snakemake -s ../Snakefile  --config database_dir='databases' input_faa=input_catalog.faa $@
