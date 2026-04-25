# GenBank Taxonomy Parsing Pipeline (R)

A reproducible R pipeline for parsing GenBank records and enriching BLAST results with organism-level taxonomy from NCBI.

## Features
- Queries NCBI using accession IDs
- Extracts organism, host, and country metadata
- Parses taxonomy into structured columns
- Handles API errors and rate limits
- Caches results for efficiency

## Tech Stack
- R
- rentrez
- tidyverse

## Usage
1. Place input files in `data/`
2. Run script in `src/`
3. Outputs saved to `output/`

## Motivation
Manual lookup of GenBank records is slow and cumbersome.  
This pipeline automates the process into a scalable workflow.

## Author
Nick M
