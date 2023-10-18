# seq-classification

Scripts for classifying eDNA reads into operational taxonomic units.

1. Create list of marine taxa using WoRMS dataset
2. Download COI sequences from GenBank
3. Download COI sequences from BOLD and find species for which the download has to be repeated due to errors in initial run.
4. Merge COI sequences from GenBank and BOLD to create one fasta file.

#### Files and directories:
- bold_download.R
    - Uses R package bold to download COI sequences from the Barcode of Life Data Systems (BOLD)
- bold_format.R
    - Formats sequences downloaded from BOLD to format acceptable to Kraken by attaching taxids from NCBI taxdump
    - Deduplicates COI and 16S sequences
    - Check taxids represented in MARES database but absent in compiled WoRMS database
- gb_parse.R
    - Creates FASTA files from GenBank flat files
    - Parses files with multiple GenBank records, extracts seqeunce information, and the COI/16S portion of sequences if mitogenomes are provided
- kraken_vis.R
    - Creates long taxonomy table from NCBI taxdump so that species 

***

Author's note: 
Project was started on 2022-07-21 and concluded on 2023-06-29. Latest code edits: 2023-10-18.

