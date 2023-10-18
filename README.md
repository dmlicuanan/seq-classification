# seq-classification

Scripts in this repository are for assigning taxonomic labels to environmental DNA (eDNA) reads. The general workflow is as follows:

1.  Create and curate a list of marine taxa 
2.  Download COI and 16S sequences from NCBI and BOLD
3.  Build Kraken databases for each gene to query reads against
4.  Perform bioinformatic processing on eDNA reads (quality checks, merging, filtering, trimming, etc.)
5.  Classify reads using Kraken
6.  Analyze and visualize resulting classifications

#### Data visualizations produced by scripts:

| <img src="images/taxalist_coverage.png" alt="Coverage of taxon lists from WoRMS and MARES" /> | 
|:--:| 
| Coverage of taxon lists from WoRMS and MARES |

| <img src="images/kraken_resolution.png" alt="Resolution of Kraken classification" width=65% />| 
|:--:| 
| Resolution of Kraken classification |

| <img src="images/bytaxid_byfamily_nmds.png" alt="nMDS based on presence-absence of taxids (left) and families (right)"> | 
|:--:| 
| nMDS based on presence-absence of taxids (left) and families (right) |

|  <img src="images/shared_families.png" alt="Variation of families in site replicates" /> | 
|:--:| 
| Variation of families in site replicates |

| <img src="images/unique_taxa.png" alt="Phyla detected in single replicates" width=65% /> | 
|:--:| 
| Phyla detected in single replicates |

| <img src="images/site_map.png" alt="Phyla captured across sites" width=65% /> | 
|:--:| 
| Phyla captured across sites |

| <img src="images/data_yield.png" alt="Data yield for each sample" width=65% /> |
|:--:| 
| Data yield for each sample |

#### Files and directories:

-   worms_list.R
    -   Creates a list of species using a World Register of Marine Species (WoRMS) dataset
    -   Extracts taxids in MARine Eukaryote Species (MARES) database that are absent in WoRMS
    -   Creates metadata table for taxids in WoRMS and MARES databases
-   bold_download.R
    -   Uses R package bold to download COI sequences from the Barcode of Life Data Systems (BOLD)
-   bold_format.R
    -   Formats sequences downloaded from BOLD to format acceptable to Kraken -- by attaching taxids from NCBI taxdump
    -   Deduplicates COI and 16S sequences
    -   Check taxids represented in MARES database but absent in compiled WoRMS database
-   ncbi_download.sh
    -   Cleans taxon list
    -   Downloads GenBank flat files (COI and 16S) from NCBI using Entrez Direct
-   gb_parse.R
    -   Creates FASTA files from GenBank flat files
    -   Parses files with multiple GenBank records, and extracts seqeunce information and the COI/16S portion of sequences if mitogenomes are provided
-   taxid_add.sh
    -   Downloads NCBI taxdump
    -   Compiles COI and 16S FASTA files
-   nextera_adapters.R
    -   Generates FASTA file containing NGS adapters to be used for trimming reads
-   ngs_cleanup.sh
    -   Checks sequence quality of raw reads using fastQC
    -   Summarizes fastQC results with multiQC
    -   Merges, trims, and filters reads
    -   Classifies COI and 16S reads using Kraken
    -   Counts reads remaining after each bioinformatic step
-   read_counts.R
    -   Processes read counts per bioinformatic step to plot data yield
-   kraken_vis.R
    -   Creates long taxonomy table from NCBI taxdump so that taxids can be linked to higher orders of classification
    -   Checks coverage of taxa lists from WoRMS and MARES databases
    -   Checks coverage of compiled COI and 16S sequence databases to be used for Kraken classification
    -   Compiles results of Kraken classification per read
    -   Visualizes resolution of Kraken classification
    -   Performs ordinations
    -   Finds taxa that are unique, shared, common across replicates
    -   Shows phyla detected in map
-   data directory
    -    Species lists
    -    List of indices used for adapters in sequencing
    -    Primer overhangs used
-   envrironments directory
    -   Contains YAML files for conda environments used
-   images directory
    -   Contains some figures produced by the scripts

------------------------------------------------------------------------

Author's note: Project was started on 2022-07-21 and concluded on 2023-06-29. Latest code edits: 2023-10-18.
