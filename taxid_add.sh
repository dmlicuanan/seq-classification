# change directory
cd /mnt/d/Documents/NGS/entrez/
cd /mnt/d/Documents/NGS/entrez/taxonomy

# create directory for taxonomy 
mkdir taxonomy
cd taxonomy/

# download taxonomy dump
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

# extract all files from archive.tar
tar -xf taxdump.tar.gz

# list of species that need taxids are in:
head /mnt/d/Documents/NGS/entrez/bold_fastas/bold_id_list
# extract list 
cut -d '|' -f 2 /mnt/d/Documents/NGS/entrez/bold_fastas/bold_id_list | sort -u




# COI compilation
# create directory for where kraken fastas will be output
mkdir /mnt/d/Documents/NGS/entrez/kraken_fastas

# change directory 
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb/wormstaxlist_fasta_cut_kraken
# compile NCBI fastas
cat *_ncbi.fasta >> /mnt/d/Documents/NGS/entrez/kraken_fastas/ncbi_coi.fa
sta

# change directory
cd /mnt/d/Documents/NGS/entrez/bold_fastas/bold_fastas_kraken/
# compile BOLD fastas
cat *_bold.fasta >> /mnt/d/Documents/NGS/entrez/kraken_fastas/bold_coi.fasta

# change directory
cd /mnt/d/Documents/NGS/entrez/kraken_fastas/
# combine NCBI and bold fastas
cat bold_coi.fasta ncbi_coi.fasta >> worms_coi.fasta

# create fasta of combining worms fasta with subset of mares fasta
cat worms_coi_deduplicated.fasta mares_no_bar_subset.fasta >> worms_mares_coi.fasta





# 16S compilation
# change directory 
cd /mnt/d/Documents/NGS/entrez/wormstaxlist_gb_16S/fasta_cut_kraken
# compile NCBI fastas
cat *_ncbi.fasta >> /mnt/d/Documents/NGS/entrez/kraken_fastas/worms_mares_16S.fasta
