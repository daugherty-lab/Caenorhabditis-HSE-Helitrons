## Usage
### 1. Create working code directory in TSCC
```mkdir /home/tscc-user/code/helitron.hses/```
### 2. go to the UCSC assembly listing, https://genome.ucsc.edu/FAQ/FAQreleases.html and pull the species names associated with the most recent assemblies and paste into textfile named UCSC_planfile.txt (Example: hg38 is the most recent human genome release, 11-25-2020)
### Needs to be revised to add to automated protocol. Try checking index files within UCSC's FTP server for an updated listing.
```touch /home/tscc-user/code/helitron.hses/UCSC_planfile.txt```

```vim /home/tscc-user/code/helitron.hses/UCSC_planfile.txt```
### 3. This version currently needs to be adjusted to allow for auto-splitting the planfile to prevent congestion of TSCC. If traffic is low, proceed by submitting the qsub script that (1) queries each species subject in the planfile, (2) downloads rmsk and genome.2bit files from UCSC, (3) converts genome.2bit to genome.fa, (4) scans genome.fa for HSE sites, (5) deletes .2bit and .fa files to minimize space usage. HSE motif scanning can be done using MEME or regex, but regex is the published approach in the manuscript.
### Consider converting this to "Fetch available rmsk files" only, and combine with submit_motifsearch.sh, submit_chromsize.sh"
```chmod u+x submit_rmsk_ucsc_default.sh```
```bash submit_rmsk_ucsc_default.sh```
### 4. Run submit_hse_helitron_counter.sh to generate summaries of HSEs within helitrons -- Needs to be checked for functionality
```chmod u+x submit_hse_helitron_counter.sh```

```bash submit_hse_helitron_counter.sh``` 
### 5. Run submit_hsepileup.sh to generate listing of query genomes and the associated counts of HSE peaks (by size)
```chmod u+x submit_hsepileup.sh```

```bash submit_hsepileup.sh``` 
### 6. Visualize overlap of HSEs and Helitrons across phylogenies
### 7. Run pybedtools to test for signifance of overlap using fisher exact test
### Extra: Tajima's D on C. briggs VCF.