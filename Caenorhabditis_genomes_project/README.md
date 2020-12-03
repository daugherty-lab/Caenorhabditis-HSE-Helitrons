# Caenorhabditis-HSE-Helitrons
HSE-Helitron overlap analyses by bvtsu for Garrigues et al., 2019

## This work addresses the following:
### 1. Does enrichment of HSEs inside helitrons (red fraction) extend to other Caenorhabditis species?

![HSE-Helitrons-in-CGP-genomes](CGP_black_red.png)
Published work with CGP genomes suggest that other Caenorhabditis species have this signature.

## Dependencies
### 1. RepeatMasker (via conda)
```conda config --add channels bioconda```

```conda install -c bioconda repeatmasker```

### 2. FIMO (MEME Suite package)

## Downloading genomes from Caenorhabditis Genome Projects:
### 1. Visit http://download.caenorhabditis.org/v1/sequence/
### 2. For species genomes of interest, select only the ".scaffolds.fa.gz" files.  Using the C. elegans and C. briggsae genomes as examples, download "Caenorhabditis_elegans_WBcel235.scaffolds.fa.gz" and "Caenorhabditis_briggsae_CB4.scaffolds.fa.gz" by using the checkbox on the top-left of each file name to select multiple files followed by clicking the download icon on the top-left portion of the webpage (next to the search icon). To conduct the entire analyses, simply select every single ".scaffolds.fa.gz" file to download. (Bonus point challenge: Figure out a way to wget these files and assign it to me in a task)
### 3. Move all locally downloaded files to a "genome" folder in your TSCC account using scp. Enter your university Active Directory password when prompted.
```scp ~/Folder-With-Genomes/*.scaffolds.fa.gz username@tscc-login.sdsc.edu:~/genomes```
### 4. Log into TSCC and enter your university Active Directory password when prompted.
```ssh username@tscc-login.sdsc.edu```
### 5. Create a working folder/ git clone the repository of interest
```mkdir ~/helitron.hses```
or
```git clone https://github.com/daugherty-lab/Caenorhabditis-HSE-Helitrons/ subproject/```
#You will be prompted for your github credentials, as this is a private repo
```cd subproject```

```git filter-branch --prune-empty --subdirectory-filter Caenorhabditis_genomes_project HEAD```

### 6. Navigate to your genomes folder and create a file from your directory list of "scaffolds.fa.gz" files.
```cd ~/genomes```
```ls | rev | cut -d . -f 2- | rev > ~/subproject/CGP_planfile.txt```

### 7. Navigate to the subproject folder and execute the bash script "submit_rmsk_CGP.sh"
```cd ~/subproject/```

```./submit_rmsk_CGP.sh``` 
#Reminder to change directory pointers
#Needs src folder containing rmsk_default_cgp.sbatch, with all directory pointers corrected

### 8. Execute the "submit_hse_helitron_counter.sh".
```./submit_hse_helitron_counter.sh``` 

### 9. Your output should now have a summary table of fractions of HSEs inside or outside of helitrons.
