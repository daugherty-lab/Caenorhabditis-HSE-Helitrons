# Caenorhabditis-HSE-Helitrons
This reposiory is a series of HSE-Helitron overlap analyses created by bvtsu for [Garrigues et al., 2019](https://elifesciences.org/articles/51139). This HPC pipeline is controlled by `bash` commands to take a list of UCSC Genome Browser, Wormbase, or Caenorhabditis Genome Project database-specific genome names (e.g. ce11 for Caenorhabditis elegans via UCSC) and uses them to create unique `Torque` jobs and run scripts and pre-built software contained in the `src/` directory.

# Installing/Prerequisites
1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html#macos-installers)
2. Install RepeatMasker via conda (Optional: These genomes already have well-curated rmsk files associated with these databases)
    ```
    conda config --add channels bioconda
    conda install -c bioconda repeatmasker
    ```
3. twoBitToFa (Tool specifically for UCSC Genome Browser genomes)
    ```
    wget -P ~/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/twoBitToFa
    ```

## Usage
Because the datascraping differs by resource, see folder-specific README.md files for additional details
