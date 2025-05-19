[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/isabela42/HyDRA">
    <!--<img src="images/HyDRA.svg" alt="Logo" width=50>-->
  </a>

  <h3 align="center">A Hybrid de novo RNA assembly pipeline</h3>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of contents</h2></summary>
  <ol>
    <li>
      <a href="#overview">Overview</a>
    </li>
    <li><a href="#pipeline-requirements">Pipeline requirements</a></li>
    <ul>
        <li><a href="#tools">Tools</a></li>
        <li><a href="#input-files">Input files</a></li>
      </ul>
    <li><a href="#usage">Usage</a></li>
    <ul>
        <li><a href="#hydra-with-docker-container-image">Docker container image</a></li>
        <li><a href="#hydra-with-terminal-stdin">Terminal stdin</a></li>
        <li><a href="#hydra-with-write-to-pbs-bash-files">Write-to-pbs BASH files</a></li>
      </ul>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## Overview

We developed HyDRA, a hybrid pipeline that integrates bulk short- and long-read RNAseq data for generating custom transcriptomes. This is achieved through (i) read treatment steps to correct sequencing errors by treating low-frequency k-mers and removing contaminants (e.g. adaptors and reads from ribosomal RNAs), (ii) steps to de novo assemble the filtered and corrected reads and further process the resulting assembly, and (iii) optional steps to discover a high-confidence set of lncRNAs supported by multiple machine-learning model predictions.

The scripts are organised as described below, where [S] stands for scripts meant to run with short-read data, [L] for scripts meant to run with long-read data, and [H] for scripts  mean to run with hybrid data.

1. Read treatment scripts

    * Script 01 [S] Quality check raw files (01S1 - FastQC; 01S2 - MultiQC)
    * Script 01 [L] Quality check raw files (01L1 - NanoPlot; 01L2 - NanoComp; 01L3 - FastQC; 01L4 - MultiQC; 01L5 - Fasta_splitter)
    * Script 02 [S] Short-reads correction (02S1 - Rcorrector) & Filter uncorrectable pair end reads (02S2 - FUPER)
    * Script 02 [L] Trim reads of adapters (02L1 - Porechop), average-quality/headcrop (02L2 - Chopper) & poly(A) tails (02L3 - Cutadapt)
    * Script 03 [S] Quality check after correcting (03S1 - FastQC; 03S2 - MultiQC)
    * Script 03 [L] Quality check after trimming (03L1 - NanoPlot; 03L2 - NanoComp; 03L3 - FastQC;  03L4 - MultiQC; 03L5 - Fasta_splitter)
    * Script 04 [S] Trim reads of adapters (04S1 - Trimmomatic) & check ends/average quality (04S2 - BBDuk)
    * Script 04 [L] Long-reads correction (04L1 - RopeBWT2; 04L2 - FMLRC2-convert; 04L3 - FMLRC2)
    * Script 05 [S] Quality check after trimming (05S1 - FastQC; 05S2 - MultiQC)
    * Script 05 [L] Read-lenght check (05L1 - Fasta_splitter)
    * Script 06 [S] Filter rRNA reads (06S1 - BBDuk) & fix read headers (06S1 - in-house script)
    * Script 06 [L] Filter rRNA reads (06L1 - BBDuk)
    * Script 07 [S] Map reads to genome (07S1 - Bowtie2 to flag contamination; 07S1 - SAMtools)
    * Script 07 [L] Map reads to genome (07L1 - Bowtie2 to flag contamination; 07L1 - SAMtools)
    * Script 08 [S] Assess strandness (08S1 - RSeQC) & Assess replicates correlation (08S2 - DeepTools)
    * Script 08 [L] Assess strandness (08L1 - RSeQC)
    * Script 09 [S] Final quality check (09S1 - FastQC; 09S2 - MultiQC)
    * Script 09 [L] Read-lenght check (09L1 - Fasta_splitter)
    * Script 10 [H] Summary of quality control & read error correction steps (10H1)

2. Hybrid de novo assembly scripts

    * Script 11 [H] Hybrid de novo transcriptome assembly (11H1 - rnaSPAdes; 11H2 - Trinity)
    * Script 12 [H] Raw assembly quality check (12H1 - BUSCO; 12H2 - TrinityStats; 12H3 - TransRate; 12H4 - TransRate)
    * Script 13 [H] Transcriptome read representation (13H1 - GMAP & Bowtie2 for rnaSPAdes, 13H2 - Bowtie2 for Trinity assembly)
    * Script 14 [H] Align transcripts to genome and assess splicing (14H1 - GMAP)
    * Script 15 [H] Read support (15H1 - rnaSPAdes based; 15H2 - Trinity based - CDHit, Minimap2, Bowtie2, SAMtools)
    * Script 16 [H] Annotation (16H1 - GMAP, sam2bed, BedTools)
    * Script 17 [H] Summary of hybrid de novo transcriptome assembly steps (17H1)

    Note that steps 11H2, 12H4, 13H2 and 15H2 were designed to create, evaluate and process a short-read-only de novo assembly. These were run during HyDRA development and are based on [Bitar et al. 2023](https://doi.org/10.1093/nar/gkad339).

3. LncRNA discovery scripts

    * Script 18 [H] Predict coding potential (18H1 - ezLncPred)
    * Script 19 [H] Identify long noncoding RNAs (19H1 - FEELnc)
    * Script 20 [H] Define lncRNAs (20H1 - Bash)
    * Script 21 [H] Retrieve metrics, annotation and filter-out protein-coding overlaps (21H1 - BedTools; 21H2 - PBLAT)
    * Script 22 [H] Summary of lncRNA discovery steps (22H1)

<!-- GETTING STARTED -->
## Pipeline requirements 

HyDRA is available in a series of BASH scripts that can be run on i) <a href="#hydra-with-docker-container-image">Docker container image</a>; ii) <a href="#hydra-with-terminal-stdin">terminal stdin</a>; or iii) <a href="#hydra-with-write-to-pbs-bash-files">write-to-pbs BASH files</a>. To run HyDRA with either one of these, you should clone the repo 
  
   ```sh
   git clone https://github.com/isabela42/HyDRA.git
   ```

### Tools

The HyDRA pipeline requires a series of tools to be installed.

Users can choose to build HyDRA's Docker container image. To do so, users will need to install Docker on their machines and run

   ```sh
   docker build -t hydra:1.0.1 .
   ```

  Users can also choose to build the container image directly from VScode. Image building time varies with machine power - in our tests it took from 566s to 1826s. Image size is 22.7GB when build on iMac.

Alternatively, users can choose to locally install all required tools ([docker/tools files](https://github.com/isabela42/HyDRA/tree/main/docker/tools)") using the following command lines:

```sh
bash docker/tools/wget_requirements.sh
bash docker/tools/conda_requirements.sh
conda run --prefix /opt/conda-envs/ezlncpred pip3 install ezlncpred
bash docker/tools/git_requirements.sh
```

Please make sure you meet all tool requirements before running the pipeline.


### Input files

In addition to paired short-read RNAseq and unpaired long-read RNAseq data, users will also need a series of additional input files throughout the pipeline. Please notice that a full description of these is available in [ipda_HyDRA.sh](https://github.com/isabela42/HyDRA/blob/main/ipda_HyDRA.sh)/[ipda_HyDRA_Docker.sh](https://github.com/isabela42/HyDRA/blob/main/ipda_HyDRA_docker.sh).

- adapters_fasta="/path/from/working/dir/to/adapters_trimmomatic.fa"
- ribosomal_rna_ref="/path/from/working/dir/to/ribosomalRNA.fa"
- genome_bowtie2_index="/path/from/working/dir/to/ref_genome_bowtie2"
- reference_gene="/path/from/working/dir/to/ref_gene.bed"
- busco_odb9="/path/from/working/dir/to/BuscoLineages/eukaryota_odb9"
- reference_fasta="/path/from/working/dir/to/ref-transcriptome.fasta"
- reference_gmap_db=reference_GMAP_genome_db_dir
- referencedir=/path/from/working/dir/to/
- transcriptome_bed="/path/from/working/dir/to/ref-transcriptome.bed"
- lncrnadb="/path/from/working/dir/to/lncRNAdb.bed" # we provide a default option you may use [InHouse_lncRNAdb.bed](https://github.com/isabela42/HyDRA/tree/main/docker/data/InHouse_lncRNAdb.bed)
- lncrnadb_fasta="/path/from/working/dir/to/lncRNAdb.fasta" #bedtools getfasta -name -fi genome.fa -bed ${lncrnadb} > lncRNAdb.fasta
- genome_proteins="/path/from/working/dir/to/protein-coding_reference-annotation.bed"
- gencodegenome="/path/from/working/dir/to/gencode.genome.fa" # e.g. <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz>
- gencodebed="/path/from/working/dir/to/gencode.v36.annotation.bed"
- transgtf="/path/from/working/dir/to/ref-transcriptome.gtf"
- ptngtf="/path/from/working/dir/to/gencode.v36.proteincoding.gtf" # Protein-coding transcripts were retrieved using grep "^#\|protein_coding"
- lncrnagtf="/path/from/working/dir/to/gencode.v36.long_noncoding_RNAs.gtf" # e.g. <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.long_noncoding_RNAs.gtf.gz>
- goodlncrnagtf="/path/from/working/dir/to/gencode.v36.confirmed_long_noncoding_RNAs.gtf" # grep -v "TEC" out of above (to be experimentally confirmed transcripts)

## Usage

### HyDRA with Docker container image

HyDRA can be run using a Docker container image. To do so, users may choose to run individual commands on the HyDRA Docker container with [ipda_HyDRA_Docker.sh](https://github.com/isabela42/HyDRA/blob/main/ipda_HyDRA_docker.sh). These commands can be parsed in any terminal running BASH. Please see the help menu with `bash ipda_HyDRA_Docker.sh -h`, replace any necessary input info and parse individual command lines in any terminal running BASH. See info on how to <a href="#tools">build HyDRA's Docker container image here</a>.

Alternatively, users can also choose to i) add <a href="#input-files">all input files</a> to /docker/data folder before ii) building container image, iii) initiate the container with `docker run -it hydra:1.0.1` and iv) run <a href="#hydra-with-terminal-stdin">terminal stdin commands</a>.

### HyDRA with terminal stdin

Users can also choose to execute each one of the commands in [ipda_HyDRA.sh](https://github.com/isabela42/HyDRA/blob/main/ipda_HyDRA.sh), replacing the variables names with respective info. Please see the help menu with `bash ipda_HyDRA.sh -h`, replace any necessary input info and parse individual command lines in any terminal running BASH. Please make sure you have all <a href="#tools">tools</a> instaled and check that you meet their requirements.

### HyDRA with write-to-pbs BASH files

HyDRA can also be run with the series of [write-to-pbs BASH files](https://github.com/isabela42/HyDRA/tree/main/write-to-pbs) to write and submit jobs on a high-performance computer (HPC) managed with a portable batch system (PBS). If you are running HyDRA with these scripts, make sure you have all of the following ready to go before you run any of them:

1. All tools are installed and your local copy of the script is updated with

    * module variables calling and loading the correct tool (either path to local version or on a HPC)
    * additional information session is updated with your paths (scripts 21H1 and 19H1, this will be updated in a next release to take inputs from `-i` main input TSV info file)

2. You will also need to prepare an input TSV file with main input info `-i`. Please see `bash ipda_HyDRA_step*-to-pbs.sh -h` for detailed info of each step requirements.

3. You know how many resources to allocate your jobs. Because HyDRA was designed to run on a PBS, it requires the user to set how much memory `-m`, walltime `-w` and CPUs `-c` each job will use. As a guideline, we provided all resources values used to develop HyDRA. Please note that these may change significantly with different input files, and run tests whenever possible not to overload your machine.

4. Provide email to receive PBS job updates `-e` and stem to name your jobs `-p`

You are now ready to run the BASH scripts of HyDRA using HPC PBS. Each of the files will create a series of PBS files (for each input file in the provided TSV input) and submit those.

```sh
bash ipda_HyDRA_stepX-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"
```

## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<!-- LICENSE -->
## License

Distributed under the MIT License. See [LICENSE][license-url] for more information.

<!-- CONTACT -->
## Contact

Please contact [Isabela Almeida](mb.isabela42@gmail.com) if you have any enquires.

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/isabela42/HyDRA.svg?style=for-the-badge
[contributors-url]: https://github.com/isabela42/HyDRA/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/isabela42/HyDRA.svg?style=for-the-badge
[forks-url]: https://github.com/isabela42/HyDRA/network/members
[stars-shield]: https://img.shields.io/github/stars/isabela42/HyDRA.svg?style=for-the-badge
[stars-url]: https://github.com/isabela42/HyDRA/stargazers
[issues-shield]: https://img.shields.io/github/issues/isabela42/HyDRA.svg?style=for-the-badge
[issues-url]: https://github.com/isabela42/HyDRA/issues
[license-shield]: https://img.shields.io/github/license/isabela42/HyDRA.svg?style=for-the-badge
[license-url]: https://github.com/isabela42/HyDRA/blob/master/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/isabela42/
