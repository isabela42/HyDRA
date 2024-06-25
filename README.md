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
    <li><a href="#pipeline-prerequisites">Pipeline prerequisites</a></li>
    <ul>
        <li><a href="#tools">Tools</a></li>
        <li><a href="#additional-information">Additional information</a></li>
      </ul>
    <li><a href="#usage">Usage</a></li>
    <ul>
        <li><a href="#run-bash-scripts-to-pbs-jobs">Run BASH scripts to PBS jobs</a></li>
        <li><a href="#run-main-bash-command-lines">Run main BASH command lines</a></li>
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
## Pipeline prerequisites

HyDRA is structured in BASH scripts that write and submit portable batch system (PBS) jobs, meaning that our scripts were designed to run in a high-performance computer (HPC) where computational tasks, or simply jobs, are allocated in a PBS system.  You can choose to clone the repo 
  
   ```sh
   git clone https://github.com/isabela42/HyDRA.git
   ```

or <a href="#run-from-comand-line">run individual commands</a>, which can be parsed in any terminal running BASH.

### Tools

The pipeline requires the following tools/versions to be installed. Some of these are dealt with as modules available in the HPC and other as locally installed. Please change to suit your needs.

Comming soon: tools and versions

### Additional information

Users of the pipeline must check script files for this section and update path to additional files required. Below you can find a list of all additional files required to run HyDRA.

Comming soon: additional files

## Usage

### Run BASH scripts to PBS jobs

Make sure you have all of the following ready to go before you run any of the HyDRA scripts:

1. All tools are installed and your local copy of the script is updated with

    * module variables calling the correct tool (either path to local version or on a HPC)
    * additional information session is updated with your reference and supporting files

2. You will also need to prepare an input file with main input info `-i`, starting with war FASTQ files for short-read RNAseq and FASTQ files for long-read RNAseq - each file specifiy which paths and stems are needed in a TSV file.

3. You know how many resources to allocate your jobs. Because HyDRA was designed to run on a PBS, it requires the user to set how much memory `-m`, walltime `-w` and CPUs `-c` each job will use. As a guideline, we provided all resources values used to develop HyDRA. Please note that these may change significantly with different input files, and run tests whenever possible not to overload your HPC.

4. Provide email to receive PBS job updates `-e` and stem to name your jobs `-p`

You are now ready to run the BASH scripts of HyDRA. Each of the files will create a series of PBS files (for each input file in the provided TSV input) and submit those.

```sh
bash script-name.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"
```

### Run main BASH command lines

Alternatively, you can execute each one of the commands, replacing the variables names with respective info.

Comming soon: main command lines.

## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

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
[license-url]: https://github.com/isabela42/HyDRA/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/isabela42/
