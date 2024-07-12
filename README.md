# Structural-Variants
## FullSV_ByPatient.py

### Overview

This code processes structural variant (SV) data from multiple sources and organizes it into individual CSV files for each patient. Each CSV file contains the full details of SVs for a single patient, with each line representing one SV, including two breakpoints. The processed data is saved in the FullSVs_ByPatient directory, which must be created before running the program.

### Prerequisites

Python 3.x 
CSV files from different sources: BCCA, CGI and St Jude SV TARGET datasets

### Usage

python3 FullSVs_ByPatient.py

### Output

The script creates a CSV file for each patient in the FullSVs_ByPatient directory. Each CSV file is named according to the unique patient ID and contains all SVs for that patient.

### License

This project is licensed under the MIT License. See the LICENSE file for details.

## SVs_near_RSSs.py

### Overview

SVs_near_RSSs.py is a script designed to analyze breakpoints of structural variants (SVs) near recombination signal sequences (RSSs) in the genome. This script creates an analysis window spanning 50 bp either side of each breakpoint. The presence of an RSS within the window is then analysed using the DNAGrep algorithm, via RSSsite.

### Requirements

To run this script, you need the following dependencies: Python 3.x 

Standard Libraries: argparse, csv, glob, os, subprocess, time, urllib.request, zipfile, psutil 

External Libraries: selenium Reference genome file in fasta format and corresponding index file are needed

### Usage

Place the input data (bed format) into the working folder and run the following command
python SV_near_RSS.py -g genome.fasta

### Output

The output of the script mainly includes: 
Folder “Bed_files” contain breakpoints of structural variants in bed format 
Folder “Fasta_files” contains sequences of upstream 50bp and downstream 50bp of each breakpoint 
Folder “Unique_RSS_Site_File” contains a csv file that includes breakpoints of structural variants near RSS

### License

This project is licensed under the MIT License. See the LICENSE file for details.

## SVs_near_RSSs_bypatient.py

### Overview

This program analyzes Structural Variants (SVs) and their proximity to Recombination Signal Sequences (RSSs) for patients. It utilizes the output of two other programs: FullSVs_BYPatient.py: Extracts SVs from various datasets and outputs a CSV file for each patient. SV_near_RSS.py: Determines SV breakpoints near RSSs and outputs CSV files with these breakpoints. The program then combines these outputs to generate detailed analyses and summaries

### Requirements

This script requires three folders:

FullSVs_ByPatient: Contains CSV files for each patient, listing SVs with the format:

Column 0: Left Chromosome Column 1: Left Position Column 2: Right Chromosome Column 3: Right Position 

Unique_RSS_Site_files: Contains CSV files listing unique RSS-adjacent sequences. 

Random_unique_sites: Contains CSV files with random distances from RSSs for background analysis.

### Usage

python3 SVs_near_RSSs_bypatient.py

### Output

The program generates the following outputs: Individual_Patient_SVs_ByRSSstatus - A folder containing individual CSV files for each patient with SVs categorized by their proximity to RSSs. 

A summary CSV file for each data source, containing metrics for each patient and background levels of RSS adjacency at specified distances.

### License

This project is licensed under the MIT License. See the LICENSE file for details.
