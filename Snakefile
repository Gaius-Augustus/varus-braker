import pandas as pd
import os
#!/usr/bin/env

import subprocess
import time
from time import sleep
import os
import urllib.request
import urllib.error
import gzip
import bz2
import tarfile
import zipfile
import rarfile
import re
import requests
import argparse
import shutil
import random
from multiprocessing import Pool, Semaphore
import configparser
from contextlib import closing
from os import path
import glob

#configfile: 'config/config.json'
#defaults

#parser = argparse.ArgumentParser(description='Run braker for the branch of species')
#parser.add_argument('--input', type=str, help='path to the input file, table should have column with species names, column with links to DNA-data, [optional] column with links to RNA-data',  default="list.txt")
# Initialize the configparser
config = configparser.ConfigParser()

# Load the config.ini file
config.read('config.ini')

#args = parser.parse_args()
#input_file_path = args.input
localrules: collect_paths, find_protein_data, download_dna, download_rna, rename_fasta

species = config.get('MAIN', 'species_info')
species_info = pd.read_csv(species,sep='\t')
clade_list = ['metazoa', 'vertebrata', 'viridiplantae', 'arthropoda', 'eukaryota', 'fungi', 'stramenopiles']

partitition = config.get('SLURM_ARGS', 'partition')
ortho_path = config.get('BRAKER', 'orthodb_path') + '/species/'
excluded = config.get('BRAKER', 'excluded')
#braker_threads = config.get('BRAKER', threads)
#varus_threads = config.get('VARUS', threads)
#braker_partition = 
main_dir = os.getcwd()

varus_path = config.get('VARUS', 'varus_path')                                                                                                                                                                                                                                                                                                                                                                          
hisat2_path = config.get('VARUS', 'hisat2_path')                                                                                                                                                                                                                                                                                                                                                                        
sratoolkit_path = config.get('VARUS', 'sratoolkit_path')                                                                                                                                                                                                                                                                                                                                                                
batchsize = str(config.get('VARUS', 'batchsize'))                                                                                                                                                                                                                                                                                                                                                                       
maxbatches = str(config.get('VARUS', 'maxbatches'))                                                                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                             
hisat2_export_path = f"export PATH={hisat2_path}:$PATH" if hisat2_path else ""                                                                                                                                                                                                                                                                                                                                          
sratoolkit_export_path = f"export PATH={sratoolkit_path}:$PATH" if sratoolkit_path else ""  

include: "rules/collect.smk"
include: "rules/repeat.smk"
include: "rules/varus.smk"
include: "rules/braker.smk"
include: "rules/busco.smk"




        
print("RNA")
print(rna_files)
print("RNA!")
sras = ""



def species_species(wildcards):
    try:
        if isinstance(wildcards, dict):
            file = f"{wildcards['species']}.tmp"
        else:
            file = f"{wildcards.species}.tmp"
        with open(file, 'a') as f:
            pass  # Just to open and close the file, equivalent to open(file, 'a').close()
        return file
    except KeyError:
        return "Error: 'species' key not found in wildcards."
    except AttributeError:
        return "Error: 'species' attribute not found in wildcards."
    except Exception as e:
        return f"An error occurred: {str(e)}"


#Generate a list of all subdirectories at depth level 1
subdirs = [d for d in glob.glob('./*/') if not glob.glob(f'{d}/*/')]  # This excludes deeper nested directories



rule all:
    input:
        #expand('paths.txt', key=species_dict.keys()),
        expand("{key}/braker.gtf", key=species_dict.keys())
    

 
rule test0:
    output:
        "{f}/temp_file"  # Temporary output file, will be renamed later
    params:
        dna_link=lambda wildcards: species_dict[wildcards.f]["DNA"],
    shell:
        """
        # Create the directory if it doesn't exist
        mkdir -p {wildcards.f}
        
        # Determine filename from URL or local path
        filename=$(basename {params.dna_link})
        
        # Check if the dna_link is a URL (HTTP, HTTPS, or FTP)
        if [[ {params.dna_link} =~ ^https?:// ]] || [[ {params.dna_link} =~ ^ftp:// ]]; then
            # It's a URL, download the file
            wget -O {wildcards.f}/${{filename}} {params.dna_link}
        else
            # It's a local file, copy it to the target location with the target filename
            cp {params.dna_link} {wildcards.f}/${{filename}}
        fi
        
        # Since Snakemake expects an output file, touch a temporary file to satisfy Snakemake's requirement.
        # This temp file can be deleted or ignored as needed.
        touch {output.temp}
        """

        # Optional: If you need to work with the actual downloaded/copied file in subsequent rules,
        # consider adding a Python script to your workflow that dynamically generates part of your Snakemake file
        # or rules based on available files and their names.
