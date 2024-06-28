localrules: collect_paths, find_protein_data, download_dna, download_rna, rename_fasta
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



partitition = config.get('SLURM_ARGS', 'partition')
species = config.get('MAIN', 'species_info')
species_info = pd.read_csv(species,sep='\t')
clade_list = ['metazoa', 'vertebrata', 'viridiplantae', 'arthropoda', 'eukaryota', 'fungi', 'stramenopiles']
ortho_path = config.get('BRAKER', 'orthodb_path') + '/species/'
excluded = config.get('BRAKER', 'excluded')
#braker_threads = config.get('BRAKER', threads)
#varus_threads = config.get('VARUS', threads)
#braker_partition = 
JOB_LIST = []
main_dir = os.getcwd()

varus_path = config.get('VARUS', 'varus_path')                                                                                                                                                                                                                                                                                                                                                                          
hisat2_path = config.get('VARUS', 'hisat2_path')                                                                                                                                                                                                                                                                                                                                                                        
sratoolkit_path = config.get('VARUS', 'sratoolkit_path')                                                                                                                                                                                                                                                                                                                                                                
batchsize = str(config.get('VARUS', 'batchsize'))                                                                                                                                                                                                                                                                                                                                                                       
maxbatches = str(config.get('VARUS', 'maxbatches'))                                                                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                             
hisat2_export_path = f"export PATH={hisat2_path}:$PATH" if hisat2_path else ""                                                                                                                                                                                                                                                                                                                                          
sratoolkit_export_path = f"export PATH={sratoolkit_path}:$PATH" if sratoolkit_path else ""  


def get_output_name(wildcards):
    archive_path = wildcards.archive_path
    for ext in ['.tar.bz2', '.tar.gz', '.zip', '.tar']:
        if archive_path.endswith(ext):
            return archive_path.rsplit(ext, 1)[0]
    return archive_path + ".unknown"

def build_fastq_info_dict(species_info):
    fastq_info_dict = {}
    for index, row in species_info.iterrows():
        species_id = row['ID']
        dna_link = row['DNA'] if pd.notna(row['DNA']) else ''
        rna_link = row['RNA'] if pd.notna(row['RNA']) else ''
        
        # Include 'ID' inside each species' dictionary
        fastq_info_dict[species_id] = {
            'DNA':  dna_link,
            'RNA':  rna_link,
            'metadata_braker_version': ''
        }
    return fastq_info_dict

def decompress_file(file_path):
    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rb") as gz_file:
            with open(file_path[:-3], "wb") as unzipped_file:
                unzipped_file.write(gz_file.read())
        file_path = file_path[:-3]
    elif file_path.endswith(".bz2") or file_path.endswith(".bzip2"):
        with bz2.open(file_path, "rb") as bz2_file:
            with open(file_path[:-4], "wb") as uncompressed_file:
                uncompressed_file.write(bz2_file.read())
        file_path = file_path[:-4]
    elif file_path.endswith(".tar"):
        with tarfile.open(file_path, "r") as tar:
            tar.extractall()
        file_path = file_path[:-4]
    elif file_path.endswith(".zip"):
        with zipfile.ZipFile(file_path, "r") as zip_ref:
            zip_ref.extractall()
        file_path = file_path[:-4]
    elif file_path.endswith(".rar"):
        with rarfile.RarFile(file_path, "r") as rar_ref:
            rar_ref.extractall()
        file_path = file_path[:-4]
    else:
        print("Uncompressed file of wrong archive format: ", file_path)
        return (file_path)
    return file_path
    
new_data_frame = build_fastq_info_dict(species_info)


new_data_frame = {key.replace(" ", "_"): value for key, value in new_data_frame.items()}

rna_files = {key: value['RNA'] for key, value in new_data_frame.items()}

# make actual lists from the file name string
for species, rnas in rna_files.items():
    if rnas != '':
        rna_files[species] = rnas.split(",")
    else:
        rna_files[species] = []
sras = ""

def rnasrc(wildcards):
    rnas = rna_files.get(wildcards.species, None)
    print(rnas)
    if not rnas:
        # Case 1: Return file path if no RNA data available
        return wildcards.species + "/VARUS.bam"
    elif isinstance(rnas, list):
        
        # Case 3: Return None to indicate no file input is needed, handle IDs directly in the shell command
        print("2: ", rnas)
        sras = ''.join(rnas)
        print("2: ", sras)
        file = 'tmp'
        open(file, 'a').close()
        return file
    else:
        # Case 2: Return file path if it's a single file
        print("3: ", rnas)
        return rnas

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

subdirs = [d for d in glob.glob('./*/') if not glob.glob(f'{d}/*/')]  

rule all:
    input:
        #expand('paths.txt', key=new_data_frame.keys()),
        expand("{key}/braker.gtf", key=new_data_frame.keys())

rule download_dna:
    output:
        "{f}/dna_file.fna"
    params:
        dna_link=lambda wildcards: new_data_frame[wildcards.f]["DNA"]     
    shell:
        """
        # Create the directory if it doesn't exist
        mkdir -p {wildcards.f}
        if [ -z "{params.dna_link}" ]; then
            # dna_link is empty, write "EMPTY" to the file
            taxon_name=$(echo "{wildcards.f}" | sed 's/_/ /g')
            cd {wildcards.f}
            /home/saenkos/datasets/datasets download genome taxon "$taxon_name" --reference
            cd -
            unzip -o {wildcards.f}/ncbi_dataset.zip -d {wildcards.f}
            # Assuming the extracted genome file could have various extensions
            genome_file=$(find {wildcards.f} -type f \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) | head -n 1)
            if [ -n "$genome_file" ]; then
                mv "$genome_file" {wildcards.f}/dna_file.fna
            fi
           
        else
        # Check if the dna_link is a URL (HTTP, HTTPS, or FTP)
            if [[ "{params.dna_link}" =~ ^https?:// ]] || [[ "{params.dna_link}" =~ ^ftp:// ]]; then
                # It's a URL, use wget to download the file
                wget -O {wildcards.f}/dna_file.fna {params.dna_link}
            else
                # It's a local file, copy it to the desired location
                cp {params.dna_link} {wildcards.f}/dna_file.fna
            fi    
        fi
        
        """

rule download_rna:
    output:
        "{f}/{f}_rna.txt",
        "{f}/temp_file_rna.tbl"
    params:
        rna_link=lambda wildcards: new_data_frame[wildcards.f]["RNA"]     
    shell:
        """
        # Create the directory if it doesn't exist
        mkdir -p {wildcards.f}
        
        
        echo {params.rna_link}
        if [ -z "{params.rna_link}" ]; then
            # rna_link is empty, write "EMPTY" to the file
            echo "EMPTY_RNA" > {output[0]}
            echo "EMPTY_RNA2" > {output[1]}
            echo "emmty"
        else
            echo "case 2";
        # Check if the rna_link is a URL (HTTP, HTTPS, or FTP)
            if [[ "{params.rna_link}" =~ ^https?:// ]] || [[ "{params.rna_link}" =~ ^ftp:// ]]; then
                # It's a URL, use wget to download the file
                wget -O {wildcards.f}/rna_file {params.rna_link}
            else
                # It's a local file, copy it to the desired location
                if [[ "{params.rna_link}" =~ ^/ ]]; then
                    cp {params.rna_link} {wildcards.f}/rna_file.fa
                else 

                    echo {params.rna_link} > {wildcards.f}/rna_file.fa
                fi
             fi     
        # Concatenate the downloaded or copied file with the output

        fi
        
        """

rule find_protein_data:
    input:
        species = species_species,
        dna_file = "{species}/dna_file.fna"
    output:
        protein_file = "{species}/proteins.fasta"  # Output path for the protein file
    params:
        excluded = excluded, #Example parameter, replace as needed
        clade_list = ['metazoa', 'vertebrata', 'viridiplantae', 'arthropoda', 'eukaryota', 'fungi', 'stramenopiles'],  # Example list of clades, replace as needed
    shell:
        """
        species_name=$(echo {wildcards.species})
        echo {wildcards.species}
        subdirs=($(ls ./ | grep -x {wildcards.species}))
        echo $subdirs

        # Fallback to fetching data
        taxon_id=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=$species_name" | grep -oP '(?<=<Id>)[^<]+')
        echo $taxon_id
        taxon_id=${{taxon_id:-2759}}
        lineage=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$taxon_id&retmode=xml" | grep -oP '(?<=<Lineage>)[^<]+')
        lineage=${{lineage:-cellular organisms; Eukaryota}}
        echo $lineage

        for clade in $(echo $lineage | tr ';' '\\n' | tac); do
            echo $clade
            echo {params.clade_list}
            clade=$(echo $clade | tr -d '[:space:]')
            if echo "{params.clade_list}" | grep -qwi $clade; then
                echo "Matching clade found: $clade"
                echo "https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/${{clade}}.fa.gz" -O {output.protein_file}.gz
                wget --timeout=30 --tries=3 -q "https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/${{clade}}.fa.gz" -O {output.protein_file}.gz
                gunzip {output.protein_file}.gz
                break
            fi
        done
        """



def check_fasta_needs_masking(path):
    UPPER_PATTERN = r'^[ACGTN\s]+$'
    LOWER_PATTERN = r'^[acgtn\s]+$'
    try:
        with open(path, 'r') as file:
            contents = file.read()
        if re.search(UPPER_PATTERN, contents):
            return "upper"
        elif re.search(LOWER_PATTERN, contents):
            return "lower"
        else:
            return "mixed"
    except UnicodeDecodeError as e:
        print("Error reading the file:", str(e))
        return "error"

def rename_fasta(input_file):
    output_file = os.path.splitext(input_file)[0] + "_renamed" + os.path.splitext(input_file)[1]
    translation_table = {}
    print("blabla")
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        count = 0
        for line in in_file:
            if line.startswith('>'):
                count += 1
                old_header = line.strip()[1:]
                new_header = f'seq{count}'
                translation_table[old_header] = new_header
                out_file.write(f'>{new_header}\n')
            else:
                out_file.write(line)
    with open(input_file + '.translation_table.txt', 'w') as tt_file:
        for old_header, new_header in translation_table.items():
            tt_file.write(f'{old_header}\t{new_header}\n')
    return output_file, translation_table

rule rename_fasta:
    input:
        fasta="{species}/dna_file.fna"
    output:
        renamed_fasta="{species}/dna_file_renamed.fna",
        translation_table="{species}/genome_translation_table.txt"
    run:
        output_file, translation_table = rename_fasta(input.fasta)
        #shell("mv {output_file} {output.renamed_fasta}")
        shell("mv {input.fasta}.translation_table.txt {output.translation_table}")

rule repeat_modeler:
    input:
        fasta="{species}/dna_file_renamed.fna"
    output:
        library="{species}/library.fa"
    log:
        "{species}/logs/repeat_modeler.log"
    threads: 8
    shell:
        """
        echo BuildDatabase -name genome_db {input.fasta}
        echo RepeatModeler -database genome_db -pa {threads} \
            1> {log} 2>&1 && echo RM_*/consensi.fa.classified > {output.library}
        """

rule repeat_masking:
    input:
        fasta="{species}/dna_file_renamed.fna",
        library="{species}/library.fa"
    output:
        softmasked="{species}/dna_file_renamed.softmasked.fna"
    params:
        extra="-xsmall"  # Ensures softmasking
    log:
        "{species}/logs/repeat_masker.log"
    threads: 4
    run:
        status = check_fasta_needs_masking(input.fasta)
        if status == "upper":
            shell(f"""
                echo RepeatMasker -lib {input.library} {params.extra} -pa {threads} \
                    -dir {output.softmasked | dirname} {input.fasta} > {log} 2>&1
            """)
        else:
            shell(f"cp {input.fasta} {output.softmasked}")

 

rule varus:
    input:
        genome="{species}/dna_file.fna"
    output:
        "{species}/VARUS.bam"
    params:
        varus_path=config["VARUS"]["varus_path"],
        hisat2_path=config["VARUS"]["hisat2_path"],
        sratoolkit_path=config["VARUS"]["sratoolkit_path"],
        batchsize=config["VARUS"]["batchsize"],
        maxbatches=config["VARUS"]["maxbatches"],
        partition="snowball",  
        genus=lambda wildcards: wildcards.species.split('_')[0],
        species=lambda wildcards: wildcards.species.split('_')[1],
    threads: 28    
    resources:
        partition="snowball",
        mem_mb=14000        
    shell:
        """
        cd {wildcards.species}
        # Write VARUS parameters to file
        echo "--batchSize {params.batchsize}" > VARUSparameters.txt
        echo "--blockSize 5000" >> VARUSparameters.txt
        echo "--components 1" >> VARUSparameters.txt
        echo "--cost 0.001" >> VARUSparameters.txt
        echo "--deleteLater 0" >> VARUSparameters.txt
        echo "--estimator 2" >> VARUSparameters.txt
        echo "--exportObservationsToFile 1" >> VARUSparameters.txt
        echo "--exportParametersToFile 1" >> VARUSparameters.txt
        echo "--fastqDumpCall fastq-dump" >> VARUSparameters.txt
        echo "--genomeDir ." >> VARUSparameters.txt
        echo "--lambda 10.0" >> VARUSparameters.txt
        echo "--lessInfo 1" >> VARUSparameters.txt
        echo "--loadAllOnce 0" >> VARUSparameters.txt
        echo "--maxBatches {params.maxbatches}" >> VARUSparameters.txt
        echo "--mergeThreshold 10" >> VARUSparameters.txt
        echo "--outFileNamePrefix ./" >> VARUSparameters.txt
        echo "--pathToParameters ./VARUSparameters.txt" >> VARUSparameters.txt
        echo "--pathToRuns ./" >> VARUSparameters.txt
        echo "--pathToVARUS {params.varus_path}/Implementation" >> VARUSparameters.txt
        echo "--profitCondition 0" >> VARUSparameters.txt
        echo "--pseudoCount 1" >> VARUSparameters.txt
        echo "--qualityThreshold 5" >> VARUSparameters.txt
        echo "--randomSeed 1" >> VARUSparameters.txt
        echo "--readParametersFromFile 1" >> VARUSparameters.txt
        echo "--runThreadN 32" >> VARUSparameters.txt
        echo "--verbosityDebug 1" >> VARUSparameters.txt

        # SLURM job script
        export PATH={params.hisat2_path}:$PATH
        export PATH={params.sratoolkit_path}:$PATH

        # Run VARUS
        MINWAIT=10
        MAXWAIT=45
        sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
        {params.varus_path}/runVARUS.pl --aligner=HISAT --readFromTable=0 --createindex=1 --latinGenus={params.genus} --latinSpecies={params.species} --speciesGenome=../{input.genome} --logfile=varus.log 2>varus.err | sbatch
        mv {params.species}/VARUS.bam ../
        cd ..
        """
        
rule braker:
    input:
        dna_path = "{species}/dna_file_renamed.softmasked.fna",
        rna_path = rnasrc,  # Adjust the path if necessary
        proteins_file_path = "{species}/proteins.fasta"
    output:
        gtf = "{species}/braker.gtf"
    params:
        genus=lambda wildcards: wildcards.species.split('_')[0],
        species=lambda wildcards: wildcards.species.split('_')[1],
        augustus_bin_path=config['BRAKER']['augustus_bin_path'],
        augustus_config_path=config['BRAKER']['augustus_config_path'],
        augustus_scripts_path=config['BRAKER']['augustus_scripts_path'],
        diamond_path=config['BRAKER']['diamond_path'],
        prothint_path=config['BRAKER']['prothint_path'],
        genemark_path=config['BRAKER']['genemark_path'],
        module_load=config['SLURM_ARGS']['module_load'],
        braker_cmd=config['BRAKER']['braker_cmd'],
        partition=config['SLURM_ARGS']['partition']
    threads: 16
    resources:
        partition="snowball",
        mem_mb=14000, 
        slurm_extra= "--output=braker.%j.%N.out --error=braker.%j.%N.err  --mail-user=saenkos@uni-greifswald.de"
    shell:
        """
        # Move to appropriate working directory
        cd {wildcards.species}  
        w_dir="{params.genus}_{params.species}_braker"
        RNA_PATH="{input.rna_path}"
        # Check if output already exists

        # Determine RNA-seq input parameters
        rna_subline=""
        if [ -z "{{input.rna_path}}" ]; then
            rna_subline=""
        elif [ "${{RNA_PATH: -3}}" == "bam" ]; then
            rna_subline="--bam=VARUS.bam"
        else 
            echo "${{RNA_PATH: -3}}"
        fi

        # Environment and other parameters setup
        # fix the bug {params.module_load}
        export LC_CTYPE=en_US.UTF-8
        export LC_ALL=en_US.UTF-8

        genemark_export=""
        if [ -n "{params.genemark_path}" ]; then
            genemark_export="export PATH={params.genemark_path}/tools:$PATH"
        fi
        # Run BRAKER command
        {params.braker_cmd} --AUGUSTUS_CONFIG_PATH={params.augustus_config_path} \
                            --AUGUSTUS_BIN_PATH={params.augustus_bin_path} \
                            --AUGUSTUS_SCRIPTS_PATH={params.augustus_scripts_path} \
                            --DIAMOND_PATH={params.diamond_path} \
                            --PROTHINT_PATH={params.prothint_path} --softmasking --useexisting \
                            --GENEMARK_PATH={params.genemark_path} --threads {threads} \
                            --species={params.genus}_{params.species} --workingdir=./$w_dir \
                            --prot_seq=../{input.proteins_file_path} --genome=./$(basename {input.dna_path}) $rna_subline

        # Clean up and move results back to the original directory
        cd -

        """

