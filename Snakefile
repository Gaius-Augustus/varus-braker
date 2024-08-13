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


print (species_info)
print ("______________")

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
    print(fastq_info_dict)
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
print (new_data_frame)

print ("-------------------")
new_data_frame = {key.replace(" ", "_"): value for key, value in new_data_frame.items()}

print(new_data_frame.keys())
print(new_data_frame.values())

rna_files = {key: value['RNA'] for key, value in new_data_frame.items()}

# make actual lists from the file name string
for species, rnas in rna_files.items():
    if rnas != '':
        rna_files[species] = rnas.split(",")
    else:
        rna_files[species] = []
        
print("RNA")
print(rna_files)
print("RNA!")
sras = ""

def rnasrc(wildcards):
    rnas = rna_files.get(wildcards.species, None)
    print(rnas)
    if not rnas:
        # Case 1: Return file path if no RNA data available
        return wildcards.species + "/VARUS.bam"
    elif isinstance(rnas, list):
        if any(':' in item or '/' in item or '\\' in item for item in rnas):
            return 1
        else:
            return 0
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


# First, we generate a list of all subdirectories at depth level 1
subdirs = [d for d in glob.glob('./*/') if not glob.glob(f'{d}/*/')]  # This excludes deeper nested directories



rule all:
    input:
        #expand('paths.txt', key=new_data_frame.keys()),
        expand("{key}/braker.gtf", key=new_data_frame.keys())
        

#rule download_file:
#    output:
#        "downloaded_file" # You can rename this to match the file you're downloading or keep it generic
#    params:
#        url = "your_download_link_here"
#    shell:
#        """
#        wget -O {output} {params.url}
#        """


#rule download:

#rule uncompress:
#    output:    

rule test0:
    output:
        "{f}/temp_file"  # Temporary output file, will be renamed later
    params:
        dna_link=lambda wildcards: new_data_frame[wildcards.f]["DNA"],
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
        

#rule varus:
#    output: ('{f}/')

rule collect_paths:
    input:
        paths = expand('{subdir}/paths.tbl', subdir=subdirs)
    output:
        'paths.txt'
    shell:
        'cat {input.paths} > {output}'

    
rule collect_file_paths:
    input:
        expand('{key}/temp_file0.tbl', key=new_data_frame.keys()),
        expand('{key}/temp_file_rna.tbl', key=new_data_frame.keys()),
        
        #"{key}/{acid}.tbl"
    output:
        paths_file="{key}/paths.tbl"
        #"{key}/paths.tbl"
    run:
        import os

        # Assuming new_data_frame is a dictionary with directory names as keys
        # Let's define or import new_data_frame here (example below is a placeholder)
        # new_data_frame = {'dir1': value1, 'dir2': value2, ...}
        print("new")
        directories = list(new_data_frame.keys())  # Convert dictionary keys to a list
        print("DIRS ", directories )
        # Open the output file as specified in the rule's output
        with open(output.paths_file, 'w') as paths_file:
            # Iterate over each directory in the list
            for directory in directories:
                # Check if the directory exists to avoid errors
                if os.path.exists(directory):
                    # Walk through each specified directory
                    for dirpath, dirnames, filenames in os.walk(directory):
                        # For each file, write its path to the paths_file
                        for filename in filenames:
                            # Construct the full path
                            file_path = os.path.join(dirpath, filename)
                            # Write the path to the file
                            paths_file.write(file_path + '\n')

rule download_dna:
    output:
        "{f}/dna_file.fna"
    params:
        dna_link=lambda wildcards: new_data_frame[wildcards.f]["DNA"]     
    shell:
        """
        # Create the directory if it doesn't exist
        mkdir -p {wildcards.f}
        echo {params.dna_link}
        echo "YOU ARE HERE"
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
            echo "case";
            # Check if the dna_link is a URL (HTTP, HTTPS, or FTP)
            if [[ "{params.dna_link}" =~ ^https?:// ]] || [[ "{params.dna_link}" =~ ^ftp:// ]]; then
            # It's a URL, use curl to download the file
                curl -o {wildcards.f}/dna_file.tmp {params.dna_link}
                # Check file type and unarchive if necessary
                if file {wildcards.f}/dna_file.tmp | grep -q 'gzip compressed'; then
                    if file {wildcards.f}/dna_file.tmp | grep -q '.tar'; then
                        cp {wildcards.f}/dna_file.tmp {wildcards.f}/dna_file.tar.gz
                        tar -xzvf {wildcards.f}/dna_file.tar.gz -C {wildcards.f}
                        #mv {wildcards.f}/*.fna {wildcards.f}/dna_file.fna
                    else 
                        cp {wildcards.f}/dna_file.tmp {wildcards.f}/dna_file.gz
                        gzip -d -c {wildcards.f}/dna_file.gz > {wildcards.f}/dna_file.fna
                        #mv {wildcards.f}/*.fna {wildcards.f}/dna_file.fna
                    fi     
                elif file {wildcards.f}/dna_file.tmp | grep -q 'Zip archive'; then
                    cp {wildcards.f}/dna_file.tmp {wildcards.f}/dna_file.zip
                    unzip {wildcards.f}/dna_file.zip -d {wildcards.f}
                    #mv {wildcards.f}/*.fna {wildcards.f}/dna_file.fna
                else
                    mv {wildcards.f}/dna_file.tmp {wildcards.f}/dna_file.fna
                fi
            else
                # It's a local file, copy it to the desired location
                cp {params.dna_link} {wildcards.f}/dna_file.fna
            fi
        fi
        """

rule download_rna:
    output:
        #"{f}/list_rna.txt",
        "{f}/temp_list_rna.tbl"
    params:
        # rna_link now returns a list of URLs or file paths
        rna_link=lambda wildcards: new_data_frame[wildcards.f]["RNA"]
    shell:
        """
        # Create the directory if it doesn't exist
        mkdir -p {wildcards.f}
        
        # Initialize temp_list_rna.tbl to store downloaded file paths
        : > {output[0]}
        
        # Process each link in the list
        for link in $(echo {params.rna_link} | tr ',' '\n'); do
            echo "Processing $link"
            if [[ "$link" =~ ^https?:// ]] || [[ "$link" =~ ^ftp:// ]]; then
                # It's a URL, use wget to download the file
                echo "curl"
                if [ -e "{wildcards.f}/$(basename $link .gz)" ]; then
                    echo 'File already exists' >&2
                else
                    curl -o {wildcards.f}/$(basename $link) $link
                fi
                
                echo "curl done"
                if file {wildcards.f}/$(basename $link) | grep -q 'gzip compressed'; then
                    if file {wildcards.f}/$(basename $link) | grep -q '.tar'; then
                        #cp {wildcards.f}/$(basename $link) {wildcards.f}/rna_file.tar.gz
                        tar -xzvf {wildcards.f}/$(basename $link) -C {wildcards.f}/$(basename $link .tar.gz)
                        #mv {wildcards.f}/*.fna {wildcards.f}/dna_file.fna
                        echo $(basename $link .tar.gz) >> {output[0]}
                        #echo "{wildcards.f}/$(basename $link .tar.gz)" >> {output[0]}
                    else 
                        #cp {wildcards.f}/$(basename $link) {wildcards.f}/rna_file.gz
                        #gzip -d -c {wildcards.f}/$(basename $link) > {wildcards.f}/$(basename $link .gz)
                        echo $(basename $link .gz) >> {output[0]}
                        #echo "{wildcards.f}/$(basename $link .gz)" >> {output[0]}
                        #mv {wildcards.f}/*.fna {wildcards.f}/dna_file.fna
                    fi     
                elif file {wildcards.f}/dna_file.tmp | grep -q 'Zip archive'; then
                    cp {wildcards.f}/dna_file.tmp {wildcards.f}/dna_file.zip
                    unzip {wildcards.f}/dna_file.zip -d {wildcards.f}
                    #mv {wildcards.f}/*.fna {wildcards.f}/dna_file.fna
                #wget -O {wildcards.f}/$(basename $link) $link
                fi
                
            else
                # It's a local file, copy it to the desired location
                if [[ "$link" =~ ^/ ]]; then
                    cp $link {wildcards.f}/$(basename $link)
                    echo "{wildcards.f}/$(basename $link)" >> {output[0]}
                else
                    echo "SRA ID: $link"
                    echo $link  >> {output[0]}
                fi
            fi
        done

        
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
        genome="{species}/dna_file_renamed.softmasked.fna"
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
        echo "PREVARUS"
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
        rna_list = "{species}/temp_list_rna.tbl", 
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
        partition=config['SLURM_ARGS']['partition'],
        rnasrc=rnasrc
        
    threads: 16
    resources:
        partition="snowball",
        mem_mb=14000 
        #slurm_extra= "--output={species}/braker.%j.%N.out --error=braker.%j.%N.err  --mail-user=saenkos@uni-greifswald.de"
    shell:
        """
        # Move to appropriate working directory
        cd {wildcards.species}  
        w_dir="{params.genus}_{params.species}_braker"
        RNA_PATH="{input.rna_list}"
        # Check if output already exists
        echo "PARAMS {params.rnasrc}" 
        # Determine RNA-seq input parameters
        rna_subline=""
        echo {input.rna_list}
        echo "599!"
        echo $w_dir
        paths=$(cat ../{input.rna_list} | tr '\\n' ',' | sed 's/,$//')
        echo $paths
        if [ -z "{{input.rna_list}}" ]; then
            rna_subline=""
            
        elif [ "${{RNA_PATH: -3}}" == "bam" ]; then
            rna_subline="--bam=VARUS.bam --rnaseq_sets_dirs="
        else 
            if [ "{params.rnasrc}" == "1" ] ; then
                echo "!!!!!!!!!!!!!!!!!!!1111111111111111111"
                cleaned_paths=$(echo "$paths" | sed -E 's/(\.[^.,]*)//g')
                echo $cleaned_paths
                rna_subline="--rnaseq_sets_ids=$cleaned_paths --rnaseq_sets_dirs=./"
            else
                rna_subline="--rnaseq_sets_ids=$paths"
            fi
        fi
        echo "605!"
        echo $rna_subline
        echo "607!"
        # Environment and other parameters setup
        echo {params.module_load}
        export LC_CTYPE=en_US.UTF-8
        export LC_ALL=en_US.UTF-8
        echo "619!"
        genemark_export=""
        if [ -n "{params.genemark_path}" ]; then
            genemark_export="export PATH={params.genemark_path}/tools:$PATH"
        fi
        id=$RANDOM

        # Prepare temporary directory for processing
        #mkdir /tmp/saenkos-$id

        echo "632!"
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
