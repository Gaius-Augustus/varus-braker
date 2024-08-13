# This rule determines the BRAKER run mode from the input
# file species_csv, automatically. It runs either BRAKER3,
# or BRAKER2 (if present including legacy proteins from 
# a previously annotated assembly of that species).
# It is not running BRAKER1, we hope that this will not be
# necessary / i.e. will resort to BRAKER2 if BRAKER1 fails.
# Failure is not caught/resumed.
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
