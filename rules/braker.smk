# This rule determines the BRAKER run mode from the input
# file species_csv, automatically. It runs either BRAKER3,
# or BRAKER2 (if present including legacy proteins from 
# a previously annotated assembly of that species).
# It is not running BRAKER1, we hope that this will not be
# necessary / i.e. will resort to BRAKER2 if BRAKER1 fails.
# Failure is not caught/resumed.
rule run_braker:
    input:
        input = "data/checkpoints_annotate/{spid}_repeats.done"
    output:
        done = "data/checkpoints_annotate/{spid}_braker.done"
    params:
        spid="{spid}",
        threads = config['SLURM_ARGS']['cpus_per_task'],
        csv = config['INPUT']['species_csv']
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        cmd_file=data/checkpoints_annotate/{params.spid}_braker.cmd
        log_file=data/checkpoints_annotate/{params.spid}_braker.log
        wdir=data/species/{params.spid}/braker
        python3 scripts/build_annotation_cmd.py -c {params.csv} -s {params.spid} -t {params.threads} -o $cmd_file -l $log_file -w $wdir
        bash $cmd_file # if you comment this line, braker will not be executed
        touch {output.done}
        """
