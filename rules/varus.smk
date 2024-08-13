

rna_files = {key: value['RNA'] for key, value in species_dict.items()}

# make actual lists from the file name string
for species, rnas in rna_files.items():
    if rnas != '':
        rna_files[species] = rnas.split(",")
    else:
        rna_files[species] = []

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