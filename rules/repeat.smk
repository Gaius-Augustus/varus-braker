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
