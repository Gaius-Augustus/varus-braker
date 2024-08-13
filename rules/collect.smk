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

def get_output_name(wildcards):
    archive_path = wildcards.archive_path
    for ext in ['.tar.bz2', '.tar.gz', '.zip', '.tar']:
        if archive_path.endswith(ext):
            return archive_path.rsplit(ext, 1)[0]
    return archive_path + ".unknown"


#convert spiceis and link list ti dictionary
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


#DEBUG
print (species_info)
print ("______________")
   
species_dict = build_fastq_info_dict(species_info)
print (species_dict)

print ("-------------------")
species_dict = {key.replace(" ", "_"): value for key, value in species_dict.items()}

#DEBUG
print(species_dict.keys())
print(species_dict.values())

rule collect_paths:
    input:
        paths = expand('{subdir}/paths.tbl', subdir=subdirs)
    output:
        'paths.txt'
    shell:
        'cat {input.paths} > {output}'

    
rule collect_file_paths:
    input:
        expand('{key}/temp_file0.tbl', key=species_dict.keys()),
        expand('{key}/temp_file_rna.tbl', key=species_dict.keys()),
        
        #"{key}/{acid}.tbl"
    output:
        paths_file="{key}/paths.tbl"
        #"{key}/paths.tbl"
    run:
        import os

        # Assuming species_dict is a dictionary with directory names as keys
        # Let's define or import species_dict here (example below is a placeholder)
        # species_dict = {'dir1': value1, 'dir2': value2, ...}
        print("new")
        directories = list(species_dict.keys())  # Convert dictionary keys to a list
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
        dna_link=lambda wildcards: species_dict[wildcards.f]["DNA"]     
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
                echo "$genome_file" > {wildcards.f}/accession_ID.txt
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
        "{f}/temp_list_rna.tbl"
    params:
        # rna_link now returns a list of URLs or file paths
        rna_link=lambda wildcards: species_dict[wildcards.f]["RNA"]
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
