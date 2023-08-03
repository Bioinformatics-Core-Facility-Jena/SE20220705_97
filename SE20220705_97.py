# conda create --name SE20220705_97
# conda activate SE20220705_97
# conda install python=3.9
# conda install numpy
# conda install -c bioconda fastqc
# conda install fastp
# conda install -c bioconda kraken2
# conda install -c bioconda trinity
# conda install -c bioconda busco
# conda install -c bioconda augustus
# conda install -c bioconda interproscan
# conda install -c conda-forge -c bioconda mmseqs2
# conda install -c bioconda hisat2
# conda install -c bioconda subread
# conda install -c bioconda bioconductor-deseq2
# conda install -c r r-kernsmooth

import os
from glob import glob
import numpy as np

workingDirectory = ''       #Enter the absolut path of your working directory here. It should contain a folder named "raw_data" containing the raw fastq files


##################################################
##### Quality Control with fastQC and fastp
##################################################

for file in glob(f'{workingDirectory}/raw_data/00_fastq/*gz'):
    os.system(f'fastqc -t 40 {file}')

os.mkdir(f'{workingDirectory}/trimmed_data/')

for file in sorted(glob(f'{workingDirectory}/raw_data/00_fastq/*R1*gz')):
    os.system(f'fastp --in1 {file} --in2 {file.replace("_R1_","_R2_")} --out1 {workingDirectory}/trimmed_data/{os.path.basename(file).replace(".fastq",".trimmed.fastq")} --out2 {workingDirectory}/trimmed_data/{os.path.basename(file).replace("_R1_","_R2_").replace(".fastq",".trimmed.fastq")} --detect_adapter_for_pe --cut_front --cut_tail --n_base_limit=1 --length_required=15 --average_qual=30 --thread=48 --html={workingDirectory}/trimmed_data/{os.path.basename(file).split("_R")[0]}.html')

##################################################
##### Check for read contaminations with Kraken2
##################################################

#make customized kraken2 library by adding two diatoms assemblies
os.mkdir(f'{workingDirectory}/additional_diatoms_genomes/')
os.system(f'wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/405/GCF_000149405.2_ASM14940v2/GCF_000149405.2_ASM14940v2_genomic.fna.gz -O {workingDirectory}/additional_diatoms_genomes/thalassiosira_pseudonana_genome.fna.gz')
os.system(f'gunzip {workingDirectory}/additional_diatoms_genomes/thalassiosira_pseudonana_genome.fna.gz')
os.system(f'curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_018806925.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_018806925.1.zip" -H "Accept: application/zip" -o {workingDirectory}/additional_diatoms_genomes/')
os.system(f'unzip {workingDirectory}/additional_diatoms_genomes/GCA_018806925.1.zip')

#extend the headers of the skeletonema_costatum contigs to make them recognizable for the kraken2 phylo linage
cache = open(f'{workingDirectory}/additional_diatoms_genomes/ncbi_dataset/data/GCA_018806925.1/GCA_018806925.1_FSU_Scostatum_1.0_genomic.fna').readlines()
with open(f'{workingDirectory}/additional_diatoms_genomes/skeletonema_costatum_contigs.fna', 'w') as outFile:
    for line in cache:
        if line[0] == '>':
            outFile.write(f'{line.strip()}|kraken:taxid|2843\n')
        else:
            outFile.write(line)

#build the custom kraken2 database based on the standard databases and the two diatom assemblies
os.mkdir(f'{workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-taxonomy --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library bacteria  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library plasmid  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library viral  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library human  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library protozoa  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'kraken2-build --add-to-library {workingDirectory}/additional_diatoms_genomes/thalassiosira_pseudonana_genome.fna --db {workingDirectory}/kraken2_db/')
os.system(f'kraken2-build --add-to-library {workingDirectory}/additional_diatoms_genomes/thalassiosira_oceanica_contigs.fna --db {workingDirectory}/kraken2_db/')
os.system(f'kraken2-build --add-to-library {workingDirectory}/additional_diatoms_genomes/skeletonema_costatum_contigs.fna --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --build --threads 48  --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --clean --threads 48  --db {workingDirectory}/kraken2_db/')

#run kraken2
for file in sorted(glob(f'{workingDirectory}/trimmed_data/*gz')):
    print(f'kraken2 --db {workingDirectory}/kraken2_db/ --threads 48 --use-names --output {workingDirectory}/kraken2_db/{os.path.basename(file)}_kraken2_output --report {workingDirectory}/kraken2_db/{os.path.basename(file)}_kraken2_report {file}')

### Generate a short basic overview

krakenOverview = {}
for file in glob(f'{workingDirectory}//kraken2_db/*report'):
    currentName = os.path.basename(file).split('_001')[0]
    krakenOverview[currentName] = {'unclassified': (0,0), 'Eukaryota': (0,0), 'Bacteria': (0,0), 'Viruses': (0,0)}
    for line in open(file).readlines():
        line = line.strip().split('\t')
        if line[-1].strip() in krakenOverview[currentName]:
            krakenOverview[currentName][line[-1].strip()] = (float(line[0]), int(line[1]))

with open(f'{workingDirectory}/kraken2_db/basic_overview.csv', 'w') as outfile:
    outfile.write('File Name\tUnclassified %\tEukaryota %\tBacteria %\tViruses %\n')
    for name in sorted(krakenOverview):
        outfile.write(f'{name}\t{krakenOverview[name]["unclassified"][0]}\t{krakenOverview[name]["Eukaryota"][0]}\t{krakenOverview[name]["Bacteria"][0]}\t{krakenOverview[name]["Viruses"][0]}\n')


##### Generate plots from kraken2 results

for file in glob(f'{workingDirectory}/kraken2_db/*_R1_*report'):
    os.system(f'ktImportTaxonomy -t 5 -m 3 -o {os.path.basename(file).split("_")[0]}_krona.html {file} {file.replace("_R1_", "_R2_")}')


##################################################
##### Trinity assemblies
##################################################
os.makedirs(f'{workingDirectory}/trinity/trinity_s_costatum/')

sCostatumSamples = ['1ASc','2ASc','3A1Sc','4ASc','9ASc','10BSc','11ASc','12BSc']

samplesFiles = list(glob(f'{workingDirectory}/trimmed_data/*gz'))
os.system(f"Trinity --seqType fq --max_memory 20G --output {workingDirectory}/trinity/trinity_s_costatum/ --CPU 48 --min_contig_length 200 --full_cleanup --left {workingDirectory}//trimmed_data/1ASc_R1_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/2ASc_R1_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/3A1Sc_R1_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/4ASc_R1_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/9ASc_R1_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/10BSc_R1_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/11ASc_R1_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/12BSc_R1_001.trimmed.fastq.gz --right {workingDirectory}//trimmed_data/1ASc_R2_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/2ASc_R2_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/3A1Sc_R2_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/4ASc_R2_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/9ASc_R2_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/10BSc_R2_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/11ASc_R2_001.trimmed.fastq.gz,{workingDirectory}//trimmed_data/12BSc_R2_001.trimmed.fastq.gz")

##### get some basic stats

contigs = []
for line in open(f'{workingDirectory}/trinity/trinity_s_costatum.Trinity.fasta'):
    if line[0] == '>':
        contigs.append(int(line.split('len=')[1].split()[0]))
print(f'S_costatum; #Transcripts {len(contigs)}; Mean Length {np.mean(contigs)}; Median Length {np.median(contigs)}')

##################################################
##### reclassify assembled transcripts with Kraken2
##################################################

os.makedirs(f'{workingDirectory}/kraken2/s_costatum/')

for file in sorted(glob(f'{workingDirectory}/trinity/*fasta')):
    currentName = os.path.basename(file).split('.')[0].split('y_')[1]
    os.system(f'kraken2 --db {workingDirectory}//kraken2_db/ --threads 48 --use-names --output {workingDirectory}/kraken2/{currentName}/{currentName}_kraken2_output --report {workingDirectory}/kraken2/{currentName}/{currentName}_kraken2_report {file}')

### Generate a short basic overview

krakenOverview = {}
for file in glob(f'{workingDirectory}//kraken2/*/*report'):
    currentName = os.path.basename(file).split('_kraken2')[0]
    krakenOverview[currentName] = {'unclassified': (0,0), 'Eukaryota': (0,0), 'Bacteria': (0,0), 'Viruses': (0,0)}
    for line in open(file).readlines():
        line = line.strip().split('\t')
        if line[-1].strip() in krakenOverview[currentName]:
            krakenOverview[currentName][line[-1].strip()] = (float(line[0]), int(line[1]))

with open(f'{workingDirectory}/kraken2/basic_overview.csv', 'w') as outfile:
    outfile.write('File Name\tUnclassified %\tEukaryota %\tBacteria %\tViruses %\n')
    for name in sorted(krakenOverview):
        outfile.write(f'{name}\t{krakenOverview[name]["unclassified"][0]}\t{krakenOverview[name]["Eukaryota"][0]}\t{krakenOverview[name]["Bacteria"][0]}\t{krakenOverview[name]["Viruses"][0]}\n')

##### Generate plots from kraken2 results

for file in glob(f'{workingDirectory}/kraken2/*/*report'):
    currentName = os.path.basename(file).split('_kraken2')[0]
    os.system(f'ktImportTaxonomy -t 5 -m 3 -o {workingDirectory}/kraken2/{currentName}_krona.html {file}')

##################################################
##### sort assembled transcripts by classification
##################################################

taxidDict = {'Bacteria': set(), 'Eukaryota': set(), 'Viruses': set(), 'Human': ('9606',)}

currentDomain = ''
for line in open(f'{workingDirectory}/kraken2/s_costatum/s_costatum_kraken2_report'):
    line = line.strip().split()
    currentDomain = line[-1] if line[-1] in taxidDict else currentDomain
    if currentDomain:
        taxidDict[currentDomain].add(line[4])

taxidDict['Stramenopiles'] = set()
collect_taxids = False
for line in open(f'{workingDirectory}/kraken2/s_costatum/s_costatum_kraken2_report'):
    line = line.strip().split()
    if collect_taxids:
        if line[3] == 'D2':
            print(line)
            break
        taxidDict['Stramenopiles'].add(line[4])
    if line[-1] == 'Stramenopiles':
        taxidDict['Stramenopiles'].add(line[4])
        collect_taxids = True
        continue

transcriptIDs = {'Bacteria': set(), 'Eukaryota': set(), 'Viruses': set(), 'Human': set(), 'Stramenopiles': set()}
for line in open(f'{workingDirectory}/kraken2/s_costatum/s_costatum_kraken2_output'):
    currentTaxID = line.split('taxid ')[1].split(')')[0]
    currentTranscriptID = line.split()[1]
    if currentTaxID in taxidDict['Bacteria']:
        transcriptIDs['Bacteria'].add(currentTranscriptID)
    elif currentTaxID in taxidDict['Viruses']:
        transcriptIDs['Viruses'].add(currentTranscriptID)
    elif currentTaxID in taxidDict['Human']:
        transcriptIDs['Human'].add(currentTranscriptID)
    elif currentTaxID in taxidDict['Stramenopiles']:
        transcriptIDs['Stramenopiles'].add(currentTranscriptID)
        transcriptIDs['Eukaryota'].add(currentTranscriptID)
    else:
        transcriptIDs['Eukaryota'].add(currentTranscriptID)

transcriptIDs['Unclassified'] = set()
for line in open(f'{workingDirectory}/kraken2/s_costatum/s_costatum_kraken2_output'):
    if line[0] == 'U':
        transcriptIDs['Unclassified'].add(line.split()[1])

sequence_switch = False
for key in taxidDict:
    with open(f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_{key}.fasta', 'w') as outfile:
        for line in open(f'{workingDirectory}/trinity/trinity_s_costatum/trinity_s_costatum.Trinity.fasta'):
            if line[0] == '>':
                if line.split()[0][1:] in transcriptIDs[key]:
                    outfile.write(line)
                    sequence_switch = True
                else:
                    sequence_switch = False
            elif sequence_switch:
                outfile.write(line)



##################################################
##### Evaluate Transcriptome Assemblies
##################################################

##### Busco evaluation

os.makedirs(f'{workingDirectory}/assembly_evaluation/s_costatum/')
os.system(f'busco -m transcriptome --in {workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles.fasta --out_path {workingDirectory}/assembly_evaluation/s_costatum/ --out s_costatum_stramenopiles -c 40 -l stramenopiles_odb10')
os.system(f'generate_plot.py -wd {workingDirectory}/assembly_evaluation/p_parvum/')

##################################################
##### Gene Prediction using Augustus
##################################################

os.system(f'augustus --strand=both --gff3=on --outfile={workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_prediction.gff3 --species=Skeletonema_costatum {workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles.fasta')

### parse augustus output into actual gff and faa files

def extract_gff_and_fasta_from_augusutus_output(augustus_output_file, gff_output_file, fasta_output_file):
    sequence_switch = False
    with open(gff_output_file, 'w') as gff_outfile, open(fasta_output_file, 'w') as fasta_outfile:
        for line in open(augustus_output_file):
            if line[0] != '#':
                gff_outfile.write(line)
            if 'start gene ' in line:
                current_gene_name = line.split()[-1].strip()
                current_sequence = ''
            if 'Evidence for and against this transcript' in line:
                sequence_switch = False
                fasta_outfile.write(f'>{current_gene_name}\n{current_sequence}\n')
            if sequence_switch:
                current_sequence += line.split()[1].split(']')[0].strip()
            if 'protein sequence =' in line:
                current_sequence += line.split('[')[1].split(']')[0].strip()
                sequence_switch = True

extract_gff_and_fasta_from_augusutus_output(f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_prediction.gff3',
                                            f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_annotation.gff3',
                                            f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_prediction.faa')

##################################################
##### Gene Annotation using InterProScan
##################################################

os.system(f'interproscan.sh -hm -goterms –pathways -appl Pfam,TIGRFAM,PANTHER,ProSiteProfiles,FunFam -f TSV, GFF3 -cpu 60 -i {workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_prediction.faa –d {workingDirectory}/trinity/trinity_s_costatum/')

##### get a simple count on the classified predicted genes

cache = {}
for line in open(f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_classification.tsv'):
    cache[line.split('\t')[0]] = line.split('\t')[1]
print(len(cache))

# 32478

##################################################
##### cluster genes based on their predicted protein sequences using mmseqs2 easy-cluster
##################################################

os.system(f'mmseqs easy-cluster {workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_prediction.faa {workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_prediction_clustered {workingDirectory}/trinity/trinity_s_costatum/tmp/ --min-seq-id 0.8 --cov-mode 1 -c 0.3 --threads 60')

s_costatum_clustered_genes = {}
for line in open(f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_prediction_clustered_cluster.tsv'):
    if line.split()[0] not in s_costatum_clustered_genes:
        s_costatum_clustered_genes[line.split()[0]] = set()
    s_costatum_clustered_genes[line.split()[0]].add(line.split()[1])

s_costatum_geneID_to_transcriptID = {}
for line in open(f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_annotation.gff3'):
    line = line.strip().split('\t')
    if line[2] == 'gene':
        s_costatum_geneID_to_transcriptID[line[-1].split('ID=')[1]] = line[0]

s_costatum_reduced_transcript_IDs = set()
for geneID in s_costatum_clustered_genes:
    s_costatum_reduced_transcript_IDs.add(s_costatum_geneID_to_transcriptID[geneID])


os.makedirs(f'{workingDirectory}/deg_analysis/s_costatum/')

# create a dict from fasta file
s_costatum_fasta_dict = {}
for line in open(f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles.fasta'):
    if line[0] == '>':
        currentID = line.split()[0][1:]
        s_costatum_fasta_dict[currentID] = ''
    else:
        s_costatum_fasta_dict[currentID] += line.strip()

with open(f'{workingDirectory}/deg_analysis/s_costatum/s_costatum_Stramenopiles_reduced.fna', 'w') as outfile:
    for transcriptID in s_costatum_reduced_transcript_IDs:
        outfile.write(f'>{transcriptID}\n{s_costatum_fasta_dict[transcriptID]}\n')

##################################################
##### Mapping Samples to the reduced transcriptomes using HISAT2
##################################################

sample_conditions = {'1': ['1ASc','2ASc','3A1Sc','4ASc'],
                     '2': ['9ASc','10BSc','11ASc','12BSc']}

##### s_costatum

### create a hisat2 index
os.system(f'hisat2-build -p 60 {workingDirectory}/deg_analysis/s_costatum/s_costatum_Stramenopiles_reduced.fna {workingDirectory}/deg_analysis/s_costatum/s_costatum_Stramenopiles_reduced.idx')

### perform mapping
os.makedirs(f'{workingDirectory}/deg_analysis/s_costatum/mapped_reads/')

### map files that belong to condition 1 and 2
for condition in ['1', '2']:
    for file in glob(f'{workingDirectory}/trimmed_data/*R1*fastq.gz'):
        currentName = os.path.basename(file).split('_')[0]
        if currentName in sample_conditions[condition]:
            os.system(f'hisat2 -p 60 -x {workingDirectory}/deg_analysis/s_costatum/s_costatum_Stramenopiles_reduced.idx -1 {file} -2 {file.replace("_R1_","_R2_")} -S {workingDirectory}/deg_analysis/s_costatum/mapped_reads/{currentName}.sam')

##################################################
##### Read counting using featureCounts
##################################################

##### s_costatum
os.makedirs(f'{workingDirectory}/deg_analysis/s_costatum/read_counts/')
os.system(f'featureCounts -T 60 -p -O -M -t gene -g ID -a {workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_annotation.gff3 -o {workingDirectory}/deg_analysis/s_costatum/read_counts/s_costatum_Stramenopiles_read_counts.txt {workingDirectory}/deg_analysis/s_costatum/mapped_reads/*sam')

### write sample reads into individual files
count_file = open(f'{workingDirectory}/deg_analysis/s_costatum/read_counts/s_costatum_Stramenopiles_read_counts.txt').readlines()[1:]
samples = [x.split('/')[-1].split('.')[0] for x in count_file[0].strip().split()[6:]]
file_handles = {sample: open(f'{workingDirectory}/deg_analysis/s_costatum/read_counts/{sample}.txt', 'w') for sample in samples}

for line in count_file[1:]:
    line = line.strip().split()
    for i, sample in enumerate(samples):
        file_handles[sample].write(f'{line[0]}\t{line[6+i]}\n')

for sample in file_handles:
    file_handles[sample].close()

##################################################
##### DEG analysis using DESeq2
##################################################

##### s_costatum

os.makedirs(f'{workingDirectory}/deg_analysis/s_costatum/DESeq2/')
os.system(f'Rscript {workingDirectory}/deseq2_comparison.r {workingDirectory}/deg_analysis/s_costatum/read_counts/s_costatum_Stramenopiles_read_counts.txt {workingDirectory}/deg_analysis/s_costatum/DESeq2/')

### add gene information to the DESeq2 output
gene_information = {}
for line in open(f'{workingDirectory}/trinity/trinity_s_costatum/s_costatum_Stramenopiles_augustus_gene_classification.tsv'):
    line = line.strip().split('\t')
    if line[0] not in gene_information:
        gene_information[line[0]] = set()
    gene_information[line[0]].add(line[5].lower() if line[5] != '-' else '')
    gene_information[line[0]].add(line[12].lower() if line[12] != '-' else '')
    if len(line) == 14:
        gene_information[line[0]].add(line[13].lower() if line[13] != '-' else '')

infile = open(f'{workingDirectory}/deg_analysis/s_costatum/DESeq2/deseq2_result.csv').readlines()[1:]
with open(f'{workingDirectory}/deg_analysis/s_costatum/DESeq2/deseq2_result.csv', 'w') as outfile:
    outfile.write('"","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","predicted function/predicted domain"\n')
    for line in infile:
        outfile.write(f'{line.strip()},{";".join(gene_information.get(line.split(",")[0][1:-1], ["unknown function"]))}\n')