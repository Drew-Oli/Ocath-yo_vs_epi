#!/usr/bin/env python3

import multiprocessing

# wild cards
SAMPLES, = glob_wildcards('data/reads/{sample}_r1.fq')

# singularity containers
trinity = 'docker://trinityrnaseq/trinityrnaseq:2.12.0'
cutadapt = 'docker://quay.io/biocontainers/cutadapt:2.5--py37h516909a_0'
salmontools = 'shub://TomHarrop/align-utils:salmontools_23eac84'
salmon = 'docker://combinelab/salmon:1.2.1'
transdecoder = 'docker://quay.io/biocontainers/transdecoder:5.5.0--pl526_2'
emapper = 'docker://quay.io/biocontainers/eggnog-mapper:2.1.2--pyhdfd78af_0'
busco = 'docker://ezlabgva/busco:v5.1.2_cv1'
bioconductor = 'shub://TomHarrop/singularity-containers:bioconductor_3.9'
star = 'docker://quay.io/biocontainers/star:2.7.9a--h9ee0642_0'
samtools = 'docker://quay.io/biocontainers/samtools:1.2--0'

# rules
rule target:
    input:    
        'output/04_deseq2/ocath-genome_deseq2_result_wannot.csv',
#        'output/011_deseq2/ocath-trinity-gg_deseq2_result_wannot.csv'

rule deseq2_trinitygg:
    input:
        quant_files = expand('output/07_salmon/quant/{sample}/quant.sf', sample = SAMPLES),
        salmon_dir = 'output/07_salmon',
        samples_file = 'output/samples.txt',
        gene_ids = 'output/06_trinity_gg/Trinity.fasta.gene_trans_map',
        emapper_annotation_file = 'output/09_emapper/Trinity.emapper.annotations'
    output:
        dds_file = 'output/011_deseq2/dds.Rds',
        deseq2_result_file = 'output/011_deseq2/ocath-trinity-gg_result.csv',
        deseq2_count_file = 'output/011_deseq2/ocath-trinity-gg_count.csv',
        deseq2_result_with_annotation_file = 'output/011_deseq2/ocath-trinity-gg_deseq2_result_wannot.csv'
    log:
        'output/logs/011_deseq2/deseq2_trinity-gg.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bioconductor
    script:
        'scripts/deseq2_trinity-gg.R'
            
rule busco:
    input:
        'output/08_transdecoder/Trinity.fasta.transdecoder.pep'
    output:
        'output/010_busco/busco/short_summary.specific.arthropoda_odb10.busco.txt'
    params:
        dir = 'output/010_busco',
        proteins = 'Trinity.fasta.transdecoder.pep',
        output = 'busco',
        lineage = 'arthropoda_odb10',
        mode = 'protein',
        cd = '../../'
    log:
        'output/logs/010_busco/busco.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        busco
    shell:
        'cp {input} {params.dir} ; '        
        'cd {params.dir} || exit 1 ; '
        'busco '
        '-i {params.proteins} '
        '-o {params.output} '
        '-l {params.lineage} '
        '-m {params.mode} '
        '-f ; '
        'cd {params.cd} ; '
        '&> {log}'

rule emapper_trinitygg:
    input:
        'output/08_transdecoder/Trinity.fasta.transdecoder.pep'
    output:
        'output/09_emapper/Trinity.emapper.annotations'
    params:
        database = 'data/ref/',
        method = 'diamond',
        filename = 'output/09_emapper/Trinity'
    log:
        'output/logs/09_emapper/emapper_trans.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        emapper
    shell:
#        'download_eggnog_data.py '
#        '-y -f --data_dir {params.database} ; '
        'emapper.py '
        '-i {input} '
        '--output {params.filename} '
        '-m {params.method} '
        '--data_dir {params.database} '
        '&> {log}'
     
rule trandecoder_trinitygg:
    input:
        'output/06_trinity_gg/Trinity.fasta'
    output:
        'output/08_transdecoder/Trinity.fasta.transdecoder.pep'
    params:
        transcripts = 'Trinity.fasta',
        dir = 'output/08_transdecoder',
        cd = '../../'
    log:
        'output/logs/08_transdecoder/transdecoder_trinity-gg.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        transdecoder
    shell:
        'cp {input} {params.dir} ; '
        'cd {params.dir} || exit 1 ; '        
        'TransDecoder.LongOrfs -t {params.transcripts} ; '
        'TransDecoder.Predict -t {params.transcripts} ; '
        'rm -rf {params.transcripts} ; '
        'cd {params.cd} ; '
        '&> {log}'

rule salmon_quant_trinitygg:
    input:
        r1 = 'output/01_cutadapt/trimmed_{sample}_r1.fq',
        r2 = 'output/01_cutadapt/trimmed_{sample}_r2.fq',
        index = 'output/07_salmon/index'
    output:
        'output/07_salmon/quant/{sample}/quant.sf'
    params:
        out_directory = 'output/07_salmon/quant/{sample}/'
    log:
        'output/logs/07_salmon/salmon_quant_{sample}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmon
    shell:
        'salmon quant '
        '--libType A '
        '--index {input.index} '
        '--mates1 {input.r1} '
        '--mates2 {input.r2} '
        '--output {params.out_directory} '
        '--threads {threads} '
        '--validateMappings '
        '--gcBias '
	'--seqBias '
        '&> {log}'

rule make_salmon_index_trinitygg:
    input:
        transcripts = 'output/06_trinity_gg/Trinity.fasta'
    output:
        directory('output/07_salmon/index')
    log:
        'output/logs/07_salmon/salmon_index.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmon
    shell:
        'salmon index '
        '--transcripts {input.transcripts} '
        '--index {output} '
        '--threads {threads} '
        '&> {log}'

rule trinitygg:
    input:
        bam = 'output/05_star/Aligned.sortedByCoord.out.bam',
        r1 = expand('data/reads/{sample}_r1.fq', sample=SAMPLES),
        r2 = expand('data/reads/{sample}_r2.fq', sample=SAMPLES)
    output:
        'output/06_trinity_gg/Trinity.fasta',
        'output/06_trinity_gg/Trinity.fasta.gene_trans_map'
    params:
        r1_line = lambda wildcards, input:
            ','.join(input.r1),
        r2_line = lambda wildcards, input:
            ','.join(input.r2),
        outdir = 'output/06_trinity'
    log:
        'output/logs/06_trinity_gg/trinity_gg.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        trinity
    shell:
        'Trinity '
        '--genome_guided_bam {input.bam} '
        '--genome_guided_max_intron 10000 '
        '--seqType fq '
        '--max_memory 500G '
        '--left {params.r1_line} '
        '--right {params.r2_line} '
        '--SS_lib_type RF '
        '--CPU {threads} '
        '--output {params.outdir} '
        '--trimmomatic '
        '&> {log}'

rule star_align:
    input:
        'output/05_star/index/genomeParameters.txt',
        r1 = expand('data/reads/{sample}_r1.fq', sample=SAMPLES),
        r2 = expand('data/reads/{sample}_r2.fq', sample=SAMPLES)
    output:
        'output/05_star/Aligned.sortedByCoord.out.bam'
    params:
        genome_index_dir = directory('output/05_star/index'),
        r1_line = lambda wildcards, input:
            ','.join(input.r1),
        r2_line = lambda wildcards, input:
            ','.join(input.r2),
        outfile = 'output/08_star/'
    threads:
         multiprocessing.cpu_count()
    log:
        'output/logs/05_star/star_align.log'
    singularity:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_index_dir} '
        '--readFilesIn {params.r1_line} {params.r2_line} '
        '--outFileNamePrefix {params.outfile} '
        '--outSAMtype BAM SortedByCoordinate '
        '&> {log}'

rule star_index:
    input: 
        genome = 'data/ref/ovalipes_catharus_genome.fasta',
        annotations = 'data/ref/ovalipes_catharus.gff3'
    output:
        'output/05_star/index/genomeParameters.txt'
    params:
        genome_index_dir = 'output/05_star/index'
    threads:
         multiprocessing.cpu_count()
    log:
        'output/logs/05_star/star_index.log'
    singularity:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.genome_index_dir} '
        '--genomeFastaFiles {input.genome} '
        '--sjdbGTFfile {input.annotations} '
        '--sjdbGTFtagExonParentTranscript Parent '
        '--sjdbOverhang 199 '
        '--genomeSAindexNbases 13 '
        '&> {log}'
          
rule deseq2:
    input:
        quant_files = expand('output/02_salmon/quant/{sample}/quant.sf', sample = SAMPLES),
        gff = 'data/ref/ovalipes_catharus.gff3',
        samples_file = 'output/samples.txt',
        salmon_dir = 'output/02_salmon',
        emapper_annotations_file = 'output/03_emapper/ovalipes_catharus.emapper.annotations'
    output:
        dds_file = 'output/04_deseq2/dds.Rds',
        deseq2_results_file = 'output/04_deseq2/ocath_deseq2_result.csv',
        deseq2_counts_file = 'output/04_deseq2/ocath_deseq2_counts.csv',
        deseq2_results_with_annotations_file = 'output/04_deseq2/ocath-genome_deseq2_result_wannot.csv',
        my_results_wtrans_ids = 'output/04_deseq2/wtrans.csv',
        tx2gene_file = 'output/04_deseq2/t2xgene.csv',
        my_results_dt_file = 'output/04_deseq2/results_dt.csv'
    log:
        'output/logs/04_deseq2/deseq2.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bioconductor
    script:
        'scripts/deseq2.R'

rule emapper:
    input:
        'data/ref/ovalipes_catharus.proteins.fa'
    output:
        'output/03_emapper/ovalipes_catharus.emapper.annotations'
    params:
        database = 'data/ref/',
        method = 'diamond',
        filename = 'output/03_emapper/ovalipes_catharus'
    log:
        'output/logs/03_emapper/emapper.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        emapper
    shell:
        'download_eggnog_data.py '
        '-y -f --data_dir {params.database} ; '
        'emapper.py '
        '-i {input} '
        '--output {params.filename} '
        '-m {params.method} '
        '--data_dir {params.database} '
        '&> {log}'
    
rule SalmonQuant:
    input:
        r1 = 'output/01_cutadapt/trimmed_{sample}_r1.fq',
        r2 = 'output/01_cutadapt/trimmed_{sample}_r2.fq',
        index = 'output/02_salmon/index'
    output:
        'output/02_salmon/quant/{sample}/quant.sf'
    params:
        out_directory = 'output/02_salmon/quant/{sample}/'
    log:
        'output/logs/02_salmon/salmon_quant_{sample}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmon
    shell:
        'salmon quant '
        '--libType A '
        '--index {input.index} '
        '--mates1 {input.r1} '
        '--mates2 {input.r2} '
        '--output {params.out_directory} '
        '--threads {threads} '
        '--validateMappings '
        '--gcBias '
	'--seqBias '
        '&> {log}'

rule makeSalmonIndex:
    input:
        gentrome = 'output/02_salmon/ref/gentrome.fa',
        decoys = 'output/02_salmon/ref/decoys.txt'
    output:
        directory('output/02_salmon/index')
    log:
        'output/logs/02_salmon/salmon_index.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmon
    shell:
        'salmon index '
        '--transcripts {input.gentrome} '
        '--index {output} '
        '--threads {threads} '
        '--decoys {input.decoys} '
        '&> {log}'

rule generateDecoyTrancriptome:
    input:
        genome = 'data/ref/ovalipes_catharus_genome.fasta',
        transcriptome = 'data/ref/ovalipes_catharus.mrna-transcripts.fa',
        annotation = 'data/ref/ovalipes_catharus.gff3'
    output:
        'output/02_salmon/ref/gentrome.fa',
        'output/02_salmon/ref/decoys.txt'
    params:
        outdir = 'output/02_salmon/ref'
    log:
        'output/logs/02_salmon/generateDecoyTrancriptome.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmontools
    shell:
        'generateDecoyTranscriptome.sh '
        '-g {input.genome} '
        '-t {input.transcriptome} '
        '-a {input.annotation} '
        '-o {params.outdir} '
	'-m /usr/local/bin/mashmap '
	'-b /usr/bin/bedtools '
        '-j {threads} '
        '&> {log}'

rule Cutadapt:
    input:
        r1 = 'data/reads/{sample}_r1.fq',
        r2 = 'data/reads/{sample}_r2.fq'
    output:
        r1 = 'output/01_cutadapt/trimmed_{sample}_r1.fq',
        r2 = 'output/01_cutadapt/trimmed_{sample}_r2.fq'
    log:
        'output/logs/01_cutadapt/{sample}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        cutadapt
    shell:
        'cutadapt '
        '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA '
        '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT '
        '-o {output.r1} '
        '-p {output.r2} '
        '{input.r1} '
        '{input.r2} '
        '--minimum-length 1 '
        '&> {log}'

