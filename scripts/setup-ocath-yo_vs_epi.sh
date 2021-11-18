#!/bin/sh

########################################
#Ovalipes catharus Y-organ Vs epidermis#
########################################

#log#
exec &> setup_ocath-yo_vs_epi.log

mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-yo_vs_epi/data/reads
reads_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-yo_vs_epi/data/reads

for (( i = 19; i <= 28; i++ ))
  do
cp /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/OG6267P2/HMJCYBCX3-6267-$i-49-01_S*_L00*_R1_001.fastq.gz $reads_dir
cp /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/OG6267P2/HMJCYBCX3-6267-$i-49-01_S*_L00*_R2_001.fastq.gz $reads_dir
cp /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/OG6267P2/HMHW2BCX3-6267-$i-49-01_S*_L00*_R1_001.fastq.gz $reads_dir
cp /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/OG6267P2/HMHW2BCX3-6267-$i-49-01_S*_L00*_R2_001.fastq.gz $reads_dir
done

zcat $reads_dir/HMJCYBCX3-6267-19-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-19-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/YO_1_r1.fq
zcat $reads_dir/HMJCYBCX3-6267-20-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-20-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/YO_2_r1.fq
zcat $reads_dir/HMJCYBCX3-6267-21-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-21-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/YO_3_r1.fq
zcat $reads_dir/HMJCYBCX3-6267-22-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-22-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/YO_4_r1.fq
zcat $reads_dir/HMJCYBCX3-6267-23-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-23-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/YO_5_r1.fq

zcat $reads_dir/HMJCYBCX3-6267-24-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-24-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/Epi_1_r1.fq
zcat $reads_dir/HMJCYBCX3-6267-25-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-25-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/Epi_2_r1.fq
zcat $reads_dir/HMJCYBCX3-6267-26-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-26-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/Epi_3_r1.fq
zcat $reads_dir/HMJCYBCX3-6267-27-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-27-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/Epi_4_r1.fq
zcat $reads_dir/HMJCYBCX3-6267-28-49-01_S*_L00*_R1_001.fastq.gz $reads_dir/HMHW2BCX3-6267-28-49-01_S*_L00*_R1_001.fastq.gz > $reads_dir/Epi_5_r1.fq

zcat $reads_dir/HMJCYBCX3-6267-19-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-19-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/YO_1_r2.fq
zcat $reads_dir/HMJCYBCX3-6267-20-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-20-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/YO_2_r2.fq
zcat $reads_dir/HMJCYBCX3-6267-21-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-21-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/YO_3_r2.fq
zcat $reads_dir/HMJCYBCX3-6267-22-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-22-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/YO_4_r2.fq
zcat $reads_dir/HMJCYBCX3-6267-23-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-23-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/YO_5_r2.fq

zcat $reads_dir/HMJCYBCX3-6267-24-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-24-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/Epi_1_r2.fq
zcat $reads_dir/HMJCYBCX3-6267-25-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-25-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/Epi_2_r2.fq
zcat $reads_dir/HMJCYBCX3-6267-26-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-26-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/Epi_3_r2.fq
zcat $reads_dir/HMJCYBCX3-6267-27-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-27-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/Epi_4_r2.fq
zcat $reads_dir/HMJCYBCX3-6267-28-49-01_S*_L00*_R2_001.fastq.gz $reads_dir/HMHW2BCX3-6267-28-49-01_S*_L00*_R2_001.fastq.gz > $reads_dir/Epi_5_r2.fq

#make sample.txt#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-yo_vs_epi/output/02_salmon
output_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-yo_vs_epi/output
echo "
sample	group
YO_1	YO
YO_2	YO
YO_3	YO
YO_4	YO
YO_5	YO
Epi_1	Epi
Epi_2	Epi	
Epi_3	Epi
Epi_4	Epi
Epi_5	Epi" > $output_dir/samples.txt


mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-yo_vs_epi/data/ref
ref_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-yo_vs_epi/data/ref
cp /Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-annotate/data/ovalipes_catharus_genome.fasta $ref_dir/ovalipes_catharus_genome.fasta
cp /Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-annotate/output/01_funannotate/annotate_results/ovalipes_catharus.gff3 $ref_dir/ovalipes_catharus.gff3
cp /Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-annotate/output/01_funannotate/annotate_results/ovalipes_catharus.mrna-transcripts.fa $ref_dir/ovalipes_catharus.mrna-transcripts.fa
cp /Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-annotate/output/01_funannotate/annotate_results/ovalipes_catharus.proteins.fa $ref_dir/ovalipes_catharus.proteins.fa


#make virtual environment in which to run snakemake#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-yo_vs_epi/yo_vs_epi_venv
venv_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ocath_genome/ocath-yo_vs_epi/yo_vs_epi_venv
python3 -m venv $venv_dir
source $venv_dir/bin/activate
pip3 install --upgrade pip
pip3 install snakemake
deactivate
 
