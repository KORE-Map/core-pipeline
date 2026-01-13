#####################################################################
### paths of RNAseq softwares

# The main variable must include the directory to be analyzed.  
# The directory specified below is an example.
main="""directory name"""
fq=${main}/fastq
tout=${main}/trimmed
qout=${main}/QC
bout=${main}/bam
eout=${main}/expression
mkdir -p ${qout}
mkdir -p ${bout}
mkdir -p ${eout}
mkdir -p ${tout}

# The following variable assignment must include each user's own index directory. 
# The variable shown below is an example.
qc_dir=/data/utils/FastQC
tg_dir=/data/utils/TrimGalore-0.6.6
rsem_dir=/data/utils/RSEM-1.3.3
star_dir=/data/utils/STAR-2.7.9a/bin/Linux_x86_64
gf=/data/Genome/REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=/data/Genome/REF/Homo_sapiens.GRCh38.84.gtf
rindex=/data/Genome/INDEX/hs/RSEM
sindex=/data/Genome/INDEX/hs/STAR99
#####################################################################
ls -a $fq | cut -d '.' -f 1 | sort -u > $main/fnameC.txt
#####################################################################
while read name
do
    echo $name
    ${tg_dir}/trim_galore \
    --adapter AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
    --phred 33 --quality 20 --length 50 --paired \
    --cores 6 --gzip --output_dir ${tout} ${fq}/${name}.R1.fq.gz ${fq}/${name}.R2.fq.gz
    echo $name
    ${star_dir}/STAR \
    --genomeDir ${sindex} \
    --readFilesIn ${tout}/${name}.R1_val_1.fq.gz ${tout}/${name}.R2_val_2.fq.gz \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --twopassMode Basic \
    --runThreadN 80 \
    --sjdbGTFfile ${gtf} \
    --genomeLoad NoSharedMemory \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped None \
    --outFileNamePrefix ${bout}/${name}.STAR.

    ${rsem_dir}/rsem-calculate-expression --time --num-threads 20 --alignments --paired-end --no-bam-output \
    ${bout}/${name}.STAR.Aligned.toTranscriptome.out.bam ${rindex}/index ${eout}/${name}.STAR.RSEM
done < $main/fnameC.txt

# The main and qc_dir variable must include the directory to be analyzed.  
# The directory specified below is an example.
qc_dir=/data/utils/FastQC
main="""directory name"""

fq=${main}/fastq
qout=${main}/QC

for file in $fq/*.fq.gz
do
    echo ${file}
    ${qc_dir}/fastqc -f fastq -o ${qout} ${file}
done