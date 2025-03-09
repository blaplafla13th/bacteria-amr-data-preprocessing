params.RawReads = "/home/linuxbrew/workspace/fastq/*_{1,2}.fastq.gz"
params.ReferenceGenome = "/home/linuxbrew/workspace/sequence.fasta"
params.cpu = 12 // Change this to maximum number of CPUs to use
process TRIMMOMATIC {
    errorStrategy "finish"
    cpus params.cpu - 1
    tag "TRIMMOMATIC on $sample"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.paired.fq.gz"), emit: paired
    path("*_summary.txt"), emit: summary


    script:
    fq_1_paired = sample + "_1.paired.fq.gz"
    fq_1_unpaired = sample + "_1.unpaired.fq.gz"
    fq_2_paired = sample + "_2.paired.fq.gz"
    fq_2_unpaired = sample + "_2.unpaired.fq.gz"
    summary_files = sample + "_summary.txt"

    """
    trimmomatic PE -threads $task.cpus -summary ${summary_files} ${reads[0]} ${reads[1]} ${fq_1_paired} ${fq_1_unpaired} ${fq_2_paired} ${fq_2_unpaired} SLIDINGWINDOW:4:25 MINLEN:90
    """
}

process BWAMEM2 {
    errorStrategy "finish"
    cpus params.cpu - 1
    tag "BWA-MEM2 in $sample"

    input:
    tuple val(sample), path(paired_reads)

    output:
    tuple val(sample), path("*.sam")

    script:
    """
    bwa-mem2 index ${params.ReferenceGenome}
    bwa-mem2 mem -t ${task.cpus} ${params.ReferenceGenome} ${paired_reads[0]} ${paired_reads[1]} > ${sample}_alignment.sam
    """
}

process SAMTOBAM {
    errorStrategy "finish"
    cpus params.cpu / 4
    tag "SAM2BAM in $sample"

    input:
    tuple val(sample), path(sam_file)

    output:
    tuple val(sample), path("*.sorted.bam"), emit: bam

    script:
    """
    samtools view -bS ${sam_file} -o ${sample}_alignment.bam -@ ${task.cpus}
    samtools sort ${sample}_alignment.bam -o ${sample}_alignment.sorted.bam -@ ${task.cpus}
    samtools index -b -o ${sample}_index.bai
    """
}

process FREEBAYES {
    errorStrategy "finish"
    cpus 1
    tag "Freebayes in $sample"

    input:
    tuple val(sample), path(sorted_bam)

    output:
    tuple val(sample), path("*.vcf")

    script:
    // task properties
    """
    freebayes -f ${params.ReferenceGenome} -p 1 ${sorted_bam} > ${sample}.vcf
    """
}

process VCF_FILTER {
    errorStrategy "finish"
    cpus 1
    tag "Qualimap in $sample"

    input:
    tuple val(sample), path(vcf_files)

    output:
    tuple val(sample), path("*recode.vcf")

    script:
    """
    vcftools --vcf ${vcf_files} --out ${sample} --minDP 20 --minQ 50 --recode --remove-filtered-all --recode-INFO-all
    """
}

process VCF_GET_SNP_AND_RENAME {
    errorStrategy "finish"
    cpus 1
    tag "Rename $sample"

    input:
    tuple val(sample), path(vcf_file)

    output:
    path("*.vcf.gz*")

    script:
    """
    echo ${sample} > ${sample}
    bcftools view -v snps ${vcf_file} | bcftools reheader -s ${sample} - > ${sample}.vcf
    bgzip ${sample}.vcf
    tabix -p vcf ${sample}.vcf.gz
    """
}

process VCF_COMBINE {
    cpus params.cpu
    tag "Combine"
    publishDir "/home/linuxbrew/workspace/final", mode:"copy"

    input:
    val(vcf_gz_files)

    output:
    path("combined.bcf")

    script:
    """
    ulimit -n 10000
    vcf_list=\$(echo "${vcf_gz_files}" | sed 's/\\[\\|,\\|\\]//g')
    bcftools merge -0 -m all \$(du -sh \$vcf_list | grep -v tbi | sort -h | awk '{if (\$1 != "4.0K") print \$2}') \
    --threads ${task.cpus} | bcftools norm -d all -m -any --threads ${task.cpus} - | sed 's/0\\/0:/0:/g'  \
    | sed 's/1\\/1:/1:/g' | sed 's/\\.\\/\\.:/0:/g'| bcftools view -c 1 -o combined.bcf -
    """
}

process PLINK_QC {
    cpus 1
    tag "PLINK QC"
    publishDir "/home/linuxbrew/workspace/final", mode:"copy"

    input:
    path(bcf)
    output:
    path("plink.log")
    path("pairwise.plink.prune.in")
    path("bacteria.prune.plink.raw")

    script:
    """
    plink --bcf ${bcf} --maf 1e-3 --geno 0.2 --mind 0.2 --make-bed --out bacteria.plink --allow-extra-chr
    mv bacteria.plink.bim bacteria.plink.bim.backup
    awk '{print(\$1,\$4"."\$6"_"\$5,\$3,\$4,\$5,\$6)}' bacteria.plink.bim.backup > bacteria.plink.bim
    plink --bfile bacteria.plink --indep-pairwise 100 10 0.5 --out pairwise.plink --allow-extra-chr
    plink --bfile bacteria.plink --extract pairwise.plink.prune.in --recode A --out bacteria.prune.plink --allow-extra-chr
    cat bacteria.plink.log pairwise.plink.log > plink.log
    """
}

process TRIMMOMATIC_LOGS {
    cpus 1
    tag "TRIMMOMATIC LOGS"
    publishDir "/home/linuxbrew/workspace/final", mode:"copy"

    input:
    val(files)

    output:
    path("trimmomatic_total.log")

    script:
    """
    for file in \$(echo "${files}" | sed 's/\\[\\|,\\|\\]//g'); do
        sample_name=\$(basename \${file} | cut -d'_' -f1)
        echo \$sample_name
        cat \$file
    done >> trimmomatic_total.log
    """
}

workflow {
    Channel.fromFilePairs(params.RawReads, checkIfExists: true, type: 'file').set { RawReads_channel }
    trim_channel = TRIMMOMATIC(RawReads_channel)
    bwamem2_channel = BWAMEM2(TRIMMOMATIC.out.paired)
    TRIMMOMATIC_LOGS(TRIMMOMATIC.out.summary.collect())
    samtobam_channel = SAMTOBAM(bwamem2_channel)
    freebayes_channel = FREEBAYES(samtobam_channel)
    filter_channel = VCF_FILTER(freebayes_channel)
    rename_channel = VCF_GET_SNP_AND_RENAME(filter_channel)
    vcf_gz_files = rename_channel.collect()
    combine_channel = VCF_COMBINE(vcf_gz_files)
    plink_channel = PLINK_QC(combine_channel)
}

workflow.onComplete {

    println(workflow.success ?
    """
    Pipeline execution summary
    ------------------------------
    Completed at: ${workflow.complete}
    Duration: ${workflow.duration}
    Success: ${workflow.success}
    workDir: ${workflow.workDir}
    exit Status: ${workflow.exitStatus}
    """
    :
    """
    Failed: ${workflow.errorReport}
    Error message: ${workflow.errorMessage}
    exit status: ${workflow.exitStatus}
    """)
}
