task genotypeTask{
    File reference_fasta
    File reference_fasta_idx
    File reference_dict
    File GVCF
    File GVCF_idx
    String sample_id

    runtime{
        docker : "us.gcr.io/broad-gatk/gatk:4.4.0.0"
        memory : "8 GB"
        cpu : "1"
        disks : "local-disk 1000 HDD"
    }


    command <<<
        /gatk/gatk --java-options  "-Xmx4g" GenotypeGVCFs --allow-old-rms-mapping-quality-annotation-data \
                                               -R ${reference_fasta} \
                                               -V ${GVCF} \
                                               -O ${sample_id}.genotyped.vcf.gz \
                                               -G StandardAnnotation \
                                               -G AS_StandardAnnotation
    >>>

    output{
        File output_vcf = "${sample_id}.genotyped.vcf.gz"
        File output_vcf_index = "${sample_id}.genotyped.vcf.gz.tbi"
    }
}

workflow genotypeGVCFs{
    File reference_fasta
    File reference_fasta_idx
    File reference_dict
    File GVCF
    File GVCF_idx
    String sample_id

    call genotypeTask{
    	input:
        	reference_fasta = reference_fasta,
            reference_fasta_idx = reference_fasta_idx,
            reference_dict = reference_dict,
            GVCF = GVCF,
            GVCF_idx = GVCF_idx,
            sample_id = sample_id
    }
    output{
        File output_vcf = genotypeTask.output_vcf
        File output_vcf_index = genotypeTask.output_vcf_index
    }
}
