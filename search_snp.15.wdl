workflow search_snp {

        File gVCF
        File gVCF_index
        Int SNP_pos
        Int chromosome

        String sample_id
		Int? disk_size_manual

        Int disk_size = select_first([disk_size_manual, ceil(size(gVCF, 'GB')*10)])


    call search_SNP {
        input:
            gz_vcf = gVCF,
            gz_vcf_index = gVCF_index,
            SNP_pos = SNP_pos,
            sample_id = sample_id,
            disk_size = disk_size,
            chromosome = chromosome

    }
    output{
    Int presence_SNP=search_SNP.snp_search_result
    }
}


    task search_SNP {

            File gz_vcf
            File gz_vcf_index
            Int chromosome
            Int SNP_pos
            String sample_id
            Int disk_size


    command <<<
            bcftools view -r chr${chromosome}:${SNP_pos} ${gz_vcf} |grep -v ^#| wc -l   
	>>>

        runtime {
      docker: "quay.io/biocontainers/bcftools:1.9--ha228f0b_3"
      disks: "local-disk " + disk_size + " HDD"

        }

    output {
        Int snp_search_result = read_int(stdout())

    }


}
