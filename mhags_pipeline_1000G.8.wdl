import "https://api.firecloud.org/ga4gh/v1/tools/HLAthena_v1:HLAthena_v1_external/versions/8/plain-WDL/descriptor" as HLAthena

workflow modified_mHag {
  
  String sample_id
  String donor_output_basename
  String host_output_basename

  File donor_input_bam
  File donor_input_bai
  File host_input_bam
  File host_input_bai

  # hg19 ref files
  File ref_fasta
  File ref_fasta_index

  # BWA files
  File ref_dict

  # Single Nucleotide Polymorphism Database - hg19
  File dbSNP_vcf
  File dbSNP_vcf_index

  # Deletion/insertion polymorphism VCF - hg19 
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  # Deepvariant input 
  File? capture_interval_list
  Int runtime_disk = 1
  String runtime_docker

  # Functotator input
  String? data_sources_tar_gz

  # bmt input 
  String tissue
  String donor_sex
  String host_sex
  String gvhd_filter

  # HLAthena input
  File hlatypes_file
  Boolean aggregate_pep
  String peptide_col_name
  Boolean exists_expr
  String expr_col_name
  Array[String] lens
  Boolean exists_ctex

  # Disks and multipliers
  Int? additional_disk
  Int preemptible_tries
  Int agg_preemptible_tries

  Float cram_disk_multiplier = 8
  Float bwa_disk_multiplier = 2.5
  Int compression_level = 2
  Float md_disk_multiplier = 2.25
  Float sort_sam_disk_multiplier = 10


  Int disk_pad = select_first([additional_disk, 20])

  # Output basenames for different tasks 
  #String? donor_output_basename = basename(donor_input_cram, ".cram")
  #String? host_output_basename = basename(host_input_cram, ".cram")

  String? donor_recalibrated_bam_basename = donor_output_basename + ".aligned.duplicates_marked.recalibrated"
  String? host_recalibrated_bam_basename = host_output_basename + ".aligned.duplicates_marked.recalibrated"


  # File sizes 
  Float dbsnp_size = size(dbSNP_vcf, "GB") + size(dbSNP_vcf_index, "GB")



  call interval_list_to_bed {
      input:
          interval_list=capture_interval_list
  }

  call deep_variant as donor_deepvariant {
    input: 
      sample=donor_output_basename,
      capture_bed=interval_list_to_bed.bed,
      bam=donor_input_bam,
      bai=donor_input_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      runtime_docker = runtime_docker
  }

    call bgzip as donor_bgzip {
        input:
            sample=donor_output_basename,
            uncompressed_vcf=donor_deepvariant.vcf,
            runtime_disk = runtime_disk
    }

  call deep_variant as host_deepvariant {
    input: 
      sample=host_output_basename,
      capture_bed=interval_list_to_bed.bed,
      bam=host_input_bam,
      bai=host_input_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      runtime_docker = runtime_docker
  }

    call bgzip as host_bgzip {
        input:
            sample=host_output_basename,
            uncompressed_vcf=host_deepvariant.vcf,
            runtime_disk = runtime_disk
    }

    call bcf_filter {
            input:
              donor_input_vcf = donor_bgzip.filtered_vcf,
              host_input_vcf = host_bgzip.filtered_vcf
    }

    # Donor annotating 
    call Funcotate as FuncotateDonor {
            input:
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    input_vcf = bcf_filter.donor_vcf,
                    input_vcf_idx = bcf_filter.donor_vcf_tbi,
                    sample = donor_output_basename,

                    data_sources_tar_gz = data_sources_tar_gz
    }

    # Host annotating 
    call Funcotate as FuncotateHost {
            input:
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    input_vcf = bcf_filter.host_vcf,
                    input_vcf_idx = bcf_filter.host_vcf_tbi,
                    sample = host_output_basename,

                    data_sources_tar_gz = data_sources_tar_gz
    }

    call bmt {
            input: 
                    donor_input_vcf = FuncotateDonor.funcotated_output_file,
                    host_input_vcf = FuncotateHost.funcotated_output_file,
                    tissue = tissue,
                    donor_sex = donor_sex,
                    host_sex = host_sex,
                    gvhd_filter = gvhd_filter
    }

    #Array[File] discordances_file = select_first([bmt.specific_discordances, bmt.gvhd_specific_discordances])

    # Donor specific_discordances = tissue_specific_discordances.txt
    call blastFileGenerationForBlast as gvl_DonorblastFileGenerationForBlast {
            input:
                    specific_discordances = bmt.specific_discordances
    }

    # Creating donor custom proteome using makeblastdb 
    # donor custom proteome against donor peptides as query 
    call blastOperations as gvl_blastOperationsDonor {
            input: 
                    proteome_filename = bmt.donor_custom_proteome,
                    makeblastdb_out = "donor_db",
                    peptide_file = gvl_DonorblastFileGenerationForBlast.peptides_for_blastp, # for blastp
                    blastp_out = donor_output_basename+"_blastp_peptides_out.csv",
                    sample = donor_output_basename
    }

    # For donor : specific_discordances = tissue_specific_discordances.txt -> donor_discordances_after_blast.txt
    call blastOperationsPostProcessing as gvl_DonorBlastOperationsPostProcessing {
            input: 
                    specific_discordances = bmt.specific_discordances,
                    input_blastp_result = gvl_blastOperationsDonor.blastp_result,
                    sample = donor_output_basename,
                    sample_type = "donor",
                    condition="GvL",
                    known_minors=bmt.known_minors_bmt
    }


    # Donor discordances after blast as input 
    call blastFileGenerationForBlast as gvl_HostblastFileGenerationForBlast {
            input:
                    specific_discordances = gvl_DonorBlastOperationsPostProcessing.discordances_after_blast
    }

    # Creating host custom proteome using makeblastdb 
    # host custom proteome against host peptides as query 
    call blastOperations as gvl_blastOperationsHost {
            input: 
                    proteome_filename = bmt.host_custom_proteome,
                    makeblastdb_out = "host_db",
                    peptide_file = gvl_HostblastFileGenerationForBlast.peptides_for_blastp, # for blastp 
                    blastp_out = host_output_basename+"_blastp_peptides_out.csv",
                    sample = host_output_basename
    }

    # For host : specific_discordances = donor_discordances_after_blast.txt with genes to remove -> host_discordances_after_blast.txt
    call blastOperationsPostProcessing as gvl_HostBlastOperationsPostProcessing {
            input: 
                    specific_discordances = gvl_DonorBlastOperationsPostProcessing.discordances_after_blast,
                    input_blastp_result = gvl_blastOperationsHost.blastp_result,
                    sample_type = "host",
                    sample = host_output_basename,
                    condition="GvL",
                    known_minors=bmt.known_minors_bmt
    }

    call HLAthena_preprocessing as gvl_HLAthena_preprocessing {
            input:
                    discordances_after_blast = gvl_HostBlastOperationsPostProcessing.allgenes_pumas
    }

    call HLAthena.parse_alleles_task {
            input : 
                    preemptible = 3,
                    memoryGB = 8,
                    diskGB = 10,

                    alleles_file = hlatypes_file
    }


    Array[String] gvl_alleles_parsed = read_lines(parse_alleles_task.alleles_parsed)
    if (aggregate_pep) {
            call HLAthena.aggregate_pep_task as gvl_aggregate_pep_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            peptide_list = gvl_HLAthena_preprocessing.HLAthena_preprocessed_peptides,
                            peptide_col_name = peptide_col_name,
                            exists_expr = exists_expr,
                            expr_col_name = expr_col_name
            }
    }

    ### Split merged peptides by length
    call HLAthena.split_peptides_len_task as gvl_split_peptides_len_task {
            input:
                preemptible = 3,
                memoryGB = 8,
                diskGB = 10,

                peptide_list = gvl_aggregate_pep_task.peptide_list_agg,
                peptide_col_name = peptide_col_name,
                lens = lens
    }

    Array[String] gvl_lens_present = read_lines(gvl_split_peptides_len_task.lens_present)

    ### Generate peptide sequences features: dummy encoding only 
    ### (blosum and fuzzy generated on demand from the dummy encoding upon prediction)
    scatter (pepfilelen in gvl_split_peptides_len_task.pep_len_files) {
            call HLAthena.featurize_encoding_task as gvl_featurize_encoding_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            peptide_len_list = pepfilelen,
                            encoding = "dummy",
                            peptide_col_name = peptide_col_name
            }
    }

    ### Predict with MS models (now includes ranks)
    scatter (allele_featfile in cross(gvl_alleles_parsed, gvl_featurize_encoding_task.feature_files)) {
            call HLAthena.predict_ms_allele_task as gvl_predict_ms_allele_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            patient = sample_id,
                            featfile = allele_featfile.right,
                            lens = gvl_lens_present,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            allele = allele_featfile.left,
                            models_path = "gs://msmodels/",
                            peptide_col_name = peptide_col_name,
                            features_all = "features_AAPos_AAPCA_LogTPM_CNN_Kidera_Gene",
                            feature_sets = "features_AAPos_AAPCA_Kidera"
            }
    }

    ### Merge predictions into ensemble model scores: MSIntrinsic, MSIntrinsicC, MSIntrinsicEC; optinally merge NetMHC scores
    scatter (len in gvl_lens_present) {
            call HLAthena.merge_mspreds_task as gvl_merge_mspreds_task {
                    input:
                            preemptible = 3,
                            memoryGB_merge = 16,
                            #diskGB_merge = 20, - PRESENT IN THE TASK BUT ERROR THAT IT IS NOT PRESENT 
                            diskGB = 10,
                            
                            patient = sample_id,
                            alleles = gvl_alleles_parsed,
                            len = len,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            mspred_files = gvl_predict_ms_allele_task.mspred_files,
                            peptide_col_name = peptide_col_name
                            #run_netmhc = run_netmhc,
                            #netmhc_EL = predict_netmhc_task.netmhc_EL,
                            #netmhc_BA = predict_netmhc_task.netmhc_BA
            }
    }

    ### Convert scores to ranks
    scatter (preds_file in gvl_merge_mspreds_task.mspreds_file_wide) {
            call HLAthena.get_ranks_task as gvl_get_ranks_task {
                    input:
                            preemptible = 3,
                            memoryGB_ranks = 16,
                            diskGB = 10,

                            patient = sample_id,
                            alleles = gvl_alleles_parsed,
                            lens = gvl_lens_present,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            #run_netmhc = run_netmhc,
                            models_path = "gs://msmodels/",
                            peptide_col_name = peptide_col_name,
                            preds_file = preds_file
            }
    }

    ### Concatenate ranks files for all lengths
    call HLAthena.ranks_concat_lens_task as gvl_ranks_concat_lens_task {
            input:
                    preemptible = 3,
                    memoryGB_assign = 16,
                    #diskGB_ranks_concat = 20,
                    diskGB = 10,

                    patient = sample_id,
                    peptide_col_name = peptide_col_name,
                    exists_ctex = exists_ctex,
                    exists_expr = exists_expr,
                    #run_netmhc = run_netmhc,
                    assign_by_ranks_or_scores = "ranks",
                    assign_threshold = "0.5",
                    assign_colors = "#8c0a82, #ab0c9f, #048edb, #00a0fa, #5F36A0, #ff6000, darkorange, #e05400, grey",
                    ranks_len_files = gvl_get_ranks_task.ranks_file
    }

    call HLAthena_postprocessing as gvl_HLAthena_postprocessing {
            input :
                    sample_predictions = gvl_ranks_concat_lens_task.sample_predictions,
                    discordances_after_blast = gvl_HostBlastOperationsPostProcessing.allgenes_pumas,
                    tissue = tissue,
                    sample_name = sample_id + "_GvL"
    }





    call merge_summary_csv_files {
      input:
        bcf_dict = bcf_filter.bcf_summary,
        bmt_gvl_dict = bmt.GvL_dict_nums,
        gvl_postBlast_dict = gvl_HostBlastOperationsPostProcessing.GvL_postBlast_dict_nums,
        gvl_postHLAthena_dict = gvl_HLAthena_postprocessing.GvL_postHLAthena_dict_nums,
#        bmt_gvhd_dict = bmt.GvHD_dict_nums,
#        gvhd_postBlast_dict = gvhd_HostBlastOperationsPostProcessing.GvHD_postBlast_dict_nums,
#        gvhd_postHLAthena_dict = gvhd_HLAthena_postprocessing.GvHD_postHLAthena_dict_nums

    }

    output {
      # BMT OUTPUT 
      File specific_discordances_file = bmt.specific_discordances 
#      File gvhd_specific_discordances_file = bmt.gvhd_specific_discordances

      # HLAthena OUTPUT
      File gvl_sample_predictions_file = gvl_ranks_concat_lens_task.sample_predictions
#      File gvhd_sample_predictions_file = gvhd_ranks_concat_lens_task.sample_predictions

      # HLAthena POST PROCESSING OUTPUT
      File GvL_postHLAthena_dict_nums_file = gvl_HLAthena_postprocessing.GvL_postHLAthena_dict_nums
      File gvl_binding_putativeMinorAntigens_file = gvl_HLAthena_postprocessing.binding_putativeMinorAntigens
#      File gvhd_binding_putativeMinorAntigens_file = gvhd_HLAthena_postprocessing.binding_putativeMinorAntigens

      # FINAL OUTPUT
      File pipeline_run_summary_csv = merge_summary_csv_files.output_summary_file_result
      File pipeline_run_summary_csv_2 = merge_summary_csv_files.output_summary_file_result_2
    }
}




# Interval list (coding regions of hg19) is converted into bed format 
task interval_list_to_bed {
    File interval_list
    String bed_path = sub(basename(interval_list), 'interval_list', 'bed')

    command <<<
    set -xeuo pipefail

    # interval lists have headers that need to be removed and are 1-indexed
    # see also https://www.biostars.org/p/84686/
    grep -v '^@' ${interval_list} \
    | awk -v OFS='\t' '{print $1, $2 - 1, $3}' \
    | sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge \
    > ${bed_path}
    >>>

    output {
        File bed = '${bed_path}'
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk 1 HDD'
        preemptible: 3
        docker: 'quay.io/biocontainers/bedtools:2.28.0--hdf88d34_0'
    }
}



# Germline Variant calling 
task deep_variant {
    String sample
    File bam
    File bai
    String model_type = 'WES'
    File ref_fasta
    File ref_fasta_index
    File? capture_bed

    Int runtime_cpus = 32 
    String runtime_docker
    Int runtime_memory = ceil(1.1 * runtime_cpus)
    Int runtime_disk_buffer = 3
    Int runtime_disk = ceil(1.15 * (size(ref_fasta, 'G') + size(bam, 'G')) + runtime_disk_buffer)
    Int runtime_preemptible = 3
    Int resource_log_interval = 10 

    command <<<
    # log resource usage for debugging purposes
    function runtimeInfo() {
        echo [$(date)]
        echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
    }
    while true;
        do runtimeInfo >> resource.log;
        sleep ${resource_log_interval};
    done &
    lscpu

    set -xeuo pipefail

    # make symbolic links to ensure BAM and index are in expected structure even after localization
    ln -s ${bai} reads.bai
    ln -s ${bam} reads.bam

    # make symbolic links to ensure reference and index are in expected structure even after localization
    ln -s ${ref_fasta} reference.fa
    ln -s ${ref_fasta_index} reference.fa.fai

    mkdir deepvariant_tmp

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=${model_type} \
        --ref=reference.fa \
        --reads=reads.bam \
        --regions=${capture_bed} \
        --intermediate_results_dir=deepvariant_tmp \
        --output_vcf=${sample}.vcf \
        --num_shards=${runtime_cpus}
    >>>

    output {
        File vcf = '${sample}.vcf'
        File resource_log = 'resource.log'
    }
    
    runtime {
        memory: '${runtime_memory} GB'
        cpu: '${runtime_cpus}'
        disks: 'local-disk ${runtime_disk} SSD'
        preemptible: '${runtime_preemptible}'
        docker: '${runtime_docker}'
    }
}


# gunzips the VCF files from deepvariant and contains only "PASS" variants 
task bgzip {
    String sample
    File uncompressed_vcf
    Int runtime_disk

    command <<<
    set -xeuo pipefail

    bcftools view -Oz -o ${sample}.vcf.gz ${uncompressed_vcf}
    bcftools index --tbi ${sample}.vcf.gz

    # create version of VCF with only PASSing variants
    bcftools view -Oz -o ${sample}_filtered.vcf.gz -f PASS ${uncompressed_vcf}
    bcftools index --tbi ${sample}_filtered.vcf.gz
    >>>

    output {
        File vcf = '${sample}.vcf.gz'
        File vcf_index = '${sample}.vcf.gz.tbi'
        File filtered_vcf = '${sample}_filtered.vcf.gz'
        File filtered_vcf_index = '${sample}_filtered.vcf.gz.tbi'
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk ${runtime_disk} HDD'
        preemptible: 3
        docker: 'quay.io/biocontainers/bcftools:1.9--ha228f0b_3'
    }
}


# Filters "PASS" and "MIS" variants 
task bcf_filter {

  # Filtering the Deepvariant VCF for donor and host files to remore ref_calls and mis 

        File donor_input_vcf
        File host_input_vcf

        command {
                bcftools filter ${donor_input_vcf} -i 'FILTER="PASS" && GT!="mis"' > donor_bcf_filtered.vcf.gz -O z 
                tabix -p vcf donor_bcf_filtered.vcf.gz

                bcftools filter ${host_input_vcf} -i 'FILTER="PASS" && GT!="mis"' > host_bcf_filtered.vcf.gz -O z
                tabix -p vcf host_bcf_filtered.vcf.gz

                donor_bcf_num=$(bcftools filter ${donor_input_vcf} -i 'FILTER="PASS" && GT!="mis"' | grep -c '^[^#]')
                host_bcf_num=$(bcftools filter ${host_input_vcf} -i 'FILTER="PASS" && GT!="mis"' | grep -c '^[^#]')

                echo "host_VCF_rows","donor_VCF_rows" > bcf_summary.csv
                echo "$host_bcf_num","$donor_bcf_num" >> bcf_summary.csv

                rm -r ${donor_input_vcf}
                rm -r ${host_input_vcf}             
        }

        output {
                File donor_vcf="donor_bcf_filtered.vcf.gz"
                File donor_vcf_tbi="donor_bcf_filtered.vcf.gz.tbi"

                File host_vcf="host_bcf_filtered.vcf.gz"
                File host_vcf_tbi="host_bcf_filtered.vcf.gz.tbi"

                File bcf_summary="bcf_summary.csv"
        }

        runtime {
          docker : "dceoy/bcftools"
        }

        meta {
                author: "Nidhi Hookeri"
                email: "nhookeri@broadinstitute.org"
                description: "Filters the vcf file."
        }
}


# Functional Annotator 
task Funcotate {
    
  # Annotating the BCF filtered VCF files for donor and host

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File input_vcf
        File input_vcf_idx
        #String reference_version
        #Boolean compress
        String sample

        #Int num_threads

        File? data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124g.tar.gz"
        #String? transcript_selection_mode = "ALL"

        String output_vcf = sample + "_funcotated" + ".vcf"
        String output_vcf_idx = output_vcf + ".idx"
        
        command <<<
                # Extract the tar.gz:
                echo "Extracting data sources tar/gzip file..."
                mkdir datasources_dir
                tar zxvf ${data_sources_tar_gz} -C datasources_dir --strip-components 1
                DATA_SOURCES_FOLDER="$PWD/datasources_dir"

                # Annotating vcf file
                gatk --java-options "-Xmx2048m" Funcotator \
                     --data-sources-path $DATA_SOURCES_FOLDER \
                     --ref-version "hg19" \
                     --output-file-format "VCF" \
                     -R ${ref_fasta} \
                     -V ${input_vcf} \
                     -O ${output_vcf} \
                     --transcript-selection-mode "ALL" \
                     --create-output-variant-index
                bgzip ${output_vcf}

                rm -r ${input_vcf}
                rm -r ${input_vcf_idx}
                rm -r ${ref_fasta}
                rm -r ${ref_fasta_index}
                rm -r ${ref_dict}
                rm -r ${data_sources_tar_gz}
        >>>

        runtime {
                docker: "us.gcr.io/broad-gatk/gatk:4.1.4.0"
                cpu: 32
        }

        output {
                File funcotated_output_file = "${output_vcf}.gz"
                File funcotated_output_file_index = "${output_vcf_idx}"
        }
}



task bmt {
        File donor_input_vcf
        File host_input_vcf
        String tissue
        String donor_sex
        String host_sex
        String gvhd_filter

        command {
                pwd
                df
                python3 /pipeline/bmt-simulation.py -donor_vcf ${donor_input_vcf} -host_vcf ${host_input_vcf} -donor_sex ${donor_sex} -host_sex ${host_sex} -tissue ${tissue} -gvhd ${gvhd_filter}
        }

        runtime {
                docker : "nidhihookeri/1000g_mhags:latest"
                cpu: 32
        }
        output {
                String result = stdout()
                #File donor_annotation = "donor_annotations.txt"
                #File donor_not_mutated = "donor_not_mutated.txt"

                #File host_annotation = "host_annotations.txt"
                #File host_not_mutated = "host_not_mutated.txt"

                #File donor_annotations_withprot = "donor_annotations_withprot.txt"
                #File donor_mutated_proteins = "donor_mutated_proteins.txt"
                #File donor_mutation_indices = "donor_mutation_indices.txt"

                #File host_annotations_withprot = "host_annotations_withprot.txt"
                #File host_mutated_proteins = "host_mutated_proteins.txt"
                #File host_mutation_indices = "host_mutation_indices.txt"

                #File donor_and_host_data_merged = "donor_and_host_data_merged.txt"

                File donor_custom_proteome = "donor_custom_proteome.fasta"
                File host_custom_proteome = "host_custom_proteome.fasta"

                #File discordances_after_bmt_simulation = "discordances_after_bmt_simulation.txt"

                File all_discordances = "all_discordances.txt"

                File specific_discordances = tissue+"_specific_discordances.txt"
                File? gvhd_specific_discordances = "GvHD_specific_discordances.txt"

                File GvL_dict_nums = "GvL_dict_nums.csv"
                File? GvHD_dict_nums = "GvHD_dict_nums.csv"

                File known_minors_bmt = "known_minors_info.csv"

        }
}


task blastFileGenerationForBlast {
        File specific_discordances

        command {
                # Python file outputs a fa file that is used in the blastOperations task
                # Function that filters the peptide_fasta_id and peptide (pep) from the tissue specific discordances file

                # 2 different specific_discordances are used - mentioned above in the call part 
                python3 /pipeline/blast_preprocessing.py -specific_discordances ${specific_discordances} 
        }

        output {
                File peptides_for_blastp="peptides_for_blast.fa"
        }

        runtime {
                docker : "nidhihookeri/known_minors_pipeline"
                cpu : 16
        }

        meta {
                author: "Nidhi Hookeri"
                email: "nhookeri@broadinstitute.org"
        }
}


task blastOperations {
        File proteome_filename
        String makeblastdb_out
        String blastp_out
        File peptide_file
        String sample

        command {
                # Creates a blast database for the donor sample custom proteome from the task bmt
                # Important to create the database and the run blastp in the same folder as result 
                makeblastdb -in ${proteome_filename} -dbtype prot -input_type fasta -out ${makeblastdb_out}
                # custome proteome db against tissue specific discordances peptide 
                blastp -query ${peptide_file} -db ${makeblastdb_out} -out ${blastp_out} -outfmt "10 sseqid qseqid slen qlen sstart send qstart qend sseq qseq evalue score length pident nident"
        }

        output {
                File phr="${makeblastdb_out}.phr"
                File pin="${makeblastdb_out}.pin"
                File psq="${makeblastdb_out}.psq"
                File blastp_result=sample+"_blastp_peptides_out.csv"
        }

        runtime {
                docker : "biocontainers/blast:v2.2.31_cv2"
                cpu : 16
        }
}


task blastOperationsPostProcessing {
        File specific_discordances
        File input_blastp_result 
        String sample_type
        String sample
        String condition
        File known_minors

        command {
                # Assigning column names to the specific_discordances + filter_blast_search 
                # For donor : specific_discordances = tissue_specific_discordances.txt -> donor_discordances_after_blast.txt
                # For host : specific_discordances = donor_discordances_after_blast.txt with genes to remove -> host_discordances_after_blast.txt
                python3 /pipeline/blast_postprocessing.py -specific_discordances ${specific_discordances} -blastp_result ${input_blastp_result} -sample ${sample} -sample_type ${sample_type} -condition ${condition} -known_minors ${known_minors}
        }

        output {
                File discordances_after_blast=sample+"_discordances_after_blast.txt"
                File? genes_with_minors="genes_with_minors.txt"
                File? allgenes_pumas="allgenes_pumas.txt"
                File? GvL_postBlast_dict_nums="GvL_postBlast_dict_nums.csv"
                File? GvHD_postBlast_dict_nums="GvHD_postBlast_dict_nums.csv"
                File? known_minors_blast="known_minors_info.csv"
        }

        runtime {
                docker :  "nidhihookeri/known_minors_pipeline"
                cpu : 16
                disks: "local-disk 15 HDD"
        }

        meta {
                author: "Nidhi Hookeri"
                email: "nhookeri@broadinstitute.org"
        }
}



task HLAthena_preprocessing {
        File discordances_after_blast 

        command <<<
                python3 /pipeline/HLAthena_preprocessing.py -discordances_after_blast ${discordances_after_blast}
        >>>

        runtime {
                docker :  "nidhihookeri/known_minors_pipeline"
        }

        output {
                File HLAthena_preprocessed_peptides="HLAthena_peptides.txt"
        }
}



task HLAthena_postprocessing {
        File sample_predictions
        File discordances_after_blast
        String tissue
        String sample_name

        command <<<
                python3 /pipeline/HLAthena_postprocessing.py -hlathena_predictions ${sample_predictions} -discordances ${discordances_after_blast} -tissue_type ${tissue} -sample_name ${sample_name}
        >>>
        
        output {
                File binding_putativeMinorAntigens=sample_name+"_binding_putativeMinorAntigens.txt"
                File GvL_postHLAthena_dict_nums="GvL_postHLAthena_dict_nums.csv"
                File? GvHD_postHLAthena_dict_nums="GvHD_postHLAthena_dict_nums.txt"
        }

        runtime{
                docker : "nidhihookeri/known_minors_pipeline"
                cpu : 32
        }
}


task merge_summary_csv_files {
    File bcf_dict
    File bmt_gvl_dict
#    File bmt_gvhd_dict
    File gvl_postBlast_dict
#    File gvhd_postBlast_dict
    File gvl_postHLAthena_dict
#    File gvhd_postHLAthena_dict

    String sample_id
    String output_summary_file = sample_id + ".csv"
    String output_summary_file_2 = sample_id + "_format2" + ".csv"

    command <<<
    paste -d, ${bcf_dict} ${bmt_gvl_dict} ${gvl_postBlast_dict} ${gvl_postHLAthena_dict}  >${output_summary_file}

    cat ${bcf_dict} ${bmt_gvl_dict} ${gvl_postBlast_dict} ${gvl_postHLAthena_dict}  > ${output_summary_file_2}
    >>>

    runtime {
      docker : "ubuntu:latest"
    }

    output {
      File output_summary_file_result="${output_summary_file}"
      File output_summary_file_result_2="${output_summary_file_2}"
    }
}