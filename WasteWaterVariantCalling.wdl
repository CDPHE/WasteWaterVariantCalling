version 1.0

workflow WasteWaterVariantCalling {

    input {
        Array[File] sorted_bam
        File covid_genome
        File voc_bed
        File voc_annotations
        Array[String] sample_id
        String out_dir
    }

    scatter (id_bam in zip(sample_id, sorted_bam)) {
        call add_RG {
            input:
                sample_id = id_bam.left,
                bam = id_bam.right
        }
        call variant_calling {
            input:
                bam = add_RG.rgbam,
                ref = covid_genome,
                sample_id = id_bam.left
        }
        call sort_vcf {
            input:
                vcf = variant_calling.vcf,
                sample_id = id_bam.left
        }
        call sample_VOCs {
            input:
                vcf = sort_vcf.sorted_vcf,
                bed = voc_bed,
                sample_id = id_bam.left
        }
        call vcf2tsv {
            input:
                vcf = sample_voc.sample_voc_vcf,
                sample_id = id_bam.left,
                bed = voc_bed
        }
        call fill_NA {
            input:
                tsv = vcf2tsv.sample_voc_tsv,
                sample_id = id_bam.left,
                voc_bed = voc_bed
        }
        call allele_freq {
            input:
                tsv = fill_NA.fill_NA_tsv,
                sample_id = id_bam.left
        }
        call reformat_tsv {
            input:
                tsv = allele_freq.allele_freq_tsv,
                sample_id = id_bam.left
        }
        call summary_prep {
            input:
                tsv = reformat_tsv.reformat_tsv_tsv,
                sample_id = id_bam.left,
                voc_annotations = voc_annotations
        }
    }
    call dashboard_tsv {
        input:
            tsv = summary_prep.sample_voc_tsv_summary,
            tsv_dash = summary_prep.sample_voc_tsv_dash,
            tsv_counts = summary_prep.sample_voc_tsv_counts,
            voc_annotations = voc_annotations
    }
    call summary_tsv {
        input:
            tsv = dashboard_tsv.voc_summary_temp
    }

    call transfer_outputs {
        input:
            variants = variant_calling.vcf,
            sorted_vcf = sort_vcf.sorted_vcf,
            sample_voc_vcf = sample_voc.sample_voc_vcf,
            sample_voc_tsv = vcf2tsv.sample_voc_tsv,
            sample_voc_tsv_summary = summary_prep.sample_voc_tsv_summary,
            sample_voc_tsv_dash = summary_prep.sample_voc_tsv_dash,
            sample_voc_tsv_counts = summary_prep.sample_voc_tsv_counts,
            voc_summary = summary_tsv.voc_summary,
            voc_dashboard = dashboard_tsv.voc_dashboard,
            voc_counts = dashboard_tsv.voc_counts,
            out_dir = out_dir
    }
    
    output {
        Array[File] addrg_bam = add_RG.rgbam
        Array[File] variants = variant_calling.vcf
        Array[File] sorted_vcf = sort_vcf.sorted_vcf
        Array[File] sample_voc_vcf = sample_voc.sample_voc_vcf
        Array[File] sample_voc_tsv = vcf2tsv.sample_voc_tsv
        Array[File] sample_voc_tsv_summary = summary_prep.sample_voc_tsv_summary
        Array[File] sample_voc_tsv_dash = summary_prep.sample_voc_tsv_dash
        Array[File] fill_NA_tsv = fill_NA.fill_NA_tsv
        Array[File] allele_freq_tsv = allele_freq.allele_freq_tsv
        Array[File] reformat_tsv_tsv = reformat_tsv.reformat_tsv_tsv
        Array[File] sample_voc_tsv_counts = summary_prep.sample_voc_tsv_counts
        File voc_summary_temp = dashboard_tsv.voc_summary_temp
        File voc_summary = summary_tsv.voc_summary
        File voc_dashboard = dashboard_tsv.voc_dashboard
        File voc_counts = dashboard_tsv.voc_counts
        String transfer_date = transfer_outputs.transfer_date
        
    }
}

task add_RG {
    input {
        String sample_id
        File bam
    }

    command <<<
                        
        samtools addreplacerg -r ID:~{sample_id} -r LB:L1 -r SM:~{sample_id} -o ~{sample_id}_addRG.bam ~{bam}
                        
    >>>

    output {
        File rgbam = "${sample_id}_addRG.bam"
    }

    runtime {
        docker: "staphb/samtools:1.10"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }
}

task variant_calling {
    input {
        String sample_id
        File bam
        File ref

    }

    command <<<
        
        freebayes -f ~{ref} --haplotype-length 0 --min-alternate-count 3 --min-alternate-fraction 0.05 --min-mapping-quality 20 --min-base-quality 20 --min-coverage 10 --use-duplicate-reads --report-monomorphic --pooled-continuous ~{bam} > ~{sample_id}_variants.vcf
        
    >>>

    output {
        File vcf = "${sample_id}_variants.vcf"
    }

    runtime {
        docker: "wgspipeline/freebayes:v0.0.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 200 SSD"
    }
}

task sort_vcf {
    input {
        String sample_id
        File vcf

    }

    command <<<
        
        bgzip -c ~{vcf} > ~{sample_id}_variants.vcf.gz
        
        tabix -p vcf ~{sample_id}_variants.vcf.gz
        
        bcftools sort -O z ~{sample_id}_variants.vcf.gz > ~{sample_id}_sorted.vcf.gz
        
    >>>

    output {
        File sorted_vcf = "${sample_id}_sorted.vcf.gz"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }
}

task sample_VOCs {
    input {
        File vcf
        File bed
        String sample_id
    }

    command <<<
        
        tabix -p vcf ~{vcf}
        
        bcftools view --regions-file ~{bed} --output-type v --output-file ~{sample_id}_voc_mutations.vcf ~{vcf}
        
    >>>

    output {
         File sample_voc_vcf = "${sample_id}_voc_mutations.vcf"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task vcf2tsv {
    input {
        File vcf
        File bed
        String sample_id
    }

    command <<<
    
        bgzip -c ~{vcf} > ~{sample_id}_voc_mutations.vcf.gz
        
        tabix -p vcf ~{sample_id}_voc_mutations.vcf.gz
        
        bcftools query --regions-file ~{bed} --format '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%RO\t%AO]\n' ~{sample_id}_voc_mutations.vcf.gz > ~{sample_id}_voc_mutations.tsv
        
    >>>

    output {
         File sample_voc_tsv = "${sample_id}_voc_mutations.tsv"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task fill_NA {
    input {
        File tsv
        String sample_id
        File voc_bed
    }

    command <<<    
        
        # create key of unique locations
        cat ~{voc_bed} | cut -f 1,2 | tr "\t" "_" | sort | uniq > keys.txt
        
        # add headers to tsv and use key to fill in missing values
        echo -e "CHROM\tPOS\tREF\t~{sample_id}_ALT\t~{sample_id}_DP\t~{sample_id}_RO\t~{sample_id}_AO" | cat - ~{tsv} | sed 's/\t/_/' | sort -t $'\t' -k1,1 > ~{sample_id}_voc_mutations_temp1.tsv
        
        # get the filled columns we want
        join -t $'\t' -e NA -a 1 -1 1 -2 1 -o "1.1,2.3,2.4,2.6" keys.txt "~{sample_id}_voc_mutations_temp1.tsv" > ~{sample_id}_voc_fill_NA.tsv
            
    >>>

    output {
         File fill_NA_tsv = "${sample_id}_voc_fill_NA.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task allele_freq {
    input {
        File tsv
        String sample_id
    }

    command <<<    
        
        # separate the comma separated alleles into separate rows (might need to fix delimiters)
        awk '{split($2,a,","); split($4,b,","); for(i in a){print $1,a[i],$3,b[i]}}' ~{tsv} > ~{sample_id}_voc_mutations_temp2.tsv
        
        # use AO and DP fields to calculate ALT allele frequency, fix delimiters, change -nan allele frequencies to NA
        awk '$3~"^NA"||$4~"^NA"{$5="NA";print;next}{$5=$4/$3}1' ~{sample_id}_voc_mutations_temp2.tsv | sed 's/ /\t/g' | awk '$5 == "-nan" {$5="NA"} 1' OFS="\t" > ~{sample_id}_voc_allele_freq.tsv
    
    >>>

    output {
         File allele_freq_tsv = "${sample_id}_voc_allele_freq.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task reformat_tsv {
    input {
        File tsv
        String sample_id
    }

    command <<<    

        # combine the rows based on matching nucl location
        
        awk '{f2[$1]=f2[$1] sep[$1] $2; 
              f3[$1]=f3[$1] sep[$1] $3;
              f4[$1]=f4[$1] sep[$1] $4;
              f5[$1]=f5[$1] sep[$1] $5; 
              sep[$1]=";"}
         END {for(k in f2) print k,f2[k],f3[k],f4[k],f5[k]}' ~{tsv} > ~{sample_id}_voc_mutations_temp3.tsv
         
        # fix delimiters, add a column containing the sample ids
        sed 's/ /\t/g' ~{sample_id}_voc_mutations_temp3.tsv | awk 'NF=NF+1{$NF="~{sample_id}"}1' > ~{sample_id}_voc_mutations_temp4.tsv
        
        # fix the column headers, convert from space to tab delimited and then sort by col1
        echo -e "CHROMPOS ~{sample_id}_ALT ~{sample_id}_DP ~{sample_id}_AO ~{sample_id}_ALTfreq sample_id" | cat - ~{sample_id}_voc_mutations_temp4.tsv | sed 's/ /\t/g' | sort -t $'\t' -k 1,1 -V > ~{sample_id}_voc_reformat.tsv
        
    >>>

    output {
         File reformat_tsv_tsv = "${sample_id}_voc_reformat.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task summary_prep {
    input {
        File tsv
        String sample_id
        File voc_annotations
    }

    command <<<    
        
        # cut the columns we want for the results summary and make output file
        cut -f2,5 ~{tsv} > ~{sample_id}_voc_mutations_forsummary.tsv
        
        # cut the columns we want for the dashboard summary
        awk '{print $6 "\t" $2 "\t" $5}' ~{tsv} > ~{sample_id}_voc_mutations_temp5.tsv
        
        # add annotations to the dashboard summary, reorder the dashboard summary columns, fix the dashboard summary column headers and make output file
        paste ~{voc_annotations} ~{sample_id}_voc_mutations_temp5.tsv | awk '{print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6}' | awk 'BEGIN{FS=OFS="\t"; print "sample_id", "AA_change", "Nucl_change", "Lineages", "ALT", "ALTfreq"} NR>1{print $1, $2, $3, $4, $5, $6}' > ~{sample_id}_voc_mutations_fordash.tsv
    
        # cut the columns we want for the counts summary
        awk '{print $6 "\t" $2 "\t" $3 "\t" $4}' ~{tsv} > ~{sample_id}_voc_mutations_temp6.tsv
        
        # add annotations to the counts summary, reorder the dashboard summary columns, fix the dashboard summary column headers and make output file
        paste ~{voc_annotations} ~{sample_id}_voc_mutations_temp6.tsv | awk '{print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7}' | awk 'BEGIN{FS=OFS="\t"; print "sample_id", "AA_change", "Nucl_change", "Lineages", "ALT", "Total_count", "ALT_count"} NR>1{print $1, $2, $3, $4, $5, $6, $7}' > ~{sample_id}_voc_mutations_counts.tsv
   
   >>>

    output {
         File sample_voc_tsv_summary = "${sample_id}_voc_mutations_forsummary.tsv"
         File sample_voc_tsv_dash = "${sample_id}_voc_mutations_fordash.tsv"
         File sample_voc_tsv_counts = "${sample_id}_voc_mutations_counts.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task dashboard_tsv {
    input {
        Array[File] tsv
        Array[File] tsv_dash
        Array[File] tsv_counts
        File voc_annotations
    }

    command <<<
        
        # concatenate the tsvs and make the dashboard summary output
        awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=' ' tsv_dash} >> voc_mutations_dashboard.tsv
        
        # concatenate the tsvs and make the dashboard summary output
        awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=' ' tsv_counts} >> voc_mutations_counts.tsv
        
        # fix delimiters in annotations file
        sed 's/ /\t/g' ~{voc_annotations} > voc_annotations.tsv
        
        # concatentate tsvs for sequencing and bioinformatics team summary file and make output
        paste voc_annotations.tsv ~{sep=' ' tsv} > voc_mutations_summary_temp.tsv

    >>>

    output {
        File voc_summary_temp = "voc_mutations_summary_temp.tsv"
        File voc_dashboard = "voc_mutations_dashboard.tsv"
        File voc_counts = "voc_mutations_counts.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
    }
}

task summary_tsv {
    input {
        File tsv
    }

    command <<<
        
        # datamash to tranpose results summary
        datamash -H transpose < ~{tsv} > voc_mutations_summary.tsv

    >>>

    output {
        File voc_summary = "voc_mutations_summary.tsv"
    }

    runtime {
        docker: "rapatsky/debian"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
    }
}

task transfer_outputs {
    input {
        Array[File] variants
        Array[File] sorted_vcf
        Array[File] sample_voc_vcf
        Array[File] sample_voc_tsv
        Array[File] sample_voc_tsv_summary
        Array[File] sample_voc_tsv_dash
        Array[File] sample_voc_tsv_counts
        File voc_summary
        File voc_dashboard
        File voc_counts
        String out_dir
        
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<
        
        gsutil -m cp ~{sep=' ' variants} ~{outdirpath}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sorted_vcf} ~{outdirpath}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_voc_vcf} ~{outdirpath}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_voc_tsv} ~{outdirpath}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_voc_tsv_summary} ~{outdirpath}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_voc_tsv_dash} ~{outdirpath}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_voc_tsv_counts} ~{outdirpath}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{voc_summary} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{voc_dashboard} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{voc_counts} ~{outdirpath}/waste_water_variant_calling/
        
        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 50 SSD"
    }
}
