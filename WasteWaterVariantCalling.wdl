version 1.0

workflow WasteWaterVariantCalling {

    input {
        Array[File] sorted_bam
        File covid_genome
        File spike_bed
        File spike_annotations
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
        call sample_spike {
            input:
                vcf = sort_vcf.sorted_vcf,
                bed = spike_bed,
                sample_id = id_bam.left
        }
        call vcf2tsv {
            input:
                vcf = sample_spike.sample_spike_vcf,
                sample_id = id_bam.left,
                bed = spike_bed
        }
        call reformat_tsv {
            input:
                tsv = vcf2tsv.sample_spike_tsv,
                sample_id = id_bam.left,
                spike_bed = spike_bed,
                spike_annotations = spike_annotations
        }
    }
    call merge_tsv {
        input:
            tsv = reformat_tsv.sample_spike_tsv_summary,
            tsv_dash = reformat_tsv.sample_spike_tsv_dash,
            spike_annotations = spike_annotations
    }
    
    call merge_vcfs {
        input:
            vcfs = sort_vcf.sorted_vcf       
    }
    call get_spike {
       input:
           vcf = merge_vcfs.merged_vcf,
           bed = spike_bed
    }
    call transfer_outputs {
        input:
            variants = variant_calling.vcf,
            sorted_vcf = sort_vcf.sorted_vcf,
            merged_vcf = merge_vcfs.merged_vcf,
            spike_mutations = get_spike.spike_mutations,
            sample_spike_vcf = sample_spike.sample_spike_vcf,
            sample_spike_tsv = vcf2tsv.sample_spike_tsv,
            sample_spike_tsv_summary = reformat_tsv.sample_spike_tsv_summary,
            sample_spike_tsv_dash = reformat_tsv.sample_spike_tsv_dash,
            spike_summary = merge_tsv.spike_summary,
            spike_dashboard = merge_tsv.spike_dashboard,
            out_dir = out_dir
    }
    
    output {
        Array[File] addrg_bam = add_RG.rgbam
        Array[File] variants = variant_calling.vcf
        Array[File] sorted_vcf = sort_vcf.sorted_vcf
        Array[File] sample_spike_vcf = sample_spike.sample_spike_vcf
        Array[File] sample_spike_tsv = vcf2tsv.sample_spike_tsv
        Array[File] sample_spike_tsv_summary = reformat_tsv.sample_spike_tsv_summary
        Array[File] sample_spike_tsv_dash = reformat_tsv.sample_spike_tsv_dash
        File merged_vcf = merge_vcfs.merged_vcf
        File spike_summary = merge_tsv.spike_summary
        File spike_dashboard = merge_tsv.spike_dashboard
        File spike_mutations = get_spike.spike_mutations
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
        
        freebayes -f ~{ref} --haplotype-length 0 --min-alternate-count 3 --min-alternate-fraction 0.05 --min-mapping-quality 20 --min-base-quality 20 --min-coverage 10 --use-duplicate-reads --report-monomorphic --pooled-continuous -n 1 ~{bam} > ~{sample_id}_variants.vcf
        
    >>>

    output {
        File vcf = "${sample_id}_variants.vcf"
    }

    runtime {
        docker: "wgspipeline/freebayes:v0.0.1"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
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

task sample_spike {
    input {
        File vcf
        File bed
        String sample_id
    }

    command <<<
        
        tabix -p vcf ~{vcf}
        
        bcftools view --regions-file ~{bed} --output-type v --output-file ~{sample_id}_spike_mutations.vcf ~{vcf}
        
    >>>

    output {
         File sample_spike_vcf = "${sample_id}_spike_mutations.vcf"
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
    
        bgzip -c ~{vcf} > ~{sample_id}_spike_mutations.vcf.gz
        
        tabix -p vcf ~{sample_id}_spike_mutations.vcf.gz
        
        bcftools query --regions-file ~{bed} --format '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%RO\t%AO]\n' ~{sample_id}_spike_mutations.vcf.gz > ~{sample_id}_spike_mutations.tsv
        
    >>>

    output {
         File sample_spike_tsv = "${sample_id}_spike_mutations.tsv"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task reformat_tsv {
    input {
        File tsv
        String sample_id
        File spike_bed
        File spike_annotations
    }

    command <<<
    
        echo -e "CHROM\tPOS\tREF\t~{sample_id}_ALT\t~{sample_id}_DP\t~{sample_id}_RO\t~{sample_id}_AO" | cat - ~{tsv} > ~{sample_id}_spike_mutations_headers.tsv    
        
        cat ~{spike_bed}  | cut -f 1,2 | tr "\t" "_" | sort | uniq > keys.txt
        
        sed 's/\t/_/' "~{sample_id}_spike_mutations_headers.tsv" | sort -t $'\t' -k1,1 > "~{sample_id}_spike_mutations_headers_temp.tsv"
        
        join -t $'\t' -e NA -a 1 -1 1 -2 1 -o "2.3,2.6,2.4" keys.txt "~{sample_id}_spike_mutations_headers_temp.tsv" > "~{sample_id}_spike_mutations_temp2.tsv"
        
        awk '{$4 = $2 / $3}1' ~{sample_id}_spike_mutations_temp2.tsv > ~{sample_id}_spike_mutations_temp3.tsv
        
        awk 'NF=NF+1{$NF="~{sample_id}"}1'  ~{sample_id}_spike_mutations_temp3.tsv >  ~{sample_id}_spike_mutations_temp4.tsv
        
        echo -e "~{sample_id}_ALT ~{sample_id}_AO ~{sample_id}_DP ~{sample_id}_ALTfreq sample_id" | cat - ~{sample_id}_spike_mutations_temp4.tsv > ~{sample_id}_spike_mutations_temp5.tsv
        
        sed 's/ /\t/g' ~{sample_id}_spike_mutations_temp5.tsv > ~{sample_id}_spike_mutations_temp6.tsv
        
        awk '$4 == "-nan" {$4="NA"} 1' OFS="\t" ~{sample_id}_spike_mutations_temp6.tsv > ~{sample_id}_spike_mutations_temp7.tsv
        
        cut -f1,4 ~{sample_id}_spike_mutations_temp7.tsv > ~{sample_id}_spike_mutations_forsummary.tsv
        
        awk '{print $5 "\t" $1 "\t" $4}' ~{sample_id}_spike_mutations_temp7.tsv > ~{sample_id}_spike_mutations_temp8.tsv
        
        paste ~{spike_annotations} ~{sample_id}_spike_mutations_temp8.tsv > ~{sample_id}_spike_mutations_temp9.tsv
        
        awk '{print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6}' ~{sample_id}_spike_mutations_temp9.tsv > ~{sample_id}_spike_mutations_temp10.tsv
        
        awk 'BEGIN{FS=OFS="\t"; print "sample_id", "AA_change", "Nucl_change", "Lineages", "ALT", "ALTfreq"} NR>1{print $1, $2, $3, $4, $5, $6}' ~{sample_id}_spike_mutations_temp10.tsv > ~{sample_id}_spike_mutations_fordash.tsv
        
        ### change counts to alt frequency - done
        ### For epis, need to add a column with sample name all the way down - done
        ### transpose tables for bioinf
        ### cat the long way for epis - done
    
    >>>

    output {
         File sample_spike_tsv_summary = "${sample_id}_spike_mutations_forsummary.tsv"
         File sample_spike_tsv_dash = "${sample_id}_spike_mutations_fordash.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task merge_tsv {
    input {
        Array[File] tsv
        Array[File] tsv_dash
        File spike_annotations
    }

    command <<<
        
        awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=' ' tsv_dash} >> spike_mutations_dashboard.tsv
        ### echo -e "sample_id\tAA_change\tNucl_change\tLineages\tALT\tALTfreq" | cat - spike_mutations_dashboard_temp.tsv > spike_mutations_dashboard.tsv
        
        sed 's/ /\t/g' ~{spike_annotations} > spike_annotations.tsv
        paste spike_annotations.tsv ~{sep=' ' tsv} > spike_mutations_summary.tsv
        
        ### datamash -H transpose < spike_mutations_summary_temp.tsv > spike_mutations_summary.tsv
        ### might need to make datamash it's own task

    >>>

    output {
        File spike_summary = "spike_mutations_summary.tsv"
        File spike_dashboard = "spike_mutations_dashboard.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task merge_vcfs {
    input {
        Array[File] vcfs
    }

    command <<<
        
        # copy files to CWD (required for tabix indexing)
        INFILES=~{write_lines(vcfs)}
        ln -s $(cat $INFILES) .
        for FN in $(cat $INFILES); do basename $FN; done > vcf_filenames.txt

        # tabix index input vcfs (must be gzipped)
        for FN in $(cat vcf_filenames.txt); do
            tabix -p vcf $FN
        done
        
        bcftools merge --missing-to-ref --force-samples --merge snps --output-type z --output merged_variants.vcf.gz --file-list vcf_filenames.txt

    >>>

    output {
        File merged_vcf = "merged_variants.vcf.gz"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task get_spike {
    input {
        File vcf
        File bed
    }

    command <<<
        
        tabix -p vcf ~{vcf}
        
        bcftools view --regions-file ~{bed} --output-type v --output-file spike_mutations.vcf ~{vcf}
        
    >>>

    output {
         File spike_mutations = "spike_mutations.vcf"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task transfer_outputs {
    input {
        Array[File] variants
        Array[File] sorted_vcf
        Array[File] sample_spike_vcf
        Array[File] sample_spike_tsv
        Array[File] sample_spike_tsv_summary
        Array[File] sample_spike_tsv_dash
        File merged_vcf
        File spike_summary
        File spike_dashboard
        File spike_mutations
        String out_dir
        
    }
    
    String outdir = sub(out_dir, "/$", "")

    command <<<
        
        gsutil -m cp ~{sep=' ' variants} ~{outdir}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sorted_vcf} ~{outdir}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_spike_vcf} ~{outdir}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_spike_tsv} ~{outdir}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_spike_tsv_summary} ~{outdir}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sample_spike_tsv_dash} ~{outdir}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{spike_summary} ~{outdir}/waste_water_variant_calling/
        gsutil -m cp ~{spike_dashboard} ~{outdir}/waste_water_variant_calling/
        gsutil -m cp ~{merged_vcf} ~{outdir}/waste_water_variant_calling/
        gsutil -m cp ~{spike_mutations} ~{outdir}/waste_water_variant_calling/
        
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
        disks: "local-disk 10 SSD"
    }
}