version 1.0

workflow WasteWaterVariantCalling {

    input {
        Array[File] sorted_bam
        File covid_genome
        File spike_bed
        Array[String] sample_id
        String out_dir
    }

    scatter (id_bam in zip(sample_id, sorted_bam)) {
        call variant_calling {
            input:
                bam = id_bam.right,
                ref = covid_genome,
                sample_id = id_bam.left
        }
        call sort_vcf {
            input:
                vcf = variant_calling.vcf,
                sample_id = id_bam.left
        }
       
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
            out_dir = out_dir
    }
    
    output {
        Array[File] variants = variant_calling.vcf
        Array[File] sorted_vcf = sort_vcf.sorted_vcf
        File merged_vcf = merge_vcfs.merged_vcf
        File spike_mutations = get_spike.spike_mutations
        String transfer_date = transfer_outputs.transfer_date
    }
}

task variant_calling {
    input {
        String sample_id
        File bam
        File ref

    }

    command <<<
        
        bcftools mpileup --threads 4 --count-orphans --ignore-RG --max-depth 1000 --fasta-ref ~{ref} --min-MQ 20 --output-type u ~{bam} | \
        bcftools call --threads 4 --keep-alts --keep-masked-ref --multiallelic-caller --output-type u | \
        bcftools filter --output-type v --soft-filter LowQual --include 'QUAL>=20 && DP>=10' > ~{sample_id}_variants.vcf
        
    >>>

    output {
        File vcf = "${sample_id}_variants.vcf"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
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
        
        bcftools sort -O z ~{vcf} > ~{sample_id}_sorted.vcf.gz
        
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
        File merged_vcf
        File spike_mutations
        String out_dir
        
    }
    
    String outdir = sub(out_dir, "/$", "")

    command <<<
        
        gsutil -m cp ~{sep=' ' variants} ~{outdir}/waste_water_variant_calling/vcfs/
        gsutil -m cp ~{sep=' ' sorted_vcf} ~{outdir}/waste_water_variant_calling/vcfs/
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