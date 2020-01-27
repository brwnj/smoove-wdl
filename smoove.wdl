workflow smoove {

    File manifest
    Array[Array[File]] sample_data = read_tsv(manifest)
    File fasta
    File fasta_index
    File bed
    File gff
    File ? ped
    File ? known_sites
    String project_id
    String exclude_chroms
    Boolean sensitive

    Int small_disk
    Int medium_disk
    Int large_disk

    scatter (sample in sample_data) {
        call smoove_call {
            input:
                sample_id = sample[0],
                alignments = sample[1],
                alignments_index = sample[2],
                fasta = fasta,
                fasta_index = fasta_index,
                bed = bed,
                exclude_chroms = exclude_chroms,
                sensitive = sensitive,
                disk_size = small_disk
        }
    }
    call smoove_merge {
        input:
            project_id = project_id,
            fasta = fasta,
            fasta_index = fasta_index,
            vcfs = smoove_call.called_vcf,
            vcf_indexes = smoove_call.called_vcf_index,
            disk_size = small_disk
    }
    scatter (sample in sample_data) {
        call smoove_genotype {
            input:
                sample_id = sample[0],
                alignments = sample[1],
                alignments_index = sample[2],
                fasta = fasta,
                fasta_index = fasta_index,
                vcf = smoove_merge.merged_vcf,
                disk_size = small_disk
        }
    }
    call smoove_square {
        input:
            project_id = project_id,
            vcfs = smoove_genotype.genotyped_vcf,
            vcf_indexes = smoove_genotype.genotyped_vcf_index,
            gff = gff,
            disk_size = small_disk
    }
    meta {
        author: "Joe Brown"
        email: "brwnjm@gmail.com"
        description: "Run @brentp's smoove to call SVs"
    }
}


task smoove_call {

    String sample_id
    File alignments
    File alignments_index
    File fasta
    File fasta_index
    File bed
    String exclude_chroms
    Boolean sensitive

    Int disk_size

    command {
        smoove call --genotype --name ${sample_id} --processes 1 \
            --fasta ${fasta} --exclude ${bed} --excludechroms "${exclude_chroms}" \
            ${alignments} 2> ${sample_id}-smoove-call.log
        bcftools stats ${sample_id}-smoove.genotyped.vcf.gz > ${sample_id}-stats.txt
    }
    runtime {
        memory: "8 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        docker: "brentp/smoove:v0.2.5"
    }

    output {
        File called_vcf = "${sample_id}-smoove.genotyped.vcf.gz"
        File called_vcf_index = "${sample_id}-smoove.genotyped.vcf.gz.csi"
        File variant_counts = "${sample_id}-stats.txt"
        File sequence_counts = "${sample_id}-smoove-call.log"
    }
}


task smoove_merge {

    String project_id
    File fasta
    File fasta_index
    Array[File] vcfs
    Array[File] vcf_indexes

    Int disk_size

    command {
        smoove merge --name ${project_id} --fasta ${fasta} ${sep=' ' vcfs}
    }

    runtime {
        memory: "16 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        docker: "brentp/smoove:v0.2.5"
    }

    output {
        File merged_vcf = "${project_id}.sites.vcf.gz"
    }
}


task smoove_genotype {

    String sample_id
    File alignments
    File alignments_index
    File fasta
    File fasta_index
    File vcf

    Int disk_size

    command {
        wget -q https://raw.githubusercontent.com/samtools/samtools/develop/misc/seq_cache_populate.pl
        perl seq_cache_populate.pl -root `pwd`/cache ${fasta} 1> /dev/null 2> err || (cat err; exit 2)
        export REF_PATH=`pwd`/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
        export REF_CACHE=xx

        samtools quickcheck -v ${alignments}
        smoove genotype --duphold --processes 1 --removepr --outdir `pwd`/ --name ${sample_id} --fasta ${fasta} --vcf ${vcf} ${alignments}
    }

    runtime {
        memory: "16 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        docker: "brentp/smoove:v0.2.5"
    }

    output {
        File genotyped_vcf = "${sample_id}-smoove.genotyped.vcf.gz"
        File genotyped_vcf_index = "${sample_id}-smoove.genotyped.vcf.gz.csi"
    }
}


task smoove_square {

    String project_id
    Array[File] vcfs
    Array[File] vcf_indexes
    File gff

    Int disk_size

    command {
        smoove paste --outdir `pwd`/ --name ${project_id} ${sep=' ' vcfs}
        smoove annotate --gff ${gff} ${project_id}.smoove.square.vcf.gz | bgzip --threads 1 -c > ${project_id}.smoove.square.anno.vcf.gz
        bcftools index ${project_id}.smoove.square.anno.vcf.gz
    }

    runtime {
        memory: "16 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        docker: "brentp/smoove:v0.2.5"
    }

    output {
        File annotated_vcf = "${project_id}.smoove.square.anno.vcf.gz"
        File annotated_vcf_index = "${project_id}.smoove.square.anno.vcf.gz.csi"
    }
}
