workflow SingleEndRnaSeqAndDgeWorkflow {
    Array[File] r1_files
    File sample_annotations
    Array[String] base_conditions
    Array[String] experimental_conditions
    String genome
    File star_index_path
    File gtf
    File bed_annotations
    String output_zip_name
    String git_repo_url
    String git_commit_hash


    scatter(fastq in r1_files){

        call assert_valid_fastq {
            input:
                r1_file = fastq        }
    }

    call assert_valid_annotations{
        input:
            r1_files = r1_files,
            sample_annotations = sample_annotations,
            base_conditions = base_conditions,
            experimental_conditions = experimental_conditions
    }
}

task assert_valid_fastq {

    File r1_file
    Int disk_size = 100

    command <<<
        python3 /opt/software/precheck/check_fastq.py -r1 ${r1_file}
    >>>

    runtime {
        docker: "docker.io/blawney/star_single_end_rnaseq:v0.0.1"
        cpu: 2
        memory: "30 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task assert_valid_annotations {

    Array[String] r1_files
    File sample_annotations
    Array[String] base_conditions
    Array[String] experimental_conditions

    Int disk_size = 10

    command <<<
        python3 /opt/software/precheck/perform_precheck.py \
            -a ${sample_annotations} \
            -r1 ${sep=" " r1_files} \
            -x ${sep=" " base_conditions} \
            -y ${sep=" " experimental_conditions}
    >>>

    runtime {
        docker: "docker.io/blawney/star_single_end_rnaseq:v0.0.1"
        cpu: 2
        memory: "3 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
