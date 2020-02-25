workflow test_multiqc {
    Array[File] star_logs
    Array[File] fc_logs
    Array[File] fastqc_zips

    call create_qc {
        input:
            star_logs = star_logs,
            fc_logs = fc_logs,
            r1_fastqc_zips = fastqc_zips
    }   

}

task create_qc {
    Array[File] star_logs
    Array[File] fc_logs
    Array[File] r1_fastqc_zips
    Array[File] dedup_metrics

    Int disk_size = 30

    command {
        multiqc .
    }

    output {
        File report = "multiqc_report.html"
    }
        
    runtime {
        docker: "docker.io/blawney/star_single_end_rnaseq:v0.0.2"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
