task create_contrast_independent_figures {

    File sample_annotations
    File raw_count_matrix
    String pca_filename
    String hc_tree_filename

    Int disk_size = 10

    command <<<
        source activate r36
        export OWD=$(pwd)
        cd /opt/software/
        Rscript contrast_independent_figures.R \
            ${raw_count_matrix} \
            ${sample_annotations} \
            $OWD/${pca_filename} \
            $OWD/${hc_tree_filename}
    >>>

    output {
        File pca = "${pca_filename}"
        File hctree = "${hc_tree_filename}"
    }

    runtime {
        docker: "docker.io/blawney/star_single_end_rnaseq:v0.0.2"
        cpu: 2
        memory: "3 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
