nextflow.enable.dsl = 2

params.samples = ""
params.metadata_csv = ""
params.percent_mito = ""
params.min_features = ""
params.max_features = ""
params.memory = ""
params.cpu = ""

process download_metadata {
    output:
    path "metadata_downloaded.csv" into metadata_ch

    script:
    """
    aws s3 cp ${params.metadata_csv} metadata_downloaded.csv
    """
}

process run_analysis {
    container '592777438441.dkr.ecr.us-east-1.amazonaws.com/scrna-seq_v1:scDblFinder_v1_default'

    memory "${params.memory}"
    cpus "${params.cpu}"

    publishDir "/mnt/workflow/pubdir"

    input:
    tuple val(sample_name), path(h5_file), path(metadata_file)

    output:
    path("${sample_name}.rds")
    path("${sample_name}.h5ad")

    script:
    """
    export R_LIBS_USER='/usr/local/lib/R/host-site-library'
    Rscript /opt/scripts/scDblFinder_v1_scale.R ${h5_file} ${metadata_file} ${sample_name} ${sample_name}.rds ${sample_name}.h5ad ${params.percent_mito} ${params.min_features} ${params.max_features}
    """
}

workflow {
    Channel
        .fromPath(params.metadata_csv)
        .set { metadata_ch }

    samples_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> tuple(row.sampleName, file(row.s3path)) }

    samples_ch
        .combine(metadata_ch)
        .set { samples_with_metadata_ch }

    run_analysis(samples_with_metadata_ch)
}
