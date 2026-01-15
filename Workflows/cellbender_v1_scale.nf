nextflow.enable.dsl = 2

params.samples = ""
params.epochs = ""
params.low_count_threshold = ""
params.z_dim = ""
params.z_layers = ""
params.training_fraction = ""
params.empty_drop_training_fraction = ""
params.fpr = ""
params.projected_ambient_count_threshold = ""
params.learning_rate = ""
params.posterior_batch_size = ""
params.debug = ""

channel
    .fromPath(params.samples)
    .splitCsv(header: ['sampleName', 's3path'], sep:',', skip: 1)
    .set{ samples_ch }

process cellbender {
    container '592777438441.dkr.ecr.us-east-1.amazonaws.com/scrna-seq_v1:CellBender_v1_default'
    
    // Request GPU for this process
    accelerator 1, type: 'nvidia-tesla-a10g'

    // Define memory requirement
    memory '32GB'

    // Request multiple CPUs for this process
    cpus 8

    publishDir "/mnt/workflow/pubdir"

    input:
    tuple val(sampleName), path(s3path)

    output:
    path "${sampleName}_output*", emit: outputs

    script:
    """
    cellbender remove-background \
      --cuda \
      --epochs ${params.epochs} \
      --low-count-threshold ${params.low_count_threshold} \
      --z-dim ${params.z_dim} \
      --z-layers ${params.z_layers} \
      --training-fraction ${params.training_fraction} \
      --empty-drop-training-fraction ${params.empty_drop_training_fraction} \
      --fpr ${params.fpr} \
      --projected-ambient-count-threshold ${params.projected_ambient_count_threshold} \
      --learning-rate ${params.learning_rate} \
      --posterior-batch-size ${params.posterior_batch_size} \
      ${params.debug ? '--debug' : ''} \
      --input ${s3path} \
      --output ${sampleName}_output.h5
    """
}

workflow {
    samples_ch
        .map { row -> tuple(row.sampleName, file(row.s3path)) }
        .set { input_files }
    cellbender(input_files)
}
