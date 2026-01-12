nextflow.enable.dsl = 2

// Define parameters
params.samples_csv = ""
params.batch_key = ""
params.base_output_name = ""
params.HVG = ""
params.continuous_covariate_keys = ""
params.categorical_covariate_keys = ""
params.gene_likelihood = ""
params.n_latent = ""
params.n_layers = ""
params.n_hidden = ""
params.dispersion = ""
params.check_val_every_n_epoch = ""
params.max_epochs = ""
params.early_stopping = ""
params.early_stopping_patience = ""
params.early_stopping_monitor = ""
params.lr = ""
params.optimizer = ""
params.dropout_rate = ""
params.cpu_model = ""
params.memory_model = ""
params.cpu_norm = ""
params.memory_norm = ""

// Create a channel for reading CSV file
Channel
    .fromPath(params.samples_csv)
    .splitCsv(header: true, sep: ',', strip: true)
    .map { row -> tuple(row.sampleName, row.s3path) }
    .set { samples_ch }

// Process 1: Download Files from S3
process DownloadFromS3 {
    container '592777438441.dkr.ecr.us-east-1.amazonaws.com/scrna-seq_v1:scVI-tools_v1_default'
 
    input:
    tuple val(sampleName), val(s3path)

    output:
    path "${sampleName}.h5ad"

    script:
    """
    aws s3 cp ${s3path} ${sampleName}.h5ad
    """
}

// Process 2: Merge Samples
process MergeFiles {
    container '592777438441.dkr.ecr.us-east-1.amazonaws.com/scrna-seq_v1:scVI-tools_v1_default'

    memory "${params.memory_model}"
    cpus "${params.cpu_model}"

    input:
    path h5ad_files

    output:
    path "merged_adata.h5ad"

    script:
    """
    python /app/scVI_merge_v1_scale.py ${h5ad_files.join(' ')}
    """
}

// Process 3: scVI Modeling
process SCVI_ModelTraining {
    container '592777438441.dkr.ecr.us-east-1.amazonaws.com/scrna-seq_v1:scVI-tools_v1_default'

    accelerator 1, type: 'nvidia-tesla-a10g'

    memory "${params.memory_model}"
    cpus "${params.cpu_model}"

    publishDir "/mnt/workflow/pubdir"

    input:
    path "merged_adata.h5ad"

    output:
    path "${params.base_output_name}_scVI_model.pt"
    path "${params.base_output_name}_partial_scVI_adata.h5ad"

    script:
    def HVG = params.HVG ? "--HVG ${params.HVG}" : "--HVG 2000"
    def batchKey = params.batch_key?.trim() ? "--batch_key '${params.batch_key}'" : ""
    def continuousCovKeys = params.continuous_covariate_keys ? "--continuous_covariate_keys '${params.continuous_covariate_keys}'" : ""
    def categoricalCovKeys = params.categorical_covariate_keys ? "--categorical_covariate_keys '${params.categorical_covariate_keys}'" : ""
    def geneLikelihood = params.gene_likelihood ? "--gene_likelihood ${params.gene_likelihood}" : "--gene_likelihood 'nb'"
    def nLatent = params.n_latent ? "--n_latent ${params.n_latent}" : "--n_latent 10"
    def nLayers = params.n_layers ? "--n_layers ${params.n_layers}" : "--n_layers 2"
    def nHidden = params.n_hidden ? "--n_hidden ${params.n_hidden}" : "--n_hidden 128"
    def dispersion = params.dispersion ? "--dispersion ${params.dispersion}" : "--dispersion 'gene'"
    def checkValEveryNEpoch = params.check_val_every_n_epoch ? "--check_val_every_n_epoch ${params.check_val_every_n_epoch}" : "--check_val_every_n_epoch 1"
    def maxEpochs = params.max_epochs ? "--max_epochs ${params.max_epochs}" : "--max_epochs 400"
    def earlyStopping = params.early_stopping ? "--early_stopping ${params.early_stopping}" : "--early_stopping True"
    def earlyStoppingPatience = params.early_stopping_patience ? "--early_stopping_patience ${params.early_stopping_patience}" : "--early_stopping_patience 50"
    def earlyStoppingMonitor = params.early_stopping_monitor ? "--early_stopping_monitor ${params.early_stopping_monitor}" : "--early_stopping_monitor 'elbo_validation'"
    def lr = params.lr ? "--lr ${params.lr}" : "--lr 0.001"
    def optimizer = params.optimizer ? "--optimizer ${params.optimizer}" : "--optimizer 'AdamW'"    
    def dropout_rate = params.dropout_rate ? "--dropout_rate ${params.dropout_rate}" : "--dropout_rate 0.1"

    """
    python /app/scVI_model_v1_scale.py --input_adata merged_adata.h5ad $HVG $batchKey $continuousCovKeys $categoricalCovKeys $geneLikelihood $nLatent $nLayers $nHidden $dispersion $checkValEveryNEpoch $maxEpochs $earlyStopping $earlyStoppingPatience $earlyStoppingMonitor $lr $optimizer $dropout_rate --output_model_path ${params.base_output_name}_scVI_model.pt --output_adata_path ${params.base_output_name}_partial_scVI_adata.h5ad
    """
}

// Process 4: Normalized Expression
process SCVI_NormalizedExpression {
    container '592777438441.dkr.ecr.us-east-1.amazonaws.com/scrna-seq_v1:scVI-tools_v1_default'

    memory "${params.memory_norm}"
    cpus "${params.cpu_norm}"

    publishDir "/mnt/workflow/pubdir"

    input:
    path "${params.base_output_name}_scVI_model.pt"
    path "${params.base_output_name}_partial_scVI_adata.h5ad"

    output:
    path "${params.base_output_name}_scVI_adata.h5ad"

    script:
    """
    python /app/scVI_normalized-expression_v1_scale.py --model_path '${params.base_output_name}_scVI_model.pt' --adata_path '${params.base_output_name}_partial_scVI_adata.h5ad' --output_adata_path '${params.base_output_name}_scVI_adata.h5ad'
    """
}

// Define the workflow
workflow {
    download = DownloadFromS3(samples_ch)
    download_files = download.collect()
    merged_adata = MergeFiles(download_files)
    model_training_output = SCVI_ModelTraining(merged_adata)
    SCVI_NormalizedExpression(model_training_output)
}
