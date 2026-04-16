nextflow.enable.dsl=2

workflow {

    qc_out         = QC()
    preprocess_out = PREPROCESS(qc_out)
    cluster_out    = CLUSTER(preprocess_out)
    deg_out        = DEG(cluster_out)

    deg_csv  = deg_out[0]
    deg_h5ad = deg_out[1]

    annot_out   = ANNOTATE(deg_h5ad)
    pathway_out = PATHWAY(deg_csv)
    traj_out    = TRAJECTORY(annot_out)
}

process QC {

    publishDir "results/qc", mode: 'copy'

    output:
    path "qc.h5ad"

    script:
    """
    python ${projectDir}/scripts/step01_qc.py --output qc.h5ad
    """
}

process PREPROCESS {

    publishDir "results/preprocess", mode: 'copy'

    input:
    path qc_file

    output:
    path "preprocessed.h5ad"

    script:
    """
    python ${projectDir}/scripts/step02_preprocessing.py \
        --input $qc_file \
        --output preprocessed.h5ad
    """
}

process CLUSTER {

    publishDir "results/cluster", mode: 'copy'

    input:
    path input_file

    output:
    path "clustered.h5ad"

    script:
    """
    python ${projectDir}/scripts/step03_clustering.py \
        --input $input_file \
        --output clustered.h5ad
    """
}

process DEG {

    publishDir "results/deg", mode: 'copy'

    input:
    path input_file

    output:
    path "deg.csv"
    path "deg.h5ad"

    script:
    """
    python ${projectDir}/scripts/step04_deg_analysis.py \
        --input $input_file \
        --out_prefix deg
    """
}


process ANNOTATE {

    publishDir "results", mode: 'copy'

    input:
    path deg_file

    output:
    path "annotated.h5ad"

    script:
    """
    python ${projectDir}/scripts/step05_celltype_annotation.py \
        --input $deg_file \
        --output annotated.h5ad
    """
}


process PATHWAY {

    publishDir "results/pathway", mode: 'copy'

    input:
    path deg_file

    output:
    path "gsea_results.csv"

    script:
    """
    python ${projectDir}/scripts/step06_pathway_enrichment.py \
        --deg $deg_file \
        --out_prefix gsea
    """
}

process TRAJECTORY {

    publishDir "results/trajectory", mode: 'copy'

    input:
    path input_file

    output:
    path "trajectory.h5ad"

    script:
    """
    python ${projectDir}/scripts/step07_trajectory_analysis.py \
        --input $input_file \
        --root_cluster 0 \
        --output trajectory.h5ad
    """
}
