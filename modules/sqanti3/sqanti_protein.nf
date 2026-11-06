process SQANTI_PROTEIN {
    tag "$meta.id"
    label 'process_long'

    container "${params.image_registry}/sqanti3_protein:${params.image_tag}"

    input:
    tuple val(meta), path(cds_gtf)
    path reference_gtf

    output:
    tuple val(meta), path("*.predicted_proteome.best_ORF_SQANTI_classification.tsv"), emit: protein_classification
    tuple val(meta), path("*_S3_PREDICTED_PROTEOME_M3_SQANTI_PROTEIN_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    exec > >(tee ${prefix}_S3_PREDICTED_PROTEOME_M3_SQANTI_PROTEIN_log.txt) 2>&1

    source /conda/miniconda3/etc/profile.d/conda.sh
    conda activate sqanti3
    export SQANTI_PATH=\$(dirname \$(which sqanti3_qc.py))

    # Add src/utilities and utilities to PYTHONPATH for cupcake and other imports
    export PYTHONPATH=\${SQANTI_PATH}/src/utilities:\${SQANTI_PATH}/utilities:\${SQANTI_PATH}:\${PYTHONPATH:-}

    # the script is in the image, already patched for the gtfToGenePred binary that v5.5.4 ships
    python /opt/sqanti3_protein/04_sqanti3_protein.py \\
        $cds_gtf \\
        $reference_gtf \\
        -d . \\
        -p $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.5.4
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.predicted_proteome.best_ORF_SQANTI_classification.tsv
    touch ${prefix}_S3_PREDICTED_PROTEOME_M3_SQANTI_PROTEIN_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
