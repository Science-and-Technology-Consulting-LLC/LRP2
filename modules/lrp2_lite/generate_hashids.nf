process GENERATE_HASHIDS {
    tag "$meta.id"
    label 'process_low'

    container "${params.image_registry}/lrp2_lite:${params.image_tag}"

    input:
    tuple val(meta), path(corrected_gtf), path(classification)
    path reference_gtf

    output:
    tuple val(meta), path("*.transcriptome.hashids_mapping.txt"), emit: hashids_mapping
    tuple val(meta), path("*.transcriptome.psl"), emit: corrected_psl
    tuple val(meta), path("*_S2_TRANSCRIPTOME_M2_GENERATE_HASHIDS_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    exec > >(tee ${prefix}_S2_TRANSCRIPTOME_M2_GENERATE_HASHIDS_log.txt) 2>&1

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript /opt/lrp2_lite/00_generate_hashids.R \\
        --basename ${prefix} \\
        --sample_gtf ${corrected_gtf} \\
        --classification ${classification} \\
        --reference_gtf ${reference_gtf} \\
        --hashlib_script /opt/lrp2_lite/hashlib_id_generator.py \\
        --output_dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcriptome.hashids_mapping.txt
    touch ${prefix}.transcriptome.psl
    touch ${prefix}_S2_TRANSCRIPTOME_M2_GENERATE_HASHIDS_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
        python: 3.11.0
    END_VERSIONS
    """
}
