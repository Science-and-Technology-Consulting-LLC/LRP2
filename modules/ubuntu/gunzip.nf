process GUNZIP {
    tag "${archive}"
    label 'process_single'

    container "ubuntu:22.04"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${gunzip}"), emit: gunzip
    tuple val("${task.process}"), val('gunzip'), eval('gunzip --version 2>&1 | head -1 | sed "s/^.*(gzip) //; s/ Copyright.*//"'), topic: versions, emit: versions_gunzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def extension = (archive.toString() - '.gz').tokenize('.')[-1]
    def name = archive.toString() - '.gz' - ".${extension}"
    def prefix = task.ext.prefix ?: name
    gunzip = prefix + ".${extension}"
    """
    # Not calling gunzip itself because it creates files
    # with the original group ownership rather than the
    # default one for that user / the work directory
    gzip \\
        -cd \\
        ${args} \\
        ${archive} \\
        > ${gunzip}
    """

    stub:
    def extension = (archive.toString() - '.gz').tokenize('.')[-1]
    def name = archive.toString() - '.gz' - ".${extension}"
    def prefix = task.ext.prefix ?: name
    gunzip = prefix + ".${extension}"
    """
    touch ${gunzip}
    """
}
