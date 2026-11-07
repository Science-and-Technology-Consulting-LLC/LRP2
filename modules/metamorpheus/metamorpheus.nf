process METAMORPHEUS {
    tag "$meta.id"
    label 'process_high_memory'

    container "${params.image_registry}/metamorpheus:${params.image_tag}"

    input:
    tuple val(meta), path(mzml_files)
    path database_fasta
    path config_toml

    output:
    tuple val(meta), path("*.psmtsv"), emit: psm_table
    tuple val(meta), path("results/**"), emit: results
    tuple val(meta), path("*_S5_PROTEOMICS_M3_METAMORPHEUS_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Rename mzML files to simple names based on meta.id
    def mzml_array = mzml_files instanceof List ? mzml_files : [mzml_files]
    def symlink_commands = []
    def renamed_files = []
    mzml_array.eachWithIndex { file, idx ->
        def simple_name = mzml_array.size() > 1 ? "${prefix}_${idx + 1}.mzML" : "${prefix}.mzML"
        symlink_commands.add("ln -s ${file} ${simple_name}")
        renamed_files.add(simple_name)
    }
    def mzml_list = renamed_files.join(' ')

    """
    exec > >(tee ${prefix}_S5_PROTEOMICS_M3_METAMORPHEUS_log.txt) 2>&1

    # Rename mzML files to simple names for cleaner MetaMorpheus output
    ${symlink_commands.join('\n    ')}

    dotnet /MetaMorpheus/CMD.dll \\
        -t $config_toml \\
        -s $mzml_list \\
        -d $database_fasta \\
        -o results/ \\
        $args

    find results -name "AllPeptides*.psmtsv" -exec cp {} ${prefix}.psmtsv \\;
    if [ ! -f "${prefix}.psmtsv" ]; then
        find results -name "*.psmtsv" -type f -print -quit -exec cp {} ${prefix}.psmtsv \\;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamorpheus: \$(dotnet /MetaMorpheus/CMD.dll --version 2>&1 | grep -oP 'MetaMorpheus version \\K[0-9.]+' || echo "1.0.5")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results
    touch ${prefix}.psmtsv
    touch results/AllPeptides.psmtsv
    touch ${prefix}_S5_PROTEOMICS_M3_METAMORPHEUS_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamorpheus: 1.0.5
    END_VERSIONS
    """
}
