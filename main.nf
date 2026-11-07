#!/usr/bin/env nextflow
/*
 * LRP2: long-read proteogenomics
 *
 * main.nf reads the sample sheet, builds the channels and calls the modules.
 * Each module is one process in one container (container/<name>/), so nothing is read from projectDir.
 *
 *   S1  PacBio Isocall        pbmm2 align, Isocall profile / prep / merge / call
 *   S2  Transcriptome         SQANTI3 QC, hash ids, filter
 *   S3  Predicted proteome    CPAT, filter CPAT, SQANTI3 protein, protein classification
 *   S4  Multisample           LeafCutter, edgeR / DRIMSeq (one comparison per pair of conditions)
 *
 * Modes (chosen by which parameters are given):
 *   multisample only : --S4_multisample_metadata --S4_custom_gtf --S4_custom_counts --S4_custom_orf_counts
 *   skip Isocall     : --S2_custom_gtf --S2_custom_counts   (non-PacBio data)
 */

include { PBMM2_ALIGN as ISOCALL_ALIGN                                  } from './modules/pbmm2/pbmm2_align'
include { ISOCALL_PROFILE                                               } from './modules/isocall/isocall_profile'
include { ISOCALL_PREP                                                  } from './modules/isocall/isocall_prep'
include { ISOCALL_MERGE                                                 } from './modules/isocall/isocall_merge'
include { ISOCALL_CALL                                                  } from './modules/isocall/isocall_call'
include { GUNZIP as GUNZIP_GTF                                          } from './modules/ubuntu/gunzip'
include { GUNZIP as GUNZIP_FASTA                                        } from './modules/ubuntu/gunzip'
include { GZIP as GZIP_GTF                                              } from './modules/ubuntu/gzip'
include { SANITIZE_PAR_IDS                                              } from './modules/lrp2_lite/sanitize_par_ids'
include { SANITIZE_TRANSCRIPT_IDS as SANITIZE_TRANSCRIPT_IDS_S2_GTF     } from './modules/ubuntu/sanitize_transcript_ids'
include { SANITIZE_TRANSCRIPT_IDS as SANITIZE_TRANSCRIPT_IDS_S2_COUNTS  } from './modules/ubuntu/sanitize_transcript_ids'
include { VALIDATE_TRANSCRIPT_IDS as VALIDATE_TRANSCRIPT_IDS_S2         } from './modules/ubuntu/validate_transcript_ids'
include { SQANTI_QC                                                     } from './modules/sqanti3/sqanti_qc'
include { GENERATE_HASHIDS                                              } from './modules/lrp2_lite/generate_hashids'
include { FILTER_TRANSCRIPTOME                                          } from './modules/lrp2_lite/filter_transcriptome'
include { CPAT_ORF                                                      } from './modules/cpat/cpat_orf'
include { FILTER_CPAT                                                   } from './modules/lrp2_lite/filter_cpat'
include { SQANTI_PROTEIN                                                } from './modules/sqanti3/sqanti_protein'
include { PROTEIN_CLASSIFICATION                                        } from './modules/lrp2_lite/protein_classification'
include { LEAFCUTTER_LONGREAD                                           } from './modules/leafcutter_longread/leafcutter_longread'
include { DIFFERENTIAL_EXPRESSION                                       } from './modules/lrp2_lite/differential_expression'

/*
 * Sample sheet: one CSV line to fields (handles "quoted,fields")
 */
def parseCsvLine(String line) {
    return line.split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/, -1).collect { it.trim().replaceAll(/^"|"$/, '') }
}

/*
 * Sample sheet text without a byte order mark or Windows line endings (Excel adds them).
 * Nextflow file() reads local paths and s3://, so the sheet can sit in a bucket.
 */
def sampleSheetText(input_file) {
    return file(input_file).getText('UTF-8')
        .replaceAll('^﻿', '')
        .replaceAll('\r\n', '\n')
        .replaceAll('\r', '\n')
}

/*
 * Sample sheet to a list of [meta, sample_path]. Required columns: sample_name, sample_path, condition, sample_type.
 */
def readSamplesheet(input_file) {
    def lines = sampleSheetText(input_file).readLines().findAll { it.trim() }
    if (lines.size() < 2) {
        error "ERROR: samplesheet '${input_file}' has no samples."
    }
    def header = parseCsvLine(lines[0])
    def required = ['sample_name', 'sample_path', 'condition', 'sample_type']
    def missing = required.findAll { !(it in header) }
    if (missing) {
        error "ERROR: samplesheet is missing required column(s): ${missing.join(', ')}. Required: ${required.join(', ')}."
    }
    def ext_pattern = ~/(?i).*\.(bam|fastq|fastq\.gz|fq|fq\.gz|raw|mzml|fasta)(\?\S*)?$/
    def result = []
    lines.drop(1).eachWithIndex { line, i ->
        def cols = parseCsvLine(line)
        def row = [:]
        header.eachWithIndex { h, j -> row[h] = j < cols.size() ? cols[j] : '' }
        def where = "samplesheet row ${i + 2} (${row.sample_name ?: 'no sample_name'})"
        required.each { c ->
            if (!row[c]) error "ERROR: ${where}: '${c}' is required."
        }
        ['sample_name', 'sample_path', 'condition'].each { c ->
            if (row[c] =~ /\s/) error "ERROR: ${where}: '${c}' must not contain spaces."
        }
        if (!(row.sample_type.toLowerCase() in ['rna', 'protein'])) {
            error "ERROR: ${where}: sample_type must be 'RNA' or 'protein', not '${row.sample_type}'."
        }
        if (!(row.sample_path ==~ ext_pattern)) {
            error "ERROR: ${where}: sample_path must end in .bam, .fastq(.gz), .fq(.gz), .raw, .mzML or .fasta."
        }
        result << [[sample_name: row.sample_name, condition: row.condition, sample_type: row.sample_type], row.sample_path]
    }
    return result
}

/*
 * Samples per condition from a metadata table (columns name|sample_name and group|condition).
 */
def samplesPerCondition(String metadata_text) {
    def lines = metadata_text.readLines().findAll { it.trim() }
    def header = parseCsvLine(lines[0])
    def group_idx = header.indexOf('group') >= 0 ? header.indexOf('group') : header.indexOf('condition')
    def name_idx  = header.indexOf('name')  >= 0 ? header.indexOf('name')  : header.indexOf('sample_name')
    def type_idx  = header.indexOf('sample_type')
    def result = [:].withDefault { [] }
    lines.drop(1).each { line ->
        def cols = parseCsvLine(line)
        if (cols.size() > group_idx && (type_idx < 0 || cols[type_idx].toLowerCase() == 'rna')) {
            result[cols[group_idx]] << cols[name_idx]
        }
    }
    return result
}

workflow {

    main:
    if (params.help) {
        log.info "nextflow run ${workflow.manifest.name} --input samplesheet.csv --outdir <OUTDIR>  (see README.md for all parameters)"
        exit 0
    }

    def multisample_only = params.S4_custom_gtf && params.S4_custom_counts && params.S4_custom_orf_counts && params.S4_multisample_metadata
    def skip_isocall     = params.S2_custom_gtf && params.S2_custom_counts

    if (skip_isocall && !(params.S2_custom_gtf && params.S2_custom_counts)) {
        error "ERROR: --S2_custom_gtf and --S2_custom_counts must be given together."
    }
    if (!multisample_only && !params.input) {
        error "ERROR: --input (the sample sheet) is required, unless --S4_multisample_metadata, --S4_custom_gtf, --S4_custom_counts and --S4_custom_orf_counts are all given."
    }
    if (params.species && !['human', 'mouse'].contains(params.species)) {
        error "ERROR: --species must be 'human' or 'mouse', not '${params.species}'."
    }
    // references: --genome (GENCODE) or both --local_reference_gtf and --local_reference_fasta
    def uses_local_refs = params.local_reference_gtf || params.local_reference_fasta
    if (uses_local_refs && !(params.local_reference_gtf && params.local_reference_fasta)) {
        error "ERROR: give --local_reference_gtf and --local_reference_fasta together, or use --genome."
    }
    if (uses_local_refs && workflow.commandLine.contains('--genome ')) {
        error "ERROR: --genome cannot be combined with --local_reference_gtf / --local_reference_fasta."
    }
    if (!uses_local_refs && params.genome && params.genomes && !params.genomes.containsKey(params.genome)) {
        error "ERROR: genome '${params.genome}' is not known. Available: ${params.genomes.keySet().join(', ')}"
    }
    if (!multisample_only && !params.gtf && !params.fasta) {
        error "ERROR: no reference. Give --genome, or --local_reference_gtf and --local_reference_fasta."
    }
    if (!multisample_only && uses_local_refs && !params.species) {
        error "ERROR: --species (human or mouse) is required with local reference files."
    }

    ch_versions = channel.empty()

    // ---------------------------------------------------------------- the sample sheet
    def samples = multisample_only ? [] : readSamplesheet(params.input)
    def rna_samples = samples.findAll { it[0].sample_type.toLowerCase() == 'rna' }
    if (samples.any { it[0].sample_type.toLowerCase() == 'protein' }) {
        error "ERROR: the sample sheet has protein rows. The proteomics steps (S5) are not wired into this version of main.nf; remove the protein rows."
    }
    if (!multisample_only && rna_samples.isEmpty()) {
        error "ERROR: no RNA samples in the sample sheet."
    }
    def rna_names = rna_samples.collect { it[0].sample_name }
    def duplicates = rna_names.findAll { n -> rna_names.count(n) > 1 }.unique()
    if (duplicates) {
        error "ERROR: duplicate sample_name values: ${duplicates.join(', ')}. Each RNA sample needs its own sample_name."
    }

    def metadata_text = multisample_only ? file(params.S4_multisample_metadata).getText('UTF-8') : sampleSheetText(params.input)
    def per_condition = samplesPerCondition(metadata_text)
    def valid_conditions = per_condition.findAll { c, s -> s.size() >= 2 }.keySet().toList().sort()
    def pairs = []
    valid_conditions.eachWithIndex { c1, i -> valid_conditions.eachWithIndex { c2, j -> if (j > i) pairs << [c1, c2] } }

    log.info "LRP2: ${multisample_only ? 'multisample only' : rna_samples.size() + ' RNA samples'}; conditions: ${per_condition.collect { c, s -> c + ' (' + s.size() + ')' }.join(', ')}"
    if (pairs.isEmpty()) {
        log.warn "Fewer than two conditions with two or more samples: the multisample comparison (S4) will not run."
    } else {
        log.info "S4 comparisons: ${pairs.collect { it[0] + '_vs_' + it[1] }.join(', ')}"
    }

    // ---------------------------------------------------------------- S1 and S2: Isocall, SQANTI3, filter
    if (!multisample_only) {

        // reference GTF: unzip if needed, fix _PAR_Y gene ids if present, gzip for Isocall
        def gtf_file = file(params.gtf)
        ch_gtf_in = channel.of([[id: 'gtf'], gtf_file])
        if (gtf_file.name.endsWith('.gz')) {
            GUNZIP_GTF(ch_gtf_in)
            ch_gtf_plain = GUNZIP_GTF.out.gunzip
        } else {
            ch_gtf_plain = ch_gtf_in
        }
        ch_gtf_plain
            .branch { meta, gtf ->
                needs_sanitization: gtf.text.contains('_PAR_Y')
                no_sanitization: true
            }
            .set { ch_gtf_branched }
        SANITIZE_PAR_IDS(ch_gtf_branched.needs_sanitization)
        ch_gtf_final = ch_gtf_branched.no_sanitization.mix(SANITIZE_PAR_IDS.out.sanitized_gtf)
        GZIP_GTF(ch_gtf_final)
        ch_gtf    = ch_gtf_final.map { _meta, f -> f }.first()
        ch_gtf_gz = GZIP_GTF.out.gzip.map { _meta, f -> f }.first()

        // reference FASTA
        def fasta_file = file(params.fasta)
        if (fasta_file.name.endsWith('.gz')) {
            GUNZIP_FASTA(channel.of([[id: 'fasta'], fasta_file]))
            ch_fasta = GUNZIP_FASTA.out.gunzip.map { _meta, f -> f }.first()
        } else {
            ch_fasta = channel.value(fasta_file)
        }

        if (!skip_isocall) {
            // S1: one alignment and profile per sample, then merge and call
            ch_rna = channel.fromList(rna_samples.collect { row ->
                [[id: row[0].sample_name, sample_name: row[0].sample_name, condition: row[0].condition, sample_type: 'RNA'], file(row[1])]
            })
            ISOCALL_ALIGN(ch_rna, ch_fasta)
            ISOCALL_PROFILE(ISOCALL_ALIGN.out.bam)
            ISOCALL_PREP(ch_gtf_gz)
            ch_all_profiles = ISOCALL_PROFILE.out.profile
                .map { _meta, profile -> profile }
                .collect()
                .map { profiles -> [[id: params.dataset_name], profiles] }
            ISOCALL_MERGE(ch_all_profiles)
            ISOCALL_CALL(ISOCALL_MERGE.out.merged_profile, ISOCALL_PREP.out.isoforms, ch_fasta)
            ch_called_gtf    = ISOCALL_CALL.out.gtf
            ch_count_matrix  = ISOCALL_CALL.out.count_matrix
            ch_versions = ch_versions.mix(ISOCALL_ALIGN.out.versions, ISOCALL_PROFILE.out.versions, ISOCALL_PREP.out.versions, ISOCALL_MERGE.out.versions, ISOCALL_CALL.out.versions)
        } else {
            // not PacBio: the transcript GTF and counts come from the user
            SANITIZE_TRANSCRIPT_IDS_S2_GTF(channel.of([[id: params.dataset_name], file(params.S2_custom_gtf)]))
            SANITIZE_TRANSCRIPT_IDS_S2_COUNTS(channel.of([[id: params.dataset_name], file(params.S2_custom_counts)]))
            VALIDATE_TRANSCRIPT_IDS_S2(SANITIZE_TRANSCRIPT_IDS_S2_GTF.out.sanitized.join(SANITIZE_TRANSCRIPT_IDS_S2_COUNTS.out.sanitized, by: 0))
            ch_called_gtf   = VALIDATE_TRANSCRIPT_IDS_S2.out.validated.map { meta, gtf, counts -> [meta, gtf] }
            ch_count_matrix = VALIDATE_TRANSCRIPT_IDS_S2.out.validated.map { meta, gtf, counts -> [meta, counts] }
            ch_versions = ch_versions.mix(SANITIZE_TRANSCRIPT_IDS_S2_GTF.out.versions, SANITIZE_TRANSCRIPT_IDS_S2_COUNTS.out.versions, VALIDATE_TRANSCRIPT_IDS_S2.out.versions)
        }

        // S2: SQANTI3 QC, hash ids, filter
        SQANTI_QC(ch_called_gtf.join(ch_count_matrix, by: 0), ch_gtf, ch_fasta)
        GENERATE_HASHIDS(SQANTI_QC.out.corrected_gtf.join(SQANTI_QC.out.classification, by: 0), ch_gtf)
        FILTER_TRANSCRIPTOME(
            SQANTI_QC.out.classification
                .join(SQANTI_QC.out.corrected_gtf, by: 0)
                .join(SQANTI_QC.out.corrected_fasta, by: 0)
                .join(GENERATE_HASHIDS.out.hashids_mapping, by: 0)
                .join(GENERATE_HASHIDS.out.corrected_psl, by: 0),
            ch_gtf,
            ch_fasta,
            file(params.sample_metadata ?: params.input)
        )

        // S3: CPAT open reading frames, filter, SQANTI3 protein, classification
        ch_transcriptome = FILTER_TRANSCRIPTOME.out.corrected_fasta_filtered
            .join(FILTER_TRANSCRIPTOME.out.corrected_gtf_filtered, by: 0)
            .join(FILTER_TRANSCRIPTOME.out.classification_filtered, by: 0)
            .join(FILTER_TRANSCRIPTOME.out.hashids_filtered, by: 0)
            .join(GENERATE_HASHIDS.out.hashids_mapping, by: 0)
        CPAT_ORF(ch_transcriptome.map { meta, fasta, gtf, classification, hashids, mapping -> [meta, fasta] })
        FILTER_CPAT(
            CPAT_ORF.out.orf_prob
                .join(CPAT_ORF.out.orf_seqs, by: 0)
                .join(ch_transcriptome, by: 0)
                .map { meta, orf_prob, orf_seqs, fasta, gtf, _classification, hashids, mapping -> [meta, orf_prob, orf_seqs, fasta, gtf, mapping] },
            ch_gtf
        )
        SQANTI_PROTEIN(FILTER_CPAT.out.cds_gtf, ch_gtf)
        PROTEIN_CLASSIFICATION(
            SQANTI_PROTEIN.out.protein_classification
                .join(FILTER_CPAT.out.cds_gtf, by: 0)
                .join(ch_transcriptome, by: 0)
                .join(FILTER_CPAT.out.all_orfs_mapped, by: 0)
                .map { meta, protein_class, cds_gtf, fasta, gtf, classification, hashids, mapping, all_orfs -> [meta, protein_class, cds_gtf, fasta, hashids, all_orfs] },
            ch_gtf
        )
        ch_versions = ch_versions.mix(SQANTI_QC.out.versions, GENERATE_HASHIDS.out.versions, FILTER_TRANSCRIPTOME.out.versions,
            CPAT_ORF.out.versions, FILTER_CPAT.out.versions, SQANTI_PROTEIN.out.versions, PROTEIN_CLASSIFICATION.out.versions)

        ch_transcripts = FILTER_TRANSCRIPTOME.out.corrected_gtf_filtered.join(FILTER_TRANSCRIPTOME.out.hashids_filtered, by: 0)
        ch_orfs        = PROTEIN_CLASSIFICATION.out.hashids_orf
    }

    // ---------------------------------------------------------------- S4: multisample comparisons
    if (multisample_only || !pairs.isEmpty()) {

        if (multisample_only) {
            ch_transcripts = channel.of([[id: params.dataset_name], file(params.S4_custom_gtf), file(params.S4_custom_counts)])
            ch_orfs        = channel.of([[id: params.dataset_name], file(params.S4_custom_orf_counts)])
            ch_metadata    = channel.value(file(params.S4_multisample_metadata))
        } else {
            // RNA rows only, with the column names the S4 scripts expect
            def rna_lines = ['name,sample_path,group,sample_type'] + rna_samples.collect { "${it[0].sample_name},${it[1]},${it[0].condition},RNA" }
            ch_metadata = channel.of(rna_lines.join('\n') + '\n')
                .collectFile(name: "${params.dataset_name}_rna_samples_metadata.csv")
                .first()
        }

        ch_pairs = channel.fromList(pairs).map { c1, c2 -> [c1, c2, "${c1}_vs_${c2}"] }

        ch_transcripts_pairs = ch_transcripts
            .map { meta, gtf, counts -> [gtf, counts] }
            .combine(ch_pairs)
            .map { gtf, counts, control, experimental, comparison_id -> [[id: comparison_id], gtf, counts, control, experimental] }
            .multiMap { meta, gtf, counts, control, experimental ->
                data: [meta, gtf, counts]
                control: control
                experimental: experimental
            }

        LEAFCUTTER_LONGREAD(
            ch_transcripts_pairs.data,
            ch_metadata,
            ch_transcripts_pairs.control,
            ch_transcripts_pairs.experimental,
            params.min_samples_per_intron,
            params.min_samples_per_group,
            params.min_usage_ratio,
            params.dataset_name,
            params.leafcutter_threads
        )

        ch_counts_pairs = ch_transcripts
            .map { meta, gtf, counts -> counts }
            .combine(ch_orfs.map { meta, orf_counts -> orf_counts })
            .combine(ch_pairs)
            .map { counts, orf_counts, control, experimental, comparison_id -> [[id: comparison_id], counts, orf_counts, control, experimental] }
            .multiMap { meta, counts, orf_counts, control, experimental ->
                data: [meta, counts, orf_counts]
                control: control
                experimental: experimental
            }

        DIFFERENTIAL_EXPRESSION(
            ch_counts_pairs.data,
            ch_metadata,
            ch_counts_pairs.control,
            ch_counts_pairs.experimental,
            params.drimseq_min_gene_expr,
            params.drimseq_min_isoform_prop
        )
        ch_versions = ch_versions.mix(LEAFCUTTER_LONGREAD.out.versions, DIFFERENTIAL_EXPRESSION.out.versions)
    }

    // ---------------------------------------------------------------- software versions
    ch_versions
        .unique()
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'lrp2_software_versions.yml', sort: true, newLine: true)

    onComplete:
    log.info(workflow.success ? "LRP2 finished: ${params.outdir}" : "LRP2 failed: ${workflow.errorMessage ?: 'see .nextflow.log'}")
}
