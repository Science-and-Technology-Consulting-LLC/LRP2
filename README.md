# LRP2: Long-Read Proteogenomics

Nextflow pipeline from PacBio Iso-Seq reads to filtered transcripts, a predicted proteome and multisample differential splicing and expression. Run it on Lifebit/CloudOS or any Nextflow with Docker.

## Layout

- `main.nf` reads the sample sheet and connects the modules with channels.
- `modules/<package>/<name>.nf`: one process per file, named for the process.
- `container/<package>/`: the image for that package (Dockerfile and context). Reference files and licences are inside their images.
- `nextflow.config`: settings, resources and the pinned image tag. `conf/` holds the GENCODE reference table and where results are published.
- `ipynb/`: Jupyter notebooks that draw the results and prepare IGV files (see Viewing the results).
- `docs/`: Sphinx documentation (`make -C docs html`).

Images are built by `.github/workflows/docker-publish.yml` and tagged with the commit SHA of the build. Set `--image_tag` to that SHA; never use `latest`.

## Stages

1. S1 PacBio Isocall: pbmm2 alignment, Isocall profile, prep, merge, call
2. S2 Transcriptome: SQANTI3 QC, hash ids, filtering
3. S3 Predicted proteome: CPAT, CPAT filter, SQANTI3 protein, protein classification
4. S4 Multisample: LeafCutter and edgeR / DRIMSeq, one comparison per pair of conditions

Proteomics (S5) is not wired yet; a sample sheet with protein rows stops with an error.

## Sample sheet

CSV with `sample_name`, `sample_path`, `condition`, `sample_type` (RNA). `sample_path` can be an `s3://` path. No spaces in `condition`. Each condition in a comparison needs at least two samples.

## Run

```bash
nextflow run main.nf --input samplesheet.csv --genome GRCh38.p14.v50
```

Reference: `--genome` (a key in `conf/gencode_references.config`) or both `--local_reference_gtf` and `--local_reference_fasta` with `--species`; not both. On Lifebit leave `outdir` as `results`.

Skip S1 with `--S2_custom_gtf` and `--S2_custom_counts`. Run S4 alone with `--S4_multisample_metadata`, `--S4_custom_gtf`, `--S4_custom_counts` and `--S4_custom_orf_counts` from an earlier run.

## Output

Under `results/`: `S1_PACBIO_ISOCALL`, `S2_TRANSCRIPTOME`, `S3_PREDICTED_PROTEOME`, `S4_MULTISAMPLE_ANALYSIS`, and `pipeline_info` (reports, trace, software versions).

## Viewing the results

The notebooks in `ipynb/` read the S3 and S4 files the pipeline wrote and do not change them. They run in an R or Python Jupyter session (on Lifebit, an interactive session with the results mounted read-only).

```bash
conda env create -f lrp2.yml
conda activate lrp2
python -m ipykernel install --user --name python3
Rscript -e 'IRkernel::installspec(name = "ir", displayname = "R")'
cd ipynb && jupyter lab
```

| Notebook | Kernel | Makes |
|---|---|---|
| `volcano.ipynb` | R | volcano plots of genes, transcript isoforms and ORFs |
| `heatmaps.ipynb` | R | clustered heatmaps of the significant genes, isoforms and ORFs, by group; isoforms and ORFs carry a known/novel colour bar |
| `gene_heatmaps.ipynb` | R | one heatmap per chosen gene (default HLA-B, ERAP1, ERAP2) of all its isoforms and ORFs |
| `specific_variants.ipynb` | R | isoforms, ORFs and splice variants found in one group and not the other: counts, heatmaps, tables in `tables/` |
| `dominant_orf.ipynb` | R | the dominant ORF (protein form) of each gene, its share, how many transcripts collapse into it, and whether it switches between groups |
| `usage.ipynb` | R | stacked bars of isoform and ORF shares per sample (DRIMSeq) |
| `splicing.ipynb` | R | volcano of splicing variants (change in PSI against cluster adjusted p) and the PSI matrix by sample, with the change between groups |
| `classification.ipynb` | R | SQANTI categories (FSM, ISM, NIC, NNC) for transcripts and proteins |
| `igv_files.ipynb` | Python | sorted, bgzip-compressed, indexed files in `igv/`, including tracks of the group-specific and significant features |

`ipynb/paths.json` is the one place that says where a run is: `run_dir`, the `comparison` name (prefix of the S4 file names), the `dataset_name` (prefix of the S3 file names) and the folders inside the run. To look at another run, change these and nothing else. Every figure is saved in `ipynb/figures/` as PNG and PDF, named with the comparison (`<comparison>_<name>.png`), and tables go to `ipynb/tables/`. The code that does the work is in `ipynb/src/`; each notebook cell only calls it.

## Credits

LRP2 builds on the original LRP pipeline (Miller et al. 2022), developed in the Sheynkman Lab. Anne Deslattes Mays and Gloria Sheynkman were the co-senior authors of that paper. We thank everyone who contributed to that work, which laid the foundation for this version.

### Development Team for LRP2

- **Megan D. Schertzer**, Sheynkman Lab - Lead developer
- **Julia T. Lewandowski**, Knowles Lab - Lead developer

- **David A. Knowles**, Development of LR LeafCutter and project support / funding 
- **Gloria Sheynkman**, Development/conceptualization of LRP and project support / funding

- **Emily F. Watts**, Sheynkman Lab - Contributions to LRP and conception of multi-sample analysis subworkflow.
- **Madison M. Mehlferber**, Sheynkman Lab - Contributor to the original LRP pipeline. Continued pipeline testing and feedback.
- **Will Rosenow**, Sheynkman Lab - Pipeline testing and feedback.
- **Scott I. Adamson**, Knowles Lab - Development of leafcutter-py. 
- **Jocelyne Bruand**, Pacific Biosciences - Development of Isocall. 
- **Elizabeth Tseng**, Pacific Biosciences - Development of Isocall. 
- **Egor Dolzhenko**, Pacific Biosciences - Lead Developer of Isocall. 


## License

This pipeline is released under the MIT License.

## Citations

### Citing LRP2

If you use LRP2 in your work, please cite:

> Schertzer, M. D., Lewandowski, J. T., et al. 2026. "LRP2: A proteogenomics pipeline for long-read informed protein isoform analysis and discovery." *bioRxiv*. doi: [10.64898/2026.05.27.728216](https://doi.org/10.64898/2026.05.27.728216)

LRP2 builds on the original LRP framework:

> Miller, R. M., Jordan, B. T., Mehlferber, M. M., et al. 2022. "Enhanced protein isoform characterization through long-read proteogenomics." *Genome Biology* 23(1): 69. doi: [10.1186/s13059-022-02624-y](https://doi.org/10.1186/s13059-022-02624-y)

