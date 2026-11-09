# ipynb

Graphics and IGV files for the results of an LRP2 run. The notebooks only read the S4 and S3 files the pipeline wrote; they do not change them.

## Environment

```bash
conda env create -f lrp2.yml
conda activate lrp2
python -m ipykernel install --user --name python3
Rscript -e 'IRkernel::installspec(name = "ir", displayname = "R")'
cd ipynb && jupyter lab
```

## Notebooks

| Notebook | Kernel | Reads | Makes |
|---|---|---|---|
| `volcano.ipynb` | R | `*_DGE_edgeR_results.txt`, `*_DTE_edgeR_results.txt` | volcano plots of genes and transcript isoforms |
| `heatmaps.ipynb` | R | `*_DTE_edgeR_normalized_CPM_matrix.txt`, `*_DE_ORF_edgeR_normalized_CPM_matrix.txt`, `*.lr_leafcutter.groups_file.txt` | clustered heatmaps by group |
| `gene_heatmaps.ipynb` | R | the isoform and ORF CPM matrices (as `heatmaps.ipynb`) | one heatmap per chosen gene (default HLA-B, ERAP1, ERAP2) of all its isoforms and ORFs, with an average-CPM bar |
| `specific_variants.ipynb` | R | `*_DTE_edgeR_raw_CPM_matrix.txt`, `*_DE_ORF_edgeR_raw_CPM_matrix.txt`, `*.subisoform_clusters.txt`, groups file | isoforms, ORFs and splice variants found in one group and not the other: counts, bar chart, heatmaps, tables in `tables/` |
| `usage.ipynb` | R | `*_DTU_transcript_DRIMSeq_summary.txt`, `*_DU_ORF_DRIMSeq_summary.txt`, groups file | stacked bars of isoform and ORF shares per sample, for the most significant genes |
| `splicing.ipynb` | R | `*.lr_leafcutter.subisoform_clusters.txt`, `*.ds_cluster_significance.txt`, `*.ds_effect_sizes.txt`, groups file | PSI matrix of splicing variants by sample, with the change between groups |
| `classification.ipynb` | R | `*.best_ORF_SQANTI_classification.tsv` | bar chart of SQANTI categories (FSM, ISM, NIC, NNC) for transcripts and proteins |
| `igv_files.ipynb` | Python | `*collapsed_high_confidence_ORF.gtf` and `.bed` | sorted, bgzip-compressed, indexed files in `igv/`, including coloured tracks of the group-specific features from `specific_variants.ipynb` |

## Where the results are

`paths.json` is the one place to set it. It holds the run's results folder (`run_dir`), the `comparison` name (the prefix of the S4 file names, for example `SpA_vs_not_SpA`), the `dataset_name` (the prefix of the S3 files), and the folders inside the run as short relative names (`leafcutter`, `de_gene`, `de_transcript`, `de_orf`, `du_transcript`, `du_orf`, `sqanti_protein`, `protein_classification`, `filter_transcriptome`). Every notebook starts with `loc <- run_locations()` (Python: `loc = run_locations()`) and uses names such as `loc$leafcutter`. To look at another run, such as SpA B27 against not B27, change `comparison` (and `dataset_name`, and `run_dir` if the run is elsewhere) in `paths.json` and nothing else.

## Saved figures

Every figure is also saved in `figures/` as PNG and PDF, named with the comparison: `figures/<comparison>_<name>.png`, for example `SpA_B27neg_vs_SpA_B27pos_volcano_genes.png`. The comparison comes from `paths.json`, so figures from different runs do not overwrite each other. `igv_files.ipynb` writes to `igv/`, whose files already carry the dataset name.

## Functions

The code that does the work is in `src/` (`volcano.R`, `heatmap.R`, `usage.R`, `splicing.R`, `classification.R`, `igv_prepare.py`); the notebook cells only call it.
