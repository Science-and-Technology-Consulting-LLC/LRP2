"""igv_prepare.py

Make the files IGV needs: sorted, compressed and indexed GTF and BED, and
indexed BAM. IGV reads a bgzip file with its .tbi index next to it.
"""

import gzip
import os
import re

import pandas as pd
import pysam


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _sort_lines(lines, chrom_col=0, start_col=3):
    """Sort data lines by chromosome name then start position; the header lines stay on top."""
    is_header = lambda l: l.startswith(("#", "track", "browser"))
    header = [l for l in lines if is_header(l)]
    body = [l for l in lines if l.strip() and not is_header(l)]
    def key(l):
        f = l.split("\t")
        return f[chrom_col], int(f[start_col])
    return header + sorted(body, key=key)


def index_annotation(path, out_dir, kind):
    """Sort `path` (kind 'gff' for GTF/GFF, 'bed' for BED), write it bgzip compressed into
    `out_dir` and index it. Returns the path of the .gz file."""
    if kind not in ("gff", "bed"):
        raise ValueError("kind must be 'gff' or 'bed'")
    os.makedirs(out_dir, exist_ok=True)
    start_col = 3 if kind == "gff" else 1
    with _open(path) as fh:
        lines = [l.rstrip("\n") + "\n" for l in fh]
    # tabix treats only '#' lines as comments, so BED 'track' and 'browser' lines are dropped
    lines = [l for l in lines if not l.startswith(("track", "browser"))]
    sorted_lines = _sort_lines(lines, 0, start_col)
    plain = os.path.join(out_dir, re.sub(r"\.gz$", "", os.path.basename(path)))
    with open(plain, "w") as out:
        out.writelines(sorted_lines)
    # tabix_index compresses with bgzip and writes plain.gz.tbi next to it
    return pysam.tabix_index(plain, preset=kind, force=True)


def psl_to_bed12(psl_path, bed_path):
    """Convert a PSL file (one aligned transcript per line) to BED12 so it can be sorted and indexed.
    Returns the number of transcripts written."""
    n = 0
    with open(psl_path) as fh, open(bed_path, "w") as out:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 21 or not f[0].isdigit():
                continue
            strand, name, chrom = f[8], f[9], f[13]
            start, end, n_blocks = int(f[15]), int(f[16]), int(f[17])
            sizes = f[18].rstrip(",")
            t_starts = [int(x) for x in f[20].rstrip(",").split(",")]
            rel = ",".join(str(x - start) for x in t_starts)
            out.write("\t".join([chrom, str(start), str(end), name, "0", strand, str(start), str(end),
                                 "0", str(n_blocks), sizes, rel]) + "\n")
            n += 1
    return n


def prepare_bam(path, out_dir):
    """Sort a BAM by position into `out_dir` and index it (writes the .bam and .bam.bai).
    Returns the path of the sorted BAM."""
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, re.sub(r"\.bam$", "", os.path.basename(path)) + ".sorted.bam")
    pysam.sort("-o", out, path)
    pysam.index(out)
    return out


def check_index(gz_path, region=None):
    """Number of records the index returns for a region such as 'chr22:1-50000000'.
    With no region, the first chromosome in the index is used."""
    with pysam.TabixFile(gz_path) as tbx:
        region = region or tbx.contigs[0]
        return sum(1 for _ in tbx.fetch(region=region))


def orf_index(gtf_path, hashids_path, de_orf_path=None):
    """One row per predicted ORF, to look up what to find in IGV. From the ORF GTF: gene, position
    (`locus` can be pasted into IGV's search box), strand, number of coding exons, coding length and the
    ORF ratio. From the hash-ids table: every transcript isoform that gives this ORF. If `de_orf_path`
    (the *_DE_ORF_edgeR_results.txt table) is given, logFC and FDR are added."""
    rows, cds = {}, {}
    with _open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "CDS"):
                continue
            attr = {k: a or b for k, a, b in re.findall(r'(\w+) (?:"([^"]*)"|([^;\s]+))', f[8])}   # values with or without quotes
            orf = attr["orf_base_id"]
            if f[2] == "transcript":
                rows[orf] = {"orf_base_id": orf, "gene": attr.get("reference_gene_name"),
                             "gene_id": attr.get("gene_id"), "chrom": f[0], "start": int(f[3]),
                             "end": int(f[4]), "strand": f[6],
                             "orf_ratio": float(attr.get("avg_orf_ratio", "nan"))}
            else:
                n, bp = cds.get(orf, (0, 0))
                cds[orf] = (n + 1, bp + int(f[4]) - int(f[3]) + 1)
    idx = pd.DataFrame(rows.values())
    idx["coding_exons"] = idx["orf_base_id"].map(lambda o: cds.get(o, (0, 0))[0])
    idx["coding_bp"] = idx["orf_base_id"].map(lambda o: cds.get(o, (0, 0))[1])
    idx["locus"] = idx["chrom"] + ":" + idx["start"].astype(str) + "-" + idx["end"].astype(str)
    iso = pd.read_csv(hashids_path, sep="\t", usecols=["orf_base_id", "orf_all_isoform_id"])
    iso["n_isoforms"] = iso["orf_all_isoform_id"].str.count(",") + 1
    idx = idx.merge(iso, on="orf_base_id", how="left")
    if de_orf_path:
        de = pd.read_csv(de_orf_path, sep="\t", usecols=["orf_base_id", "logFC", "FDR"])
        idx = idx.merge(de, on="orf_base_id", how="left")
    return idx.sort_values(["gene", "start"]).reset_index(drop=True)


def write_specific_bed(table_path, out_path, level, bed_in=None):
    """Write a BED of the features marked specific to a group in a table made by specific.R,
    coloured by group (the table's `colour` column). level is 'isoforms' or 'orfs' (then `bed_in` is the
    BED12 of the transcripts or ORFs, whose lines are kept) or 'variants' (one line per intron of each
    splicing variant, from the table itself). The result is not indexed: pass it to index_annotation.
    Returns the number of lines written."""
    tab = pd.read_csv(table_path, sep="\t", dtype=str, keep_default_na=False)
    tab = tab[tab["specific_to"] != "not specific"]
    n = 0
    with open(out_path, "w") as out:
        if level == "variants":
            for r in tab.to_dict("records"):
                for s, e in zip(r["intron_starts"].split(","), r["intron_ends"].split(",")):
                    out.write("\t".join([r["chr"], s, e, r["gene"] + " [" + r["specific_to"] + "]", "0", r["strand"],
                                         s, e, r["colour"]]) + "\n")
                    n += 1
        else:
            info = {r["id"]: r for r in tab.to_dict("records")}
            with _open(bed_in) as fh:
                for line in fh:
                    if line.startswith(("track", "browser", "#")):
                        continue
                    f = line.rstrip("\n").split("\t")
                    key = f[3].split("|")[0]
                    if level == "isoforms":
                        key = re.sub(r"_ENSG[^_]*$", "", key)   # the BED name is <isoform>_<gene id>
                    if key in info and len(f) >= 12:
                        f[3] = key + " [" + info[key]["specific_to"] + "]"
                        f[8] = info[key]["colour"]
                        out.write("\t".join(f) + "\n")
                        n += 1
    return n
