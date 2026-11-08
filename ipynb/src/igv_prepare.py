"""igv_prepare.py

Make the files IGV needs: sorted, compressed and indexed GTF and BED, and
indexed BAM. IGV reads a bgzip file with its .tbi index next to it.
"""

import gzip
import os
import re

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


def check_index(gz_path, region):
    """Number of records the index returns for a region such as 'chr22:1-50000000'."""
    with pysam.TabixFile(gz_path) as tbx:
        return sum(1 for _ in tbx.fetch(region=region))
