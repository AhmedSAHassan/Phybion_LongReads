#!/usr/bin/env python3
"""
diagnose_novel_junctions.py

Automates the IGV-by-eye checks we were doing manually:

  1. Exact match      - does short-read data (STAR SJ.out.tab) support the
                         junction at its EXACT novel coordinates?
  2. Nearby junction   - is there a DIFFERENT, well-supported SR junction
                         within +/- NEARBY_WINDOW bp? (boundary-shift signal)
  3. Locus activity    - is there ANY meaningful SR splicing happening within
                         +/- LOCUS_WINDOW bp, regardless of exact position?
                         (distinguishes "SR blind spot" from "SR active but
                         doesn't support this junction")
  4. LR consistency    - how many of your long-read BAM files independently
                         support the exact junction coordinates?
  5. Mappability       - optional: average mappability score in the region,
                         if you supply a bigWig track.

Output: one TSV row per junction, with a plain-language verdict column.

--------------------------------------------------------------------------
EDIT THE CONFIG SECTION BELOW, THEN RUN:

    python3 diagnose_novel_junctions.py

Requires: samtools on PATH, and the pysam Python package.
    pip install pysam --break-system-packages   (if not already installed)

Optional mappability check requires pyBigWig:
    pip install pyBigWig --break-system-packages
--------------------------------------------------------------------------
"""

import glob
import os
import re
import sys
from collections import defaultdict

import pysam

# ============================================================
# CONFIG - edit these paths/values for your setup
# ============================================================

GTF_PATH = "../novel_26_transcripts.gtf"

# The 19 transcripts with filter_result == "Isoform" from your
# classification.txt. Edit this list if your filter-passed set changes.
TRANSCRIPT_IDS = [
    "transcript25470.chr1.nic",
    "transcript34121.chr1.nic",
    "transcript61117.chr1.nnic",
    "transcript61119.chr1.nnic",
    "transcript30521.chr2.nnic",
    "transcript926.chr3.nnic",
    "transcript32230.chr3.nnic",
    "transcript29638.chr6.nnic",
    "transcript31443.chr6.nnic",
    "transcript21664.chr9.nnic",
    "transcript26863.chr12.nnic",
    "transcript26870.chr12.nnic",
    "transcript26885.chr12.nnic",
    "transcript33002.chr12.nnic",
    "transcript8598.chr15.nnic",
    "transcript13826.chr15.nnic",
    "transcript20940.chr17.nnic",
    "transcript39964.chr19.nnic",
    "transcript13844.chr22.nnic",
]

# Glob pattern matching your per-sample STAR SJ.out.tab files
SJ_GLOB = "../Star_aligned_sub/*/SJ.out.tab"

# Your long-read (ONT) BAM files - one per sample/replicate, sorted+indexed
LR_BAMS = sorted(glob.glob("../minimap_align/sorted_BAMs/*.bam"))  # <-- edit this glob

# Optional: path to a mappability bigWig (e.g. UCSC CRG Alignability,
# GEM, or Umap/Bismap hg38 track). Leave as None to skip this check.
MAPPABILITY_BW = None  # e.g. "../reference/hg38.k100.Umap.MultiTrackMappability.bw"

OUTPUT_TSV = "novel_junction_diagnostics.tsv"

# Thresholds
NEARBY_WINDOW = 500        # bp, for "boundary-shift" search
LOCUS_WINDOW = 5000        # bp, for "is this locus SR-active at all"
MIN_READS_FOR_ACTIVITY = 5  # min unique reads (any sample) to count as "active"
SMALL_OFFSET_THRESH = 100   # bp, offset below which we call it "boundary shift"
MIN_SAMPLES_FOR_CONFIRMED = 2  # min # SR samples w/ >=1 read to call "confirmed"

# ============================================================
# Step 1: parse junctions from GTF for the transcripts of interest
# ============================================================

def parse_gtf_junctions(gtf_path, transcript_ids):
    wanted = set(transcript_ids)
    exons = defaultdict(list)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9 or x[2] != "exon":
                continue
            m = re.search(r'transcript_id\s+"([^"]+)"', x[8])
            if not m:
                continue
            tx = m.group(1)
            if tx not in wanted:
                continue
            exons[tx].append((x[0], int(x[3]), int(x[4]), x[6]))

    junctions = []  # (tx, chrom, intron_start, intron_end, strand, junction_number)
    for tx, exon_list in exons.items():
        exon_list = sorted(exon_list, key=lambda z: z[1])
        if len(exon_list) < 2:
            continue
        for i, (e1, e2) in enumerate(zip(exon_list, exon_list[1:]), 1):
            if e1[0] != e2[0]:
                continue
            junctions.append((tx, e1[0], e1[2] + 1, e2[1] - 1, e1[3], i))

    found = {j[0] for j in junctions}
    missing = wanted - found
    if missing:
        print(f"WARNING: no multi-exon junctions found for: {sorted(missing)}",
              file=sys.stderr)

    return junctions


# ============================================================
# Step 2: load all SR SJ.out.tab files into memory, indexed by chrom
# ============================================================

def load_sj_data(sj_glob):
    sj_files = sorted(glob.glob(sj_glob))
    if not sj_files:
        print(f"ERROR: no SJ.out.tab files matched {sj_glob}", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(sj_files)} SJ.out.tab files.", file=sys.stderr)

    # per_sample[sample][chrom] = list of (start, end, unique_reads)
    per_sample = {}
    for sjfile in sj_files:
        sample = os.path.basename(os.path.dirname(sjfile))
        by_chrom = defaultdict(list)
        with open(sjfile) as f:
            for line in f:
                x = line.rstrip("\n").split("\t")
                if len(x) < 7:
                    continue
                chrom, start, end = x[0], int(x[1]), int(x[2])
                unique_reads = int(x[6])
                by_chrom[chrom].append((start, end, unique_reads))
        for chrom in by_chrom:
            by_chrom[chrom].sort()
        per_sample[sample] = by_chrom
    return per_sample


def exact_sr_support(sj_data, chrom, start, end):
    """Return dict sample -> unique_reads for an EXACT coordinate match."""
    out = {}
    for sample, by_chrom in sj_data.items():
        n = 0
        for s, e, reads in by_chrom.get(chrom, []):
            if s == start and e == end:
                n = reads
                break
        out[sample] = n
    return out


def nearest_sr_junction(sj_data, chrom, start, end, window):
    """
    Find the best-supported DIFFERENT junction within +/- window bp of the
    novel junction (start or end shifted). Returns (offset_start, offset_end,
    best_reads, best_sample) or None if nothing found.
    """
    best = None  # (total_reads, offset_start, offset_end, sample)
    for sample, by_chrom in sj_data.items():
        for s, e, reads in by_chrom.get(chrom, []):
            if s == start and e == end:
                continue  # that's the exact match, not "nearby"
            if abs(s - start) <= window or abs(e - end) <= window:
                if best is None or reads > best[0]:
                    best = (reads, s - start, e - end, sample)
    if best is None:
        return None
    reads, off_s, off_e, sample = best
    return {"offset_start": off_s, "offset_end": off_e,
            "reads": reads, "sample": sample}


def sr_locus_activity(sj_data, chrom, start, end, window, min_reads):
    """Is there ANY junction with >= min_reads within +/- window of the locus?"""
    lo, hi = start - window, end + window
    max_reads = 0
    for sample, by_chrom in sj_data.items():
        for s, e, reads in by_chrom.get(chrom, []):
            if lo <= s <= hi or lo <= e <= hi:
                max_reads = max(max_reads, reads)
    return max_reads >= min_reads, max_reads


# ============================================================
# Step 3: long-read support - count spliced reads matching the junction
# ============================================================

def lr_support(lr_bams, chrom, start, end, tolerance=2):
    """
    For each LR BAM, count reads with a spliced (N-op) alignment whose
    intron boundaries match (start,end) within `tolerance` bp.
    Returns dict bam_name -> read_count.
    """
    out = {}
    region_start = max(0, start - 50)
    region_end = end + 50
    for bam_path in lr_bams:
        name = os.path.basename(bam_path)
        count = 0
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                for read in bam.fetch(chrom, region_start, region_end):
                    if read.is_unmapped or read.cigartuples is None:
                        continue
                    ref_pos = read.reference_start
                    for op, length in read.cigartuples:
                        if op in (0, 2, 7, 8):  # M, D, =, X consume ref
                            ref_pos += length
                        elif op == 3:  # N = intron gap
                            intron_start = ref_pos + 1
                            intron_end = ref_pos + length
                            if (abs(intron_start - start) <= tolerance and
                                    abs(intron_end - end) <= tolerance):
                                count += 1
                            ref_pos += length
        except (OSError, ValueError) as ex:
            print(f"WARNING: could not read {bam_path} at {chrom}:{start}-{end}: {ex}",
                  file=sys.stderr)
            count = -1
        out[name] = count
    return out


# ============================================================
# Step 4: optional mappability check
# ============================================================

def mappability_score(bw_path, chrom, start, end, flank=50):
    if bw_path is None:
        return None
    try:
        import pyBigWig
    except ImportError:
        print("pyBigWig not installed - skipping mappability check.", file=sys.stderr)
        return None
    try:
        bw = pyBigWig.open(bw_path)
        val = bw.stats(chrom, max(0, start - flank), end + flank, type="mean")[0]
        bw.close()
        return round(val, 3) if val is not None else None
    except Exception as ex:
        print(f"WARNING: mappability lookup failed for {chrom}:{start}-{end}: {ex}",
              file=sys.stderr)
        return None


# ============================================================
# Verdict logic
# ============================================================

def make_verdict(exact_total_reads, exact_n_samples, nearby, locus_active,
                  lr_n_samples, lr_total_bams):
    if exact_n_samples >= MIN_SAMPLES_FOR_CONFIRMED:
        return "SR-CONFIRMED: real junction, supported at exact coordinates"

    if nearby is not None and abs(nearby["offset_start"]) <= SMALL_OFFSET_THRESH \
            and abs(nearby["offset_end"]) <= SMALL_OFFSET_THRESH and nearby["reads"] >= MIN_READS_FOR_ACTIVITY:
        return (f"LIKELY BOUNDARY-SHIFT of known SR junction "
                f"(offset {nearby['offset_start']}/{nearby['offset_end']} bp, "
                f"{nearby['reads']} reads)")

    if not locus_active:
        return "INCONCLUSIVE: SR blind spot in whole locus - check mappability before concluding"

    if lr_n_samples >= 3:
        return "CANDIDATE NOVEL: SR active in locus but absent at this junction, consistent LR support"

    return "WEAK: low LR consistency and no SR support - treat cautiously"


# ============================================================
# Main
# ============================================================

def main():
    print("Parsing GTF junctions...", file=sys.stderr)
    junctions = parse_gtf_junctions(GTF_PATH, TRANSCRIPT_IDS)
    print(f"Found {len(junctions)} junctions across {len(TRANSCRIPT_IDS)} transcripts.",
          file=sys.stderr)

    print("Loading SR SJ.out.tab data...", file=sys.stderr)
    sj_data = load_sj_data(SJ_GLOB)

    if not LR_BAMS:
        print("WARNING: no LR BAM files found via LR_BAMS glob - "
              "LR columns will be empty. Edit LR_BAMS in the config section.",
              file=sys.stderr)

    lr_bam_names = [os.path.basename(b) for b in LR_BAMS]

    rows = []
    header = [
        "transcript", "junction_number", "chrom", "start", "end", "strand",
        "SR_exact_total_reads", "SR_exact_n_samples_supporting",
        "SR_nearby_offset_start", "SR_nearby_offset_end", "SR_nearby_reads",
        "SR_locus_active", "SR_locus_max_reads",
        "LR_n_bams_supporting", "LR_n_bams_total", "LR_total_reads",
    ] + [f"LR_reads__{name}" for name in lr_bam_names] + [
        "mappability_score", "verdict",
    ]

    for i, (tx, chrom, start, end, strand, jnum) in enumerate(junctions, 1):
        print(f"[{i}/{len(junctions)}] {tx} junction_{jnum} "
              f"{chrom}:{start}-{end}", file=sys.stderr)

        exact = exact_sr_support(sj_data, chrom, start, end)
        exact_total = sum(exact.values())
        exact_n_samples = sum(1 for v in exact.values() if v > 0)

        nearby = nearest_sr_junction(sj_data, chrom, start, end, NEARBY_WINDOW)

        locus_active, locus_max_reads = sr_locus_activity(
            sj_data, chrom, start, end, LOCUS_WINDOW, MIN_READS_FOR_ACTIVITY)

        if LR_BAMS:
            lr = lr_support(LR_BAMS, chrom, start, end)
            lr_n_supporting = sum(1 for v in lr.values() if v and v > 0)
            lr_total_reads = sum(v for v in lr.values() if v and v > 0)
            # per-bam breakdown, in the same order as lr_bam_names
            lr_per_bam = [lr.get(name, "NA") for name in lr_bam_names]
        else:
            lr_n_supporting = 0
            lr_total_reads = 0
            lr_per_bam = ["NA"] * len(lr_bam_names)

        mapp = mappability_score(MAPPABILITY_BW, chrom, start, end)

        verdict = make_verdict(exact_total, exact_n_samples, nearby,
                                locus_active, lr_n_supporting, len(LR_BAMS))

        rows.append([
            tx, jnum, chrom, start, end, strand,
            exact_total, exact_n_samples,
            nearby["offset_start"] if nearby else "NA",
            nearby["offset_end"] if nearby else "NA",
            nearby["reads"] if nearby else "NA",
            locus_active, locus_max_reads,
            lr_n_supporting, len(LR_BAMS), lr_total_reads,
        ] + lr_per_bam + [
            mapp if mapp is not None else "NA",
            verdict,
        ])

    with open(OUTPUT_TSV, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    print(f"\nDone. Wrote {len(rows)} rows to {OUTPUT_TSV}", file=sys.stderr)


if __name__ == "__main__":
    main()