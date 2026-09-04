"""
Build a real transcriptome-background gene list for enrichment analysis:
genes reliably detected (CPM >= 1 in >= 20% of samples) across ALL real
TCGA-BRCA Primary Tumor RNA-seq (STAR - Counts) profiles from GDC.

This is a one-time derivation; the output (transcriptome_background_genes.txt)
is what the notebook actually loads, so the notebook doesn't need network
access or the (deleted) full raw-counts cohort to reproduce the enrichment step.

Downloads in resumable batches (cached under data/.background_cache/) since
the full cohort is ~1,100 samples / ~4.6GB -- a batch that already completed
and was cached is skipped on a re-run, so an interrupted run can pick back up
without re-downloading everything.
"""
import io
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).parent
CACHE_DIR = HERE / ".background_cache"
OUT_PATH = HERE / "transcriptome_background_genes.txt"
CPM_THRESHOLD = 1.0
DETECTION_FRACTION = 0.20
BATCH_SIZE = 50


def get_all_sample_file_ids():
    filt = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
        ],
    }
    params = {"filters": json.dumps(filt), "fields": "file_id", "size": "2000", "format": "json"}
    r = requests.get("https://api.gdc.cancer.gov/files", params=params, timeout=60)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    return [h["file_id"] for h in hits]


def download_batch(file_ids, max_retries=4):
    for attempt in range(max_retries):
        try:
            resp = requests.post(
                "https://api.gdc.cancer.gov/data",
                json={"ids": file_ids},
                headers={"Content-Type": "application/json"},
                timeout=300,
            )
            resp.raise_for_status()
            return resp.content
        except (requests.RequestException,) as e:
            wait = 5 * (2 ** attempt)
            print(f"    batch download failed ({e}); retrying in {wait}s ({attempt + 1}/{max_retries})")
            import time
            time.sleep(wait)
    raise RuntimeError(f"Batch download failed after {max_retries} retries")


def parse_batch(tar_bytes):
    """Returns (detected_counts: Series indexed by gene_name, n_parsed: int)."""
    detected_counts = None
    n_parsed = 0
    with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tf:
        for member in tf.getmembers():
            if not member.name.endswith(".tsv"):
                continue
            content = tf.extractfile(member).read().decode("utf-8")
            df = pd.read_csv(io.StringIO(content), sep="\t", skiprows=1)
            df = df[~df["gene_id"].str.startswith("N_")]
            df = df[df["gene_type"] == "protein_coding"]
            counts = df.set_index("gene_name")["unstranded"]
            counts = counts[~counts.index.duplicated(keep="first")]
            cpm = counts / counts.sum() * 1e6
            detected = (cpm >= CPM_THRESHOLD).astype(int)

            if detected_counts is None:
                detected_counts = detected
            else:
                detected_counts = detected_counts.add(detected, fill_value=0)
            n_parsed += 1
    return detected_counts, n_parsed


def main():
    CACHE_DIR.mkdir(exist_ok=True)

    print("Querying all TCGA-BRCA Primary Tumor STAR-counts file IDs...")
    file_ids = get_all_sample_file_ids()
    print(f"Got {len(file_ids)} file IDs total.")

    batches = [file_ids[i:i + BATCH_SIZE] for i in range(0, len(file_ids), BATCH_SIZE)]
    print(f"Split into {len(batches)} batches of up to {BATCH_SIZE} files each.")

    total_detected = None
    total_parsed = 0
    total_genes_seen = None

    for i, batch_ids in enumerate(batches):
        cache_file = CACHE_DIR / f"batch_{i:04d}.parquet"
        meta_file = CACHE_DIR / f"batch_{i:04d}.json"

        if cache_file.exists() and meta_file.exists():
            detected = pd.read_parquet(cache_file)["detected"]
            n_parsed = json.load(open(meta_file))["n_parsed"]
            print(f"  batch {i + 1}/{len(batches)}: cached ({n_parsed} samples)")
        else:
            print(f"  batch {i + 1}/{len(batches)}: downloading {len(batch_ids)} files...")
            max_batch_attempts = 4
            for attempt in range(max_batch_attempts):
                try:
                    tar_bytes = download_batch(batch_ids)
                    detected, n_parsed = parse_batch(tar_bytes)
                    break
                except (tarfile.TarError, EOFError) as e:
                    if attempt == max_batch_attempts - 1:
                        raise
                    wait = 5 * (2 ** attempt)
                    print(f"    batch response corrupted/truncated ({e}); re-downloading in {wait}s "
                          f"({attempt + 1}/{max_batch_attempts})")
                    import time
                    time.sleep(wait)
            detected.to_frame("detected").to_parquet(cache_file)
            json.dump({"n_parsed": n_parsed}, open(meta_file, "w"))
            print(f"    parsed {n_parsed} samples, {len(detected)} genes")

        if total_detected is None:
            total_genes_seen = set(detected.index)
            total_detected = detected.reindex(sorted(total_genes_seen), fill_value=0)
        else:
            total_genes_seen |= set(detected.index)
            total_detected = total_detected.reindex(sorted(total_genes_seen), fill_value=0).add(
                detected.reindex(sorted(total_genes_seen), fill_value=0), fill_value=0
            )
        total_parsed += n_parsed

    print(f"\nParsed {total_parsed} samples total, {len(total_genes_seen)} protein-coding genes seen.")
    detection_frac = total_detected / total_parsed
    background_genes = detection_frac[detection_frac >= DETECTION_FRACTION].index.tolist()
    print(f"Genes with CPM>={CPM_THRESHOLD} in >={DETECTION_FRACTION:.0%} of {total_parsed} samples: "
          f"{len(background_genes)}")

    header = [
        "# transcriptome_background_genes.txt",
        "#",
        "# WHAT: gene-symbol list used as the statistical background for GO/KEGG enrichment",
        "#       (gp.enrich(..., background=...)) in TCGA_BRCA_histological_type_pipeline.ipynb",
        "#       and TCGA_BRCA_ML_pipeline.ipynb, in place of Enrichr's default genome-wide",
        "#       background -- see 'Gene Ontology & KEGG Enrichment: Methodology' in either",
        "#       notebook for why a genome-wide background is inappropriate for an RNA-seq-",
        "#       derived gene list.",
        "#",
        "# HOW GENERATED: data/build_transcriptome_background.py",
        f"#   - Source: ALL {total_parsed} TCGA-BRCA Primary Tumor RNA-seq samples "
        "(STAR - Counts workflow) available in GDC (api.gdc.cancer.gov) at generation time --",
        "#     not a subsample, so this is deterministic/reproducible up to GDC's own data updates.",
        f"#   - Filter: protein-coding genes with CPM >= {CPM_THRESHOLD:.1f} in "
        f">= {DETECTION_FRACTION:.0%} of the {total_parsed} tumors",
        f"#   - Result: {len(background_genes)} of {len(total_genes_seen)} protein-coding genes seen "
        "passed the filter",
        f"#   - Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "#   - Re-run `python3 data/build_transcriptome_background.py` to regenerate "
        "(re-queries GDC live; downloads are cached per-batch under data/.background_cache/ "
        "so an interrupted run resumes instead of restarting)",
        "#",
        "# FORMAT: one gene symbol per line below this header, sorted alphabetically. "
        "Lines starting with '#' are header comments -- notebook loaders skip them.",
        "",
    ]

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(header))
        f.write("\n".join(sorted(background_genes)))
    print(f"Saved to {OUT_PATH}")


if __name__ == "__main__":
    main()
