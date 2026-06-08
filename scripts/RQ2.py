#!/usr/bin/env python
# coding: utf-8

import argparse
from pathlib import Path
import utils

RQ = "RQ2"
ROOT = Path.cwd()


def parse_args():
    parser = argparse.ArgumentParser(description="Compute RQ2 GO-term stability scores.")
    parser.add_argument("--day_dir", type=str, required=True, help="Day directory produced by the pipeline, e.g. 2026_06_03")
    return parser.parse_args()


def save_scores(output_dir: Path, dataset_id: str, scores):
    output_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(output_dir / f"{RQ}_{dataset_id}_q5_gontosim_scores.csv", index=False)


def save_summary_tables(output_dir: Path, dataset_id: str, summary):
    if summary.empty:
        return
    summary.to_csv(output_dir / f"{RQ}_{dataset_id}_gontosim_summary.csv", index=False)
    for approach in utils.APPROACH_SPECS:
        for go_category in utils.GO_CATEGORIES:
            subset = summary[(summary["approach"] == approach) & (summary["go_category"] == go_category)]
            if not subset.empty:
                subset.to_csv(output_dir / f"{RQ}_{dataset_id}_{approach}_{go_category}_summary.csv", index=False)


def process_dataset(dataset_dir: Path, day_dir: str):
    dataset_id = dataset_dir.name
    output_dir = ROOT / "ImgOut" / day_dir / RQ / dataset_id / "INPUT"
    run_dirs = utils.discover_run_dirs(dataset_dir, baseline_only=False)
    print(f"DATASET: {dataset_id}")
    print(f"Run directories: {len(run_dirs)}")
    scores = utils.compute_rq2_scores(dataset_id, run_dirs)
    summary = utils.summarize_rq2_scores(scores)
    save_scores(output_dir, dataset_id, scores)
    save_summary_tables(output_dir, dataset_id, summary)
    print(f"Saved {len(scores)} scores for dataset {dataset_id}")


def main():
    args = parse_args()
    day_dir = args.day_dir
    day_base = ROOT / day_dir
    dataset_dirs = utils.discover_dataset_dirs(day_base, baseline_only=False)
    if not dataset_dirs:
        print(f"No dataset folders with FEA results found in: {day_base}")
        return
    print("Discovered datasets:")
    for dataset_dir in dataset_dirs:
        print(f"  - {dataset_dir.name}")
    for dataset_dir in dataset_dirs:
        process_dataset(dataset_dir, day_dir)
    print("\nRQ2 GO-term stability analysis completed.")


if __name__ == "__main__":
    main()