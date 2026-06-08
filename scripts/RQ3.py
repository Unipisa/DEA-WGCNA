#!/usr/bin/env python
# coding: utf-8

import argparse
from pathlib import Path
import utils


RQ = "RQ3"
ROOT = Path.cwd()


def parse_args():
    parser = argparse.ArgumentParser(description="Compute RQ3 DEA/WGCNA GO-term coherence tables.")
    parser.add_argument(
        "--day_dir",
        type=str,
        required=True,
        help="Day directory produced by the pipeline, e.g. 2026_06_03",
    )
    return parser.parse_args()


def save_similarity_tables(output_dir: Path, dataset_id: str, go_category: str, mean_matrix, count_matrix):
    output_dir.mkdir(parents=True, exist_ok=True)
    mean_matrix.to_csv(output_dir / f"{RQ}_{dataset_id}_{go_category}_cross_mean_similarity.csv")
    count_matrix.to_csv(output_dir / f"{RQ}_{dataset_id}_{go_category}_cross_informative_run_counts.csv")


def process_dataset(dataset_dir: Path, day_dir: str):
    dataset_id = dataset_dir.name
    output_dir = ROOT / "ImgOut" / day_dir / RQ / dataset_id / "INPUT"
    run_dirs = utils.discover_run_dirs(dataset_dir, baseline_only=True)
    print(f"DATASET: {dataset_id}")
    print(f"Run directories: {len(run_dirs)}")
    for go_category in utils.GO_CATEGORIES:
        print(f"GO category: {go_category}")
        terms_dea = utils.collect_terms_by_run_and_param(run_dirs=run_dirs, approach="DEA", go_category=go_category)
        terms_wgcna = utils.collect_terms_by_run_and_param(run_dirs=run_dirs, approach="WGCNA", go_category=go_category)
        dea_params = sorted({param for param_map in terms_dea.values() for param in param_map.keys()}, key=utils.APPROACH_SPECS["DEA"]["sort_key"],)
        wgcna_params = sorted({param for param_map in terms_wgcna.values() for param in param_map.keys()}, key=utils.APPROACH_SPECS["WGCNA"]["sort_key"],)
        if not dea_params or not wgcna_params:
            print("  No DEA/WGCNA parameters found; skipping.")
            continue
        mean_matrix, count_matrix = utils.compute_cross_approach_similarity_matrix(terms_dea, terms_wgcna, dea_params, wgcna_params)
        save_similarity_tables(output_dir, dataset_id, go_category, mean_matrix, count_matrix)


def main():
    args = parse_args()
    day_dir = args.day_dir
    day_base = ROOT / day_dir
    dataset_dirs = utils.discover_dataset_dirs(day_base, baseline_only=True)
    if not dataset_dirs:
        print(f"No dataset folders with FEA results found in: {day_base}")
        return
    print("Discovered datasets:")
    for dataset_dir in dataset_dirs:
        print(f"  - {dataset_dir.name}")
    for dataset_dir in dataset_dirs:
        process_dataset(dataset_dir, day_dir)
    print("\nRQ3 DEA/WGCNA coherence table generation completed.")


if __name__ == "__main__":
    main()