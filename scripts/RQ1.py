#!/usr/bin/env python
# coding: utf-8

import argparse
from pathlib import Path
import utils

RQ = "RQ1"
ROOT = Path.cwd()

def parse_args():
    parser = argparse.ArgumentParser(description="Compute RQ1 intra-approach GO-term similarity tables.")
    parser.add_argument(
        "--day_dir",
        type=str,
        required=True,
        help="Day directory produced by the pipeline, e.g. 2026_06_03",
    )
    return parser.parse_args()


def save_similarity_tables(output_dir: Path, dataset_id: str, approach: str, go_category: str, mean_matrix, count_matrix, annotation):
    output_dir.mkdir(parents=True, exist_ok=True)
    mean_matrix.to_csv(output_dir / f"{RQ}_{dataset_id}_{approach}_{go_category}_mean_similarity.csv")
    count_matrix.to_csv(output_dir / f"{RQ}_{dataset_id}_{approach}_{go_category}_informative_run_counts.csv")
    annotation.to_csv(output_dir / f"{RQ}_{dataset_id}_{approach}_{go_category}_lower_triangle_annotation.csv")


def process_dataset(dataset_dir: Path, day_dir: str):
    dataset_id = dataset_dir.name
    output_dir = ROOT / "ImgOut" / day_dir / RQ / dataset_id / "INPUT"
    run_dirs = utils.discover_run_dirs(dataset_dir, baseline_only=True)
    print(f"DATASET: {dataset_id}")
    print(f"Run directories: {len(run_dirs)}")

    for approach, spec in utils.APPROACH_SPECS.items():
        print(f"\nApproach: {approach}")
        for go_category in utils.GO_CATEGORIES:
            print(f"GO category: {go_category}")
            terms_by_run = utils.collect_terms_by_run_and_param(run_dirs=run_dirs, approach=approach, go_category=go_category)
            ordered_params = sorted({param for param_map in terms_by_run.values() for param in param_map.keys()}, key=spec["sort_key"],)
            if not ordered_params:
                print("  No parameters found; skipping.")
                continue

            mean_matrix, count_matrix = utils.compute_intra_parameter_similarity_matrix(terms_by_run=terms_by_run, ordered_params=ordered_params)
            annotation = utils.build_lower_triangle_annotation(mean_matrix)
            save_similarity_tables(output_dir, dataset_id, approach, go_category, mean_matrix, count_matrix, annotation)


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
    print("\nRQ1 similarity table generation completed.")


if __name__ == "__main__":
    main()