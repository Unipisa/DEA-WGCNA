#!/usr/bin/env python3

import argparse
import os
import re
import time
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from suds.client import Client
from goatools.obo_parser import GODag
from goatools.godag_plot import plot_gos

DAVID_WSDL = "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl"
DAVID_ENDPOINT = (
    "https://davidbioinformatics.nih.gov/webservice/services/"
    "DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/"
)
EASE_THRESHOLD = 0.1
MIN_GENE_COUNT = 2
BENJAMINI_THRESHOLD = 0.05
EMPTY_RESULT_COLUMNS = [
    "Category",
    "Term",
    "Count",
    "%",
    "Pvalue",
    "Genes",
    "List Total",
    "Pop Hits",
    "Pop Total",
    "Fold Enrichment",
    "Bonferroni",
    "Benjamini",
    "FDR",
]

RUN_FOLDER_REGEX = re.compile(r"(?:^|_)k\d+(?:_|$)", re.IGNORECASE)
DEA_FILE_REGEX = re.compile(r"^deg5_(\d+)\.txt$", re.IGNORECASE)
WGCNA_FILE_REGEX = re.compile(r"^wgcna_genes_(\d+)_(\d+)_\.csv$", re.IGNORECASE)

def empty_result_frame():
    return pd.DataFrame(columns=EMPTY_RESULT_COLUMNS)

def write_empty_results(out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    empty_result_frame().to_csv(out_path)

def load_gene_list(input_path):
    if not os.path.exists(input_path):
        return []
    try:
        if input_path.endswith(".csv"):
            df = pd.read_csv(input_path)
            gene_col = None
            for col in ["x", "0", "gene", "Gene", "ENTREZID", "entrez", "Entrez"]:
                if col in df.columns:
                    gene_col = col
                    break
            if gene_col is None and len(df.columns) > 0:
                gene_col = df.columns[0]
            genes = df[gene_col].dropna().astype(str).tolist() if gene_col is not None else []
        else:
            df = pd.read_csv(input_path, header=None, names=["gene"])
            genes = df["gene"].dropna().astype(str).tolist()
    except pd.errors.EmptyDataError:
        return []
    cleaned = []
    for gene in genes:
        value = gene.strip()
        if not value:
            continue
        if value.lower() in {"x", "gene", "entrez", "entrezid"}:
            continue
        cleaned.append(value)
    return cleaned

def fetch_david_results(genes, method, run_id, david_email):
    client = Client(DAVID_WSDL)
    client.wsdl.services[0].setlocation(DAVID_ENDPOINT)
    client.service.authenticate(david_email)

    input_ids = ",".join(genes)
    list_name = f"list_{method}_{run_id}"
    client.service.addList(input_ids, "ENTREZ_GENE_ID", list_name, 0)
    client.service.setCategories("GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT")
    chart_report = client.service.getChartReport(EASE_THRESHOLD, MIN_GENE_COUNT)

    if not chart_report:
        return empty_result_frame()

    records = []
    for record in chart_report:
        records.append({
            "Category": str(record.categoryName),
            "Term": str(record.termName),
            "Count": str(record.listHits),
            "%": str(record.percent),
            "Pvalue": str(record.ease),
            "Genes": str(record.geneIds),
            "List Total": str(record.listTotals),
            "Pop Hits": str(record.popHits),
            "Pop Total": str(record.popTotals),
            "Fold Enrichment": str(record.foldEnrichment),
            "Bonferroni": str(record.bonferroni),
            "Benjamini": str(record.benjamini),
            "FDR": str(record.afdr),
        })

    df_result = pd.DataFrame(records)
    if df_result.empty:
        return empty_result_frame()
    return df_result[df_result["Benjamini"].astype(float) < BENJAMINI_THRESHOLD].reset_index(drop=True)


def find_overlap(folder, output_file, method):
    fname = os.path.join(folder, output_file)
    if not os.path.exists(fname):
        return

    df = pd.read_csv(fname, index_col=0)
    if df.empty:
        return

    if len(df) > 20:
        df = df.sort_values("Benjamini").head(20).reset_index(drop=True)

    df["Genes"] = df["Genes"].astype(str).fillna("")
    df["GeneList"] = df["Genes"].apply(lambda x: set(x.split(",")) if x else set())
    df["NumGenes"] = df["GeneList"].apply(len)

    if len(df) < 2:
        return

    intersection_matrix = pd.DataFrame(0, index=df["Term"], columns=df["Term"])
    for i in range(len(df)):
        for j in range(i, len(df)):
            intersection_count = len(df["GeneList"].iloc[i].intersection(df["GeneList"].iloc[j]))
            intersection_matrix.iloc[i, j] = intersection_count
            intersection_matrix.iloc[j, i] = intersection_count

    intersection_matrix.index = [
        f"{str(term).split('~')[0]} (n={num})"
        for term, num in zip(df["Term"], df["NumGenes"])
    ]
    intersection_matrix.columns = intersection_matrix.index
    mask = np.triu(np.ones(intersection_matrix.shape), k=1)

    plt.figure(figsize=(max(8, len(df) * 0.6), max(6, len(df) * 0.5)))
    sns.heatmap(intersection_matrix, mask=mask, cmap="coolwarm", annot=True, fmt="g", cbar_kws={"label": "Number of Common Genes"})
    plt.title(f"Gene Overlap Between Enriched Terms ({method.upper()})")
    plt.tight_layout()
    plt.savefig(os.path.join(folder, f"pathway_intersection_{method}.pdf"), format="pdf", bbox_inches="tight")
    plt.close()


def plot_subtree(folder, output_file, method):
    fname = os.path.join(folder, output_file)
    if not os.path.exists(fname):
        return
    obo_file = "go-basic.obo"
    df = pd.read_csv(fname, index_col=0)
    if df.empty:
        return

    go_dag = GODag(obo_file, optional_attrs={"relationship"})
    for cat, df_cat in df.groupby("Category"):
        goids = [str(term).split("~")[0] for term in df_cat["Term"].values if "~" in str(term)]
        if not goids:
            continue
        output_path = os.path.join(folder, f"{method}_dag_plot_{cat.replace('GOTERM_', '').replace('_DIRECT', '')}.png")
        try:
            plot_gos(fout_png=output_path, goids=goids, obo_dag=go_dag, title="")
        except Exception:
            print(f"Error in plotting DAG for category {cat} in {output_file}")
            continue


def find_dataset_dirs(base_dir):
    base_dir = os.path.abspath(base_dir)
    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"Base directory not found: {base_dir}")
    dataset_dirs = []
    for item in sorted(os.listdir(base_dir)):
        item_path = os.path.join(base_dir, item)
        if not os.path.isdir(item_path):
            continue
        if find_run_folders(item_path):
            dataset_dirs.append(item_path)
    return dataset_dirs


def find_run_folders(dataset_dir):
    folders = []
    for item in sorted(os.listdir(dataset_dir)):
        item_path = os.path.join(dataset_dir, item)
        if not os.path.isdir(item_path):
            continue
        if RUN_FOLDER_REGEX.search(item) and "removed" in item.lower():
            folders.append(item_path)
    return folders


def collect_parameterized_inputs(run_dir):
    inputs = []
    for file_name in sorted(os.listdir(run_dir)):
        dea_match = DEA_FILE_REGEX.match(file_name)
        if dea_match:
            threshold = dea_match.group(1)
            inputs.append({
                "method": "dea",
                "label": threshold,
                "input_file": file_name,
                "output_file": f"dea_go_{threshold}.csv",
                "plot_stem": f"dea_{threshold}",
            })
            continue

        wgcna_match = WGCNA_FILE_REGEX.match(file_name)
        if wgcna_match:
            power = wgcna_match.group(1)
            min_module = wgcna_match.group(2)
            label = f"{power}_{min_module}"
            inputs.append({
                "method": "wgcna",
                "label": label,
                "input_file": file_name,
                "output_file": f"wgcna_go_{label}.csv",
                "plot_stem": f"wgcna_{label}",
            })
    return inputs


def ensure_output_dir(run_dir, output_subdir):
    output_dir = os.path.join(run_dir, output_subdir)
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def run_single_fea(run_dir, output_dir, entry, david_email, run_id):
    input_path = os.path.join(run_dir, entry["input_file"])
    output_path = os.path.join(output_dir, entry["output_file"])

    print(f"\n{'=' * 60}")
    print(f"Method: {entry['method'].upper()} | Parameter: {entry['label']}")
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    print(f"{'=' * 60}")

    genes = load_gene_list(input_path)
    print(f"Number of genes: {len(genes)}")
    log_file = os.path.join(output_dir, "gene_counts.txt")
    
    with open(log_file, "a") as f:
        print(f"{run_id},{entry['method']},{entry['label']},{len(genes)}", file=f)

    if not genes:
        print("No genes found after parsing input file.")
        write_empty_results(output_path)
        return False, output_path

    try:
        df_result = fetch_david_results(genes=genes, method=entry["method"], run_id=run_id, david_email=david_email)

    except Exception as exc:
        error_log = os.path.join(output_dir, "david_errors.log")
        with open(error_log, "a") as f:
            print(f"{run_id},{entry['method']},{entry['label']},{exc}", file=f)
        print(f"DAVID API error: {exc}")
        return False, output_path

    if df_result.empty:
        print(f"No enrichment results passed the Benjamini threshold (< {BENJAMINI_THRESHOLD}).")
        write_empty_results(output_path)
        return False, output_path

    df_result.to_csv(output_path)
    print(f"Results saved: {output_path}")
    return True, output_path


def create_plots(output_dir, output_file, plot_stem, method):
    local_output_file = os.path.basename(output_file)
    find_overlap(output_dir, local_output_file, plot_stem)
    plot_subtree(output_dir, local_output_file, plot_stem)



def main():
    parser = argparse.ArgumentParser(
        description="Run DAVID FEA on robustness results for each dataset and parameter."
    )
    parser.add_argument(
        "--david_email",
        type=str,
        required=True,
        help="Registered DAVID account email address",
    )
    parser.add_argument(
        "--base_dir",
        type=str,
        required=True,
        help="Base output directory containing dataset folders, e.g. 2026_06_03",
    )
    parser.add_argument(
        "--output_subdir",
        type=str,
        default="FEA",
        help="Subdirectory created inside each run folder to store results",
    )
    parser.add_argument(
        "--skip_plots",
        action="store_true",
        help="Skip overlap heatmaps and DAG plots",
    )
    args = parser.parse_args()

    dataset_dirs = find_dataset_dirs(args.base_dir)
    if not dataset_dirs:
        print(f"No dataset directories with run folders found under {os.path.abspath(args.base_dir)}")
        return

    print("Resolved dataset directories:")
    for dataset_dir in dataset_dirs:
        print(f"  - {dataset_dir}")

    global_run_id = 0
    for dataset_dir in dataset_dirs:
        run_folders = find_run_folders(dataset_dir)
        if not run_folders:
            print(f"\nNo run folders found in dataset directory: {dataset_dir}")
            continue

        print(f"\n{'#' * 70}")
        print(f"DATASET: {os.path.basename(dataset_dir)}")
        print(f"Found {len(run_folders)} run folder(s)")
        print(f"{'#' * 70}")

        for run_dir in run_folders:
            print(f"\nProcessing run folder: {run_dir}")
            entries = collect_parameterized_inputs(run_dir)
            if not entries:
                print("No DEA/WGCNA parameterized input files found.")
                continue

            output_dir = ensure_output_dir(run_dir, args.output_subdir)
            print(f"Discovered {len(entries)} input files. Results will be saved in: {output_dir}")

            for entry in entries:
                success, output_path = run_single_fea(
                    run_dir=run_dir,
                    output_dir=output_dir,
                    entry=entry,
                    david_email=args.david_email,
                    run_id=global_run_id,
                )
                if success and not args.skip_plots:
                    create_plots(
                        output_dir=output_dir,
                        output_file=output_path,
                        plot_stem=entry["plot_stem"],
                        method=entry["method"],
                    )
                time.sleep(1)
                global_run_id += 1

    print(f"\n{'=' * 60}")
    print("INTRA-DATASET ENRICHMENT ANALYSIS COMPLETE!")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
