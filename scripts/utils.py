import pandas as pd
import numpy as np
import re
from pathlib import Path
from utils_similarity import Similarity_of_Set_of_GOTerms

### CONFIGURATION ### 

FEA_DIR = 'FEA'
GO_CATEGORY_MAP = {
    "GOTERM_BP_DIRECT": "BP",
    "GOTERM_CC_DIRECT": "CC",
    "GOTERM_MF_DIRECT": "MF",
}

GO_CATEGORIES = ["BP", "CC", "MF"]

APPROACH_SPECS = {
    "DEA": {
        "pattern": re.compile(r"^dea_go_(\d+)\.csv$", re.IGNORECASE),
        "sort_key": lambda label: int(label),
        "display_name": "FC",
        "param_name": "FC_TH",
        "label": lambda match: str(int(match.group(1)) / 100),

    },
    "WGCNA": {
        "pattern": re.compile(r"^wgcna_go_(\d+)_(\d+)\.csv$", re.IGNORECASE),
        "label": lambda match: f"{int(match.group(1)) / 100:.2f}",
        "sort_key": lambda label: tuple(int(x) for x in label.split("_")),
        "display_name": "GSMM",
        "param_name": "GSMM_TH",

    },
}


### PATH UTILS ### 

def discover_run_dirs(dataset_dir: Path, baseline_only: bool = False) -> list[Path]:
    run_dirs = []
    if not dataset_dir.exists():
        return run_dirs
    for item in sorted(dataset_dir.iterdir()):
        if not item.is_dir():
            continue
        name = item.name.lower()
        if "removed" not in name or "_k" not in name:
            continue
        if baseline_only and "_k0_" not in name:
            continue
        if (item / FEA_DIR).exists():
            run_dirs.append(item)
    return run_dirs


def discover_dataset_dirs(day_base: Path, baseline_only: bool = False) -> list[Path]:
    if not day_base.exists():
        raise FileNotFoundError(f"Day directory not found: {day_base}")
    dataset_dirs = []
    for item in sorted(day_base.iterdir()):
        if not item.is_dir():
            continue
        if discover_run_dirs(item, baseline_only=baseline_only):
            dataset_dirs.append(item)
    return dataset_dirs

### ENRICHMENT UTILS ### 

def read_enrichment_file(csv_path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(csv_path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["Category", "Term"])
    except Exception:
        return pd.DataFrame(columns=["Category", "Term"])
    if "Term" in df.columns:
        return df
    try:
        df = pd.read_csv(csv_path, index_col=0)
    except Exception:
        return pd.DataFrame(columns=["Category", "Term"])
    return df if "Term" in df.columns else pd.DataFrame(columns=["Category", "Term"])


def extract_go_ids(csv_path: Path, go_category: str):
    df = read_enrichment_file(csv_path)
    if df.empty or "Category" not in df.columns or "Term" not in df.columns:
        return set()
    subset = df[df["Category"].map(GO_CATEGORY_MAP).eq(go_category)].copy()
    if subset.empty:
        return set()
    go_ids = (subset["Term"].dropna().astype(str).map(lambda x: x.split("~")[0].strip()))
    return {go_id for go_id in go_ids if go_id.startswith("GO:")}


def get_param_label(approach: str, match) -> str:
    if approach == "DEA":
        return match.group(1)
    if approach == "WGCNA":
        return f"{match.group(1)}_{match.group(2)}"
    raise ValueError(f"Unknown approach: {approach}")


def collect_terms_by_run_and_param(run_dirs: list[Path], approach: str, go_category: str,) -> dict[str, dict[str, set[str]]]:
    spec = APPROACH_SPECS[approach]
    terms_by_run = {}

    for run_dir in run_dirs:
        fea_dir = run_dir / FEA_DIR
        param_sets = {}

        for csv_path in sorted(fea_dir.glob("*.csv")):
            match = spec["pattern"].match(csv_path.name)

            if not match:
                continue

            param_label = get_param_label(approach, match)
            param_sets[param_label] = extract_go_ids(csv_path, go_category)

        if param_sets:
            terms_by_run[run_dir.name] = param_sets

    return terms_by_run

### SIMILARITY UTILS ### 


def compute_intra_parameter_similarity_matrix(terms_by_run: dict, ordered_params: list[str], method: str = "GOntoSim") -> tuple[pd.DataFrame, pd.DataFrame]:
    mean_matrix = pd.DataFrame(np.nan, index=ordered_params, columns=ordered_params, dtype=float,)
    count_matrix = pd.DataFrame(0, index=ordered_params, columns=ordered_params, dtype=int)

    for i, param_i in enumerate(ordered_params):
        for j, param_j in enumerate(ordered_params):
            if i < j:
                continue
            similarities = []
            for param_map in terms_by_run.values():
                set_i = param_map.get(param_i, set())
                set_j = param_map.get(param_j, set())
                if not set_i or not set_j:
                    continue
                score = Similarity_of_Set_of_GOTerms(sorted(set_i), sorted(set_j), method)
                if pd.notna(score):
                    similarities.append(float(score))
            if similarities:
                mean_matrix.loc[param_i, param_j] = float(np.mean(similarities))
                count_matrix.loc[param_i, param_j] = len(similarities)
    return mean_matrix, count_matrix

def compute_cross_approach_similarity_matrix(terms_by_run_a: dict, terms_by_run_b: dict, ordered_params_a: list[str], ordered_params_b: list[str], method: str = "GOntoSim") -> tuple[pd.DataFrame, pd.DataFrame]:
    mean_matrix = pd.DataFrame(np.nan, index=ordered_params_a, columns=ordered_params_b, dtype=float)
    count_matrix = pd.DataFrame(0, index=ordered_params_a, columns=ordered_params_b, dtype=int)

    for param_i in ordered_params_a:
        for param_j in ordered_params_b:
            similarities = []
            for param_map_a in terms_by_run_a.values():
                set_i = param_map_a.get(param_i, set())
                for param_map_b in terms_by_run_b.values():
                    set_j = param_map_b.get(param_j, set())
                    if not set_i or not set_j:
                        continue
                    score = Similarity_of_Set_of_GOTerms(sorted(set_i), sorted(set_j), method,)
                    if pd.notna(score):
                        similarities.append(float(score))
            if similarities:
                mean_matrix.loc[param_i, param_j] = float(np.mean(similarities))
                count_matrix.loc[param_i, param_j] = len(similarities)
    return mean_matrix, count_matrix

def build_lower_triangle_annotation(mean_matrix: pd.DataFrame) -> pd.DataFrame:
    annotation = pd.DataFrame("", index=mean_matrix.index, columns=mean_matrix.columns)
    for row_label in mean_matrix.index:
        for col_label in mean_matrix.columns:
            value = mean_matrix.loc[row_label, col_label]
            annotation.loc[row_label, col_label] = "NA" if pd.isna(value) else f"{value:.3f}"
    return annotation


def extract_k(run_dir: Path) -> int | None:
    match = re.search(r"_k(\d+)_", run_dir.name)
    return int(match.group(1)) if match else None


def compute_rq2_scores(dataset_id: str, run_dirs: list[Path]) -> pd.DataFrame:
    rows = []
    for approach in APPROACH_SPECS:
        for go_category in GO_CATEGORIES:
            terms_by_run = collect_terms_by_run_and_param(run_dirs, approach, go_category)
            ref_run = next((run for run in terms_by_run if extract_k(Path(run)) == 0), None)
            if ref_run is None:
                continue
            ref_terms_by_param = terms_by_run[ref_run]
            for run_name, current_terms_by_param in terms_by_run.items():
                k_value = extract_k(Path(run_name))
                if k_value is None or k_value == 0:
                    continue
                for param_label, current_terms in current_terms_by_param.items():
                    ref_terms = ref_terms_by_param.get(param_label, set())
                    if not current_terms or not ref_terms:
                        continue
                    score = Similarity_of_Set_of_GOTerms(sorted(current_terms), sorted(ref_terms), "GOntoSim")
                    if pd.isna(score):
                        continue
                    rows.append({"dataset": dataset_id, "approach": approach, "go_category": go_category, "k": k_value, "run": run_name, "param": param_label, "similarity": float(score), "n_terms_current": len(current_terms), "n_terms_reference": len(ref_terms)})
    return pd.DataFrame(rows)


def summarize_rq2_scores(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    summary = scores.groupby(["dataset", "approach", "go_category", "k", "param"], as_index=False).agg(mean_similarity=("similarity", "mean"), std_similarity=("similarity", "std"), n_runs=("similarity", "count"), mean_terms_current=("n_terms_current", "mean"), mean_terms_reference=("n_terms_reference", "mean"))
    summary["std_similarity"] = summary["std_similarity"].fillna(0)
    return summary