import re
import os
import numpy as np
import pandas as pd


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
        "param_name": r"$FC^{threshold}$",
        "param":['0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.1', '1.2', '1.3', '1.4', '1.5'],
        "param_int":[50,60,70,80,90,100,110,120,130,140,150],
    },
    "WGCNA": {
        "pattern": re.compile(r"^wgcna_go_(\d+)_(\d+)\.csv$", re.IGNORECASE),
        "sort_key": lambda label: tuple(int(x) for x in label.split("_")),
        "param_name":  r"$WGCNA^{threshold}$",
        "param": ['0.60', '0.65', '0.70', '0.75', '0.80', '0.85', '0.90', '0.95'],
        "param_int":[60,65,70,75,80,85,90,95],    
    },
}


### GENE SETS UTILS ### 

def dea_thr(fname):  # deg5_100.txt -> 100
    return int(re.findall(r"\d+", fname)[1])

def wgcna_thr(fname):  # wgcna_genes_60_60_.csv -> (60,60)
    nums = re.findall(r"\d+", fname)
    return int(nums[0]), int(nums[1])

def read_files(folder):
    dea_files = [f for f in os.listdir(folder) if re.match(r"deg5_\d+\.txt$", f)]
    wgcna_files = [f for f in os.listdir(folder) if re.match(r"wgcna_genes_\d+_\d+.*\.csv$", f)]
    dea_files = sorted(dea_files, key=dea_thr)
    wgcna_files = sorted(wgcna_files, key=wgcna_thr)
    return dea_files, wgcna_files

def read_gene_set(path: str) -> set:
    ext = os.path.splitext(path)[1].lower()

    if ext == ".txt":
        try:
            s = pd.read_csv(path, header=None, dtype=str).iloc[:, 0].dropna().astype(str)
        except:
            return set()
        return set(s.tolist())

    if ext == ".csv":
        df = pd.read_csv(path, dtype=str)
        col = df.columns[0]

        s = df[col].dropna().astype(str)
        return set(s.tolist())
    raise ValueError(f"Unsupported extension: {ext}")

def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0 
    return len(a & b) / len(a | b)


def compute_intra_jaccard(sets_list):
    n = len(sets_list)
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            val = jaccard(sets_list[i], sets_list[j])
            M[i, j] = val
            M[j, i] = val
    return M

def compute_cross_jaccard(w_sets, d_sets):
    M = np.zeros((len(w_sets), len(d_sets)), dtype=float)
    for i, w in enumerate(w_sets):
        for j, d in enumerate(d_sets):
            M[i, j] = jaccard(w, d)
    return M
