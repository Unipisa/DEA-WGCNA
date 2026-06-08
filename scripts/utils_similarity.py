# Modified version of GOntoSim Similarity.py https://github.com/cbrl-nuces/GOntoSim/blob/main/Similarity.py


from goatools.semantic import common_parent_go_ids, min_branch_length, semantic_distance, semantic_similarity
from goatools.obo_parser import GODag
from goatools.base import get_godag
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim
from goatools.semantic import TermCounts, get_info_content
from goatools.associations import read_associations
from goatools.semantic import deepest_common_ancestor

import numpy as np
import pandas as pd
from collections import namedtuple
from collections import defaultdict
import os
import tempfile

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn import metrics
from scipy.cluster.hierarchy import dendrogram 
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import squareform

import sys
import random

Gene = namedtuple('Gene', 'GeneName GOAnnotations')
_term_similarity_cache = {}
_set_similarity_cache = {}

def _resolve_go_basic_obo():
	candidates = [
		os.environ.get("GO_BASIC_OBO"),
		os.path.join(os.path.dirname(os.path.abspath(__file__)), "go-basic.obo"),
		os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "go-basic.obo"),
		os.path.join(os.getcwd(), "go-basic.obo"),
		"go-basic.obo",
	]
	for candidate in candidates:
		if candidate and os.path.exists(candidate):
			return candidate
	return "go-basic.obo"

def resolve_obo_file(folder):
	candidates = [
		os.path.join(folder, "go-basic.obo"),
		os.path.join(os.path.dirname(__file__), "go-basic.obo"),
		os.path.join(os.path.dirname(folder), "go-basic.obo"),
		_resolve_go_basic_obo(),
	]
	for candidate in candidates:
		if candidate and os.path.exists(candidate):
			return candidate
	return None

def _normalize_similarity_method(method):
	if method is None:
		return "GOntoSim"
	method_str = str(method).strip()
	method_key = method_str.lower()
	alias_map = {
		"bp": "GOntoSim",
		"cc": "GOntoSim",
		"mf": "GOntoSim",
		"biological_process": "GOntoSim",
		"cellular_component": "GOntoSim",
		"molecular_function": "GOntoSim",
		"gontosim": "GOntoSim",
		"goontosim": "GOntoSim",
		"wang": "wang",
		"baseline": "Baseline",
		"baseline_lca": "Baseline_LCA",
		"baseline_desc": "Baseline_Desc",
		"baseline_desc_only": "Baseline_Desc_only",
		"gogo": "GOGO",
	}
	return alias_map.get(method_key, method_str)

def _sanitize_go_terms(go_terms):
	clean_terms = []
	seen = set()
	for go_term in go_terms:
		if pd.isna(go_term):
			continue
		go_id = str(go_term).strip()
		if not go_id.startswith("GO:"):
			continue
		if go_id not in go or go_id in seen:
			continue
		clean_terms.append(go_id)
		seen.add(go_id)
	return clean_terms

def _term_similarity_cache_key(go_id1, go_id2, method):
	return (method, tuple(sorted((go_id1, go_id2))))

def _set_similarity_cache_key(set1, set2, method):
	set_pair = sorted((tuple(set1), tuple(set2)))
	return (method, tuple(set_pair))

go = get_godag(_resolve_go_basic_obo(), optional_attrs={'relationship'})
is_a = 0.8
part_of = 0.6
regulates = 0.7
negatively_regulates = 0.7
positively_regulates = 0.7
reg = 'reg0.7'
c = 0.67

def all_common_parent_go_ids(goids, godag):
	rec = godag[goids[0]]
	candidates = rec.get_all_upper()
	candidates.update({goids[0]})
	for goid in goids[1:]:
		rec = godag[goid]
		parents = rec.get_all_upper()
		parents.update({goid})
		candidates.intersection_update(parents)
	return candidates

def lowest_common_ancestor(goterms, godag):
	return max(all_common_parent_go_ids(goterms, godag), key=lambda t: godag[t].depth)

def all_paths_to_top(term, godag):
		""" Returns all possible paths to the root node"""
		if term not in godag:
			sys.stderr.write("Term %s not found!\n" % term)
			return
		def _all_paths_to_top_recursive(rec):
			if rec.level == 0:
				return [[rec]]
			paths = []
			parents = rec.get_goterms_upper()
			for parent in parents:
				top_paths = _all_paths_to_top_recursive(parent)
				for top_path in top_paths:
					top_path.append(rec)
					paths.append(top_path)
			return paths

		go_term = godag[term]
		return _all_paths_to_top_recursive(go_term)
    
    
def all_paths_to_top_wang(term, godag, optional_relationships):
		""" Returns all possible paths to the root node"""
		if term not in godag:
			sys.stderr.write("Term %s not found!\n" % term)
			return
		def _all_paths_to_top_recursive(rec):
			if rec.level == 0:
				return [[rec]]
			paths = []
			parents = rec.get_goterms_upper_rels(optional_relationships)
			for parent in parents:
				top_paths = _all_paths_to_top_recursive(parent)
				for top_path in top_paths:
					top_path.append(rec)
					paths.append(top_path)
			return paths
		go_term = godag[term]
		return _all_paths_to_top_recursive(go_term)
def Semantic_Value(go_id, go, method):
	if method == 'wang':
		# calculates all paths to top (with all relationships)
		optional_relationships = {'part_of'}
		all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
		S_values = list()
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			path.reverse()
			for idx, term in enumerate(path):
				if idx == 0:
					S_values[index].append((go_id, 1))
				if idx < len(path)-1:
					if term.relationship != {}:
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
							else:
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						else:
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
					else:
						S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
		return final_values(S_values,'max')
	elif method == 'GOGO':
		optional_relationships = {'part_of'}
		all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
		S_values = list()
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			path.reverse()
			for idx, term in enumerate(path):
				if idx == 0:
					S_values[index].append((go_id, 1))
				if idx < len(path)-1:
					if term.relationship != {}:
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
								weight = (1 / (c + len(go[path[idx+1].item_id].children))) + part_of
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight))
							else:
								weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
						else:
							weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
					else:
						weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
						S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
		return final_values(S_values,'max')
	else:
		all_all_paths = all_paths_to_top(go_id, go)
		S_values = list()
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			path.reverse()
			for idx, term in enumerate(path):
				if idx == 0:
					S_values[index].append((go_id, 1))
				if idx < len(path)-1 and path[idx+1].item_id != 'GO:0003674' and path[idx+1].item_id != 'GO:0005575' and path[idx+1].item_id != 'GO:0008150':
					if term.relationship != {}:
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
							else: 
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'regulates' in term.relationship:
							if path[idx+1] in term.relationship['regulates']:
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * regulates))
							else: 
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'negatively_regulates' in term.relationship:
							if path[idx+1] in term.relationship['negatively_regulates']:
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * negatively_regulates))
							else: 
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'positively_regulates' in term.relationship:
							if path[idx+1] in term.relationship['positively_regulates']:
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * positively_regulates))
							else: 
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
						else: 
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					else: 
						S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
				if (term.item_id == 'GO:0003674' or term.item_id == 'GO:0005575' or term.item_id == 'GO:0008150'):
						S_values[index].append((term.item_id, 0))
		if method == 'Baseline':
			return final_values(S_values, 'max')
		if method == 'Baseline_LCA':
			svaluesfinal = final_values(S_values, 'max')
			S_values_Modified = list()
			for index, path in enumerate(S_values):
				S_values_Modified.append([])
				for idx, term1 in enumerate(path):
					ind = [value for (term, value) in svaluesfinal if term == term1[0]]
					S_values_Modified[index].append((path[idx][0],ind[0]))
			SumOfNodesOnEachPath = list()
			for index, path in enumerate(S_values_Modified):
				SumOfNodesOnEachPath.append(sum(x[1] for x in path))
			maxPath = SumOfNodesOnEachPath.index(max(SumOfNodesOnEachPath))
			return S_values_Modified[maxPath], SumOfNodesOnEachPath[maxPath]

def final_values(S_values, isMax):
	unique_terms_s_values = []
	for path in S_values:
		for term in path:
			unique_terms_s_values.append(term)
	unique_terms_s_values = sorted(unique_terms_s_values, key=lambda x: x[0])
	_s_values = {}
	for y, x in unique_terms_s_values: 
		if y in _s_values: 
			_s_values[y].append((y,x)) 
		else: 
			_s_values[y] = [(y, x)]
	final_s_values = []
	if(isMax == 'max'):
		for node in _s_values:
			final_s_values.append(max(_s_values[node]))
	elif (isMax == 'min'):
		for node in _s_values:
			final_s_values.append(min(_s_values[node]))
	return final_s_values
def intersection(lst1, lst2):
	da = {v:k for v, k in lst1}
	db = {v:k for v, k in lst2} 
	return [(da[k],db[k]) for k in da.keys() & db.keys()]
def common_children_go_ids(goids, godag):
	rec = godag[goids[0]]
	candidates = rec.get_all_lower()
	candidates.update({goids[0]})
	for goid in goids[1:]:
		rec = godag[goid]
		children = rec.get_all_lower()
		children.update({goid})
		candidates.intersection_update(children)
	return candidates

def highest_common_descendant(goterms, godag):
	common_children = common_children_go_ids(goterms, godag)
	if len(common_children) != 0:
		# take depth attribute instead of depth to accomodate all relationships
		return min(common_children, key=lambda t: godag[t].depth)
	else:
		return 0
def all_paths_to_bottom(term, godag, x):
	""" Returns all possible paths to the leaf nodes x levels below"""
	if term not in godag:
		sys.stderr.write("Term %s not found!\n" % term)
		return

	def _all_paths_to_bottom_recursive(rec):

		if rec.depth == godag[term].depth+x:
			return [[rec]]
		else:
			paths = []
			children = rec.get_goterms_lower()
			for child in children:
				bottom_paths = _all_paths_to_bottom_recursive(child)
				for bottom_path in bottom_paths:
					bottom_path.append(rec)
					paths.append(bottom_path)
			return paths
	go_term = godag[term]
	return _all_paths_to_bottom_recursive(go_term)
def Downward_Semantic_Value(go_id, go, x):
	all_all_paths = all_paths_to_bottom(go_id, go, x)
	S_values = list()
	for index, path in enumerate(all_all_paths):
		S_values.append([])
		path.reverse()
		for idx, term in enumerate(path):
			if idx == 0:
				S_values[index].append((go_id, 1))
			if idx < len(path)-1:
				if term.relationship != {}:
					if 'part_of' in term.relationship:
						if path[idx+1] in term.relationship['part_of']:
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
						else: 
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'regulates' in term.relationship:
						if path[idx+1] in term.relationship['regulates']:
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * regulates))
						else: 
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'negatively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['negatively_regulates']:
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * negatively_regulates))
						else: 
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'positively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['positively_regulates']:
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * positively_regulates))
						else: 
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
				else: 
					S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
	return final_values(S_values, 'max')

def _plot_go_subgraph(go_id1, go_id2, method, lca=None, hcd=None, sv_a=None, sv_b=None):
	try:
		from goatools.obo_parser import GODag
		from goatools.godag_plot import plot_gos
		import matplotlib.image as mpimg
		import matplotlib.pyplot as plt
	except ImportError:
		return

	folder = os.path.dirname(os.path.abspath(__file__))
	obo_file = resolve_obo_file(folder)
	if obo_file is None:
		return

	goids = []
	for go_id in (go_id1, go_id2, lca, hcd):
		if go_id and go_id in go:
			goids.append(go_id)
	for semantic_values in (sv_a, sv_b):
		if semantic_values:
			goids.extend(term for term, _ in semantic_values if term in go)

	goids = sorted(set(goids))
	if not goids:
		return

	go_dag = GODag(obo_file, optional_attrs={"relationship"})
	with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_file:
		output_path = tmp_file.name

	try:
		plot_gos(fout_png=output_path, goids=goids, obo_dag=go_dag, title="")
		image = mpimg.imread(output_path)
		plt.figure(figsize=(14, 10))
		plt.imshow(image)
		plt.title(f"GO subgraph for {go_id1} vs {go_id2} ({method})")
		plt.axis("off")
		plt.tight_layout()
		plt.show()
	except Exception:
		return
	finally:
		if os.path.exists(output_path):
			os.remove(output_path)

def Similarity_of_Two_GOTerms(go_id1, go_id2, go, method, show_plot=False):
	method = _normalize_similarity_method(method)
	if go_id1 not in go or go_id2 not in go:
		return np.nan
	cache_key = None
	if not show_plot:
		cache_key = _term_similarity_cache_key(go_id1, go_id2, method)
		if cache_key in _term_similarity_cache:
			return _term_similarity_cache[cache_key]
	if method == 'Baseline_LCA':
		if go_id1 == go_id2:
			lca = go_id1
		else:
			lca = lowest_common_ancestor((go_id1, go_id2), go)

		sim_lca = Semantic_Value(lca, go, method)
		sim1 = Semantic_Value(go_id1, go, method)
		sim2 = Semantic_Value(go_id2, go, method)
		sum_sim1_sim2 = (sim1[1] + sim2[1])
		sim = ((sim_lca[1]*2)/sum_sim1_sim2)
		if show_plot:
			_plot_go_subgraph(go_id1, go_id2, method, lca=lca, sv_a=sim1[0], sv_b=sim2[0])
		if cache_key is not None:
			_term_similarity_cache[cache_key] = sim
		return sim
	
	elif method == 'GOntoSim':
		if go_id1 == go_id2:
			hcd = go_id1
		else:
			hcd = highest_common_descendant((go_id1, go_id2), go)
		sv_a = []
		sv_b = []
		if hcd != 0:
			hcd_depth = go[hcd].depth
			go1_depth = go[go_id1].depth
			go2_depth = go[go_id2].depth
			x = hcd_depth - go1_depth
			y = hcd_depth - go2_depth
			sv_a = Downward_Semantic_Value(go_id1, go, x)
			sv_b = Downward_Semantic_Value(go_id2, go, y)
			intersecting_terms = intersection(sv_a,sv_b)
			numerator = sum([x for t in intersecting_terms for x in t])
			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b)
			sim_down = (numerator/denominator)
		else:
			sim_down = 0


		sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline_LCA', show_plot=False)
		if go_id1 == go_id2:
			if show_plot:
				_plot_go_subgraph(go_id1, go_id2, method, hcd=hcd, sv_a=sv_a, sv_b=sv_b)
		sim = (sim_down*0.5) + (sim_upper*0.5)
		if cache_key is not None:
			_term_similarity_cache[cache_key] = sim
		return sim
	elif method == 'wang' or method == 'Baseline' or method == 'GOGO':
		sv_a = Semantic_Value(go_id1, go, method)
		sv_b = Semantic_Value(go_id2, go, method)
		intersecting_terms = intersection(sv_a,sv_b)
		numerator = sum([x for t in intersecting_terms for x in t])
		denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b)
		Similarity = (numerator/denominator)
		if show_plot:
			_plot_go_subgraph(go_id1, go_id2, method, sv_a=sv_a, sv_b=sv_b)
		if cache_key is not None:
			_term_similarity_cache[cache_key] = Similarity
		return Similarity
	elif method == 'Baseline_Desc':
		hcd = highest_common_descendant((go_id1, go_id2), go)
		sv_a = []
		sv_b = []
		if hcd != 0:
			hcd_depth = go[hcd].depth
			go1_depth = go[go_id1].depth
			go2_depth = go[go_id2].depth
			x = hcd_depth - go1_depth
			y = hcd_depth - go2_depth
			sv_a = Downward_Semantic_Value(go_id1, go, x)
			sv_b = Downward_Semantic_Value(go_id2, go, y)
			intersecting_terms = intersection(sv_a,sv_b)
			numerator = sum([x for t in intersecting_terms for x in t])
			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b)
			sim_down = (numerator/denominator)
		else:
			sim_down = 0
		sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline', show_plot=False)
		sim = (sim_down*0.5) + (sim_upper*0.5)
		if show_plot:
			_plot_go_subgraph(go_id1, go_id2, method, hcd=hcd, sv_a=sv_a, sv_b=sv_b)
		if cache_key is not None:
			_term_similarity_cache[cache_key] = sim
		return sim
	elif method == 'Baseline_Desc_only':
		hcd = highest_common_descendant((go_id1, go_id2), go)
		sv_a = []
		sv_b = []
		if hcd != 0:
			hcd_depth = go[hcd].depth
			go1_depth = go[go_id1].depth
			go2_depth = go[go_id2].depth
			x = hcd_depth - go1_depth
			y = hcd_depth - go2_depth
			sv_a = Downward_Semantic_Value(go_id1, go, x)
			sv_b = Downward_Semantic_Value(go_id2, go, y)
			intersecting_terms = intersection(sv_a,sv_b)
			numerator = sum([x for t in intersecting_terms for x in t])
			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b)
			sim_down = (numerator/denominator)
		else:
			sim_down = 0
		if show_plot:
			_plot_go_subgraph(go_id1, go_id2, method, hcd=hcd, sv_a=sv_a, sv_b=sv_b)
		if cache_key is not None:
			_term_similarity_cache[cache_key] = sim_down
		return sim_down
def Similarity_of_Set_of_GOTerms(set1, set2, method):
	method = _normalize_similarity_method(method)
	set1 = _sanitize_go_terms(set1)
	set2 = _sanitize_go_terms(set2)
	if len(set1) == 0 or len(set2) == 0:
		return np.nan
	cache_key = _set_similarity_cache_key(set1, set2, method)
	if cache_key in _set_similarity_cache:
		return _set_similarity_cache[cache_key]

	Sim1 = []
	Sim2 = []
	for idx, goterm in enumerate(set1):
		Sim1.append([])
		for goid in set2:
			similarity = Similarity_of_Two_GOTerms(goterm, goid, go, method)
			if not pd.isna(similarity):
				Sim1[idx].append((goterm, goid, similarity))
	
	for idx, goterm in enumerate(set2):
		Sim2.append([])
		for goid in set1:
			similarity = Similarity_of_Two_GOTerms(goterm, goid, go, method)
			if not pd.isna(similarity):
				Sim2[idx].append((goterm, goid, similarity))
	sem1 = [max(sim_list, key=lambda x: x[2]) for sim_list in Sim1 if len(sim_list) > 0]
	sem2 = [max(sim_list, key=lambda x: x[2]) for sim_list in Sim2 if len(sim_list) > 0]
	if len(sem1) == 0 or len(sem2) == 0:
		return np.nan

	similarity = (sum(x[2] for x in sem1) + sum(x[2] for x in sem2)) / (len(set1) + len(set2))
	similarity = round(similarity, 3)
	_set_similarity_cache[cache_key] = similarity
	return similarity
def Similarity_Matrix(genes, method):
	sim_matrix = [[np.nan for _ in genes] for _ in genes]
	for row_idx, gene_a in enumerate(genes):
		for col_idx in range(row_idx, len(genes)):
			gene_b = genes[col_idx]
			similarity = Similarity_of_Set_of_GOTerms(gene_a[1], gene_b[1], method)
			sim_matrix[row_idx][col_idx] = similarity
			sim_matrix[col_idx][row_idx] = similarity
	return sim_matrix

def Agglomerative_Clustering(pathway, Genes, n_clusters, method):
	data = Similarity_Matrix(Genes, method)
	data_matrix = []
	for row in data:
		data_matrix.append([(lambda x: 1-x)(x) for x in row])
	my_data = pd.DataFrame(data = data_matrix, index=[x[0] for x in Genes],columns=[x[0]for x in Genes])
	return Agglomerative(my_data,Genes, pathway, n_clusters)

def Agglomerative( data, Genes, pathway, n_clusters):
	model = AgglomerativeClustering(n_clusters, affinity='precomputed', linkage='complete').fit_predict(data)
	GeneNames = [gene.GeneName for gene in Genes]
	return model.tolist()

def ARI_Score(label_pred, labels_true):
	return adjusted_rand_score(label_pred, labels_true)

def AMI_Score(label_pred, labels_true):
	return adjusted_mutual_info_score(label_pred, labels_true, average_method='arithmetic')

def cont_matrix(y_true, y_pred):
	contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
	return(contingency_matrix)

def purity_score(y_true, y_pred):
	contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
	return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix)
