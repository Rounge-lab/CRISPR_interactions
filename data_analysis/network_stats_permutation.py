#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Network weights, connectivity and randomness

"""

import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Dict, Tuple, List, Optional

wdir='PATH_TO_MANUS_FOLDER'

nodes=pd.read_csv('/'.join([wdir,'results/Nodes_network_crisprfree.csv']))
edges=pd.read_csv('/'.join([wdir,'results/Edges_network_crisprfree.csv']))

def validate_edges_df(df,b_col,v_col, w_col)-> pd.DataFrame:
  
    missing = [c for c in [b_col, v_col, w_col] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    out = df[[b_col, v_col, w_col]].copy()
    out = out.dropna(subset=[b_col, v_col, w_col])

    # Enforce integers
    out[w_col] = pd.to_numeric(out[w_col], errors="raise").astype(np.int64)
    if (out[w_col] <= 0).any():
        bad = out.loc[out[w_col] <= 0].head(5)
        raise ValueError(f"All {w_col} must be positive integers. Examples:\n{bad}")

    # Collapse duplicates
    out = out.groupby([b_col, v_col], as_index=False)[w_col].sum()

    return out

def compute_network_statistics(df,b_col,v_col, w_col)-> Dict:
    """
    Compute U, Mk, mean multiplicity, HHI, Gini
    from a weighted bipartite edge table.
    """
    counts = df[w_col].to_numpy(dtype=float)
    N = counts.sum()
    U = len(counts)

    # M_k
    Mk = dict(pd.Series(counts.astype(int)).value_counts().sort_index())

    # HHI
    p = counts / N
    HHI = float(np.sum(p * p))

    # Gini
    x = np.sort(counts)
    n = len(x)
    i = np.arange(1, n + 1)
    Gini = float((2 * np.sum(i * x) / (n * np.sum(x))) - (n + 1) / n)

    return {
        "N": int(N),
        "U": int(U),
        "mean_mult": float(N / U),
        "HHI": HHI,
        "Gini": Gini,
        "Mk": Mk,
    }

def expand_to_observations(df,b_col,v_col, w_col)-> pd.DataFrame:
    """
    Expand weighted edges into observation-level rows.
    """
    rows = []
    for _, r in df.iterrows():
        rows.extend([(r[b_col], r[v_col])] * int(r[w_col]))

    return pd.DataFrame(rows, columns=[b_col, v_col])

def bipartite_edge_swap(obs_df,b_col,v_col, rng, n_swaps=1)-> pd.DataFrame:
    """
    Perform degree-preserving swaps on observation-level bipartite edges.
    """

    obs = obs_df.copy().reset_index(drop=True)
    N = len(obs)

    for _ in range(n_swaps):
        i, j = rng.integers(0, N, size=2)
        if i == j:
            continue

        b1, v1 = obs.loc[i, [b_col, v_col]]
        b2, v2 = obs.loc[j, [b_col, v_col]]

        if b1 == b2 or v1 == v2:
            continue

        obs.loc[i, v_col] = v2
        obs.loc[j, v_col] = v1

    return obs

def aggregate_observations(obs_df,b_col,v_col,w_col) -> pd.DataFrame:
    """
    Convert observation-level edges back to weighted pairs.
    """
    return (
        obs_df
        .groupby([b_col, v_col], as_index=False)
        .size()
        .rename(columns={"size": w_col})
    )

def generate_nullB_network(df,b_col,v_col, w_col,rng, burnin_swaps_per_edge) -> pd.DataFrame:
    """
    Generate one Null B network preserving all margins exactly.
    """

    obs = expand_to_observations(df,b_col,v_col, w_col)
    N = len(obs)
    n_swaps = int(burnin_swaps_per_edge * N)

    obs = bipartite_edge_swap(obs, b_col,v_col, rng, n_swaps)
    nullB=aggregate_observations(obs, b_col, v_col,w_col)
    
    return nullB

def nullB_permutation_test(df, b_col, v_col, w_col, seed, n_perm, burnin_swaps_per_edge = 5.0):
    """
    Run Null B permutations and collect statistics.
    """
    rng = np.random.default_rng(seed)
    
    print('\n')
    print('Calculating observed network stats')
    observed = compute_network_statistics(df,b_col, v_col, w_col)
    
    print(f"Observed U: {observed['U']}")
    print(f"Observed N: {observed['N']}")
    print(f"Observed mean_mult: {observed['mean_mult']}")
    print(f"Observed HHI: {observed['HHI']}")
    print(f"Observed Gini: {observed['Gini']}")
    print('----------------------------------------------')
      
    rows = []
    Mk_accum = defaultdict(list)
    
    print('\n')
    print('Generating null networks.....')
    for x in range(n_perm):
        
        print(f'Permutation number {x}...')
        df_null = generate_nullB_network(df,b_col,v_col, w_col,rng,burnin_swaps_per_edge)
        st = compute_network_statistics(df_null,b_col,v_col, w_col)

        rows.append({k: st[k] for k in ["U", "HHI", "Gini","mean_mult"]})
        for k, v in st["Mk"].items():
            Mk_accum[k].append(v)

    perm_df = pd.DataFrame(rows)

    pvals = {
        "p_U_less": (1 + (perm_df["U"] <= observed["U"]).sum()) / (1 + n_perm),
        "p_HHI_greater": (1 + (perm_df["HHI"] >= observed["HHI"]).sum()) / (1 + n_perm),
        "p_Gini_greater": (1 + (perm_df["Gini"] >= observed["Gini"]).sum()) / (1 + n_perm),
        "p_mult_greater": (1 + (perm_df["mean_mult"] >= observed["mean_mult"]).sum()) / (1 + n_perm),
    }

    Mk_mean = {k: float(np.mean(v)) for k, v in Mk_accum.items()}

    return perm_df, {"observed": observed, "pvals": pvals, "Mk_null_mean": Mk_mean}

#set the parameters
seed=103
n_perm=200

#Calculate for MGEs
clean_edges=validate_edges_df(edges,'MAG','Taxon','NumInteractions')
perm_stats, obs_stats=nullB_permutation_test(clean_edges,'MAG','Taxon','NumInteractions',seed,n_perm)

perm_stats.to_csv('/'.join([wdir,'results/Edges_network_stats_permutation_crisprfree.csv']), index=False)
pd.DataFrame([obs_stats['observed']]).round(5).to_csv('/'.join([wdir,'results/Edges_network_stats_observed_crisprfree.csv']), index=False)
pd.DataFrame([obs_stats['pvals']]).round(5).to_csv('/'.join([wdir,'results/Edges_network_stats_Pvals_crisprfree.csv']), index=False)

print('Done')

#Plot Mk differences

mkobs=pd.DataFrame([obs_stats['observed']['Mk']]).T
mkobs.columns=['Observed']
mkmean=pd.DataFrame([obs_stats['Mk_null_mean']]).T
mkmean.columns=['Nullmean']
mk=mkobs.merge(mkmean, left_index=True, right_index=True, how='outer')
mk=mk.reset_index()
mk=mk.rename(columns={'index':'multiplicity'})
mk=mk.fillna(0)
mk.to_csv('/'.join([wdir,'results/Multiplicity_Mk_obs_vs_null.csv']), index=False)

toplot=mk.melt(id_vars=['multiplicity'], value_name='Count', var_name='ObsExp')
toplot=toplot.query('Count>0')
fig=sb.barplot(toplot.query('multiplicity<=160'), x='multiplicity', y='Count',hue='ObsExp')
fig.set(yscale='log')
plt.xticks(rotation=90, fontsize=9)
plt.savefig('/'.join([wdir,'results/Multiplicity_Mk_obs_vs_null.pdf']))
