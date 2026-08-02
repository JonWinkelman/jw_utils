"""Weighted categorical mutual information and permutation tests."""

import numpy as np
import pandas as pd


def sequence_matrix(sequences):
    """Convert equal-length strings into a two-dimensional character array."""
    sequences = pd.Series(sequences).astype(str)
    lengths = sequences.str.len().unique()
    if len(lengths) != 1:
        raise ValueError(f"Sequences must have equal lengths; observed {lengths}")
    return np.asarray([list(sequence) for sequence in sequences])


def weighted_mutual_information(x, y, weights):
    """Calculate weighted categorical MI in bits, excluding invalid states."""
    x = np.asarray(x)
    y = np.asarray(y)
    weights = np.asarray(weights, dtype=float)
    valid = (
        np.isin(x, list("ACGT"))
        & ~np.isin(y, ["-", ".", "X", "?"])
        & np.isfinite(weights)
        & (weights > 0)
    )
    if valid.sum() < 2:
        return 0.0
    _, x_codes = np.unique(x[valid], return_inverse=True)
    _, y_codes = np.unique(y[valid], return_inverse=True)
    nx = x_codes.max() + 1
    ny = y_codes.max() + 1
    joint = np.bincount(
        x_codes * ny + y_codes,
        weights=weights[valid],
        minlength=nx * ny,
    ).reshape(nx, ny)
    joint = joint / joint.sum()
    expected = joint.sum(axis=1, keepdims=True) * joint.sum(axis=0, keepdims=True)
    occupied = joint > 0
    return float(
        np.sum(joint[occupied] * np.log2(joint[occupied] / expected[occupied]))
    )


def weighted_mi_matrix(dna_matrix, protein_matrix, weights):
    """Calculate protein-position by DNA-position weighted MI."""
    return np.asarray(
        [
            [
                weighted_mutual_information(
                    dna_matrix[:, dna_column],
                    protein_matrix[:, protein_column],
                    weights,
                )
                for dna_column in range(dna_matrix.shape[1])
            ]
            for protein_column in range(protein_matrix.shape[1])
        ]
    )


def shuffled_indices(rng, groups):
    """Return indices shuffled independently inside each supplied group."""
    groups = np.asarray(groups)
    indices = np.arange(len(groups))
    shuffled = indices.copy()
    for group in pd.unique(groups):
        members = indices[groups == group]
        shuffled[members] = rng.permutation(members)
    return shuffled


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg correction retaining the input array shape."""
    p_values = np.asarray(p_values)
    flat = p_values.ravel()
    order = np.argsort(flat)
    ranked = flat[order]
    adjusted = ranked * len(flat) / np.arange(1, len(flat) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result.reshape(p_values.shape)


def permutation_mi_test(
    dna_matrix,
    protein_matrix,
    weights,
    groups,
    n_permutations=500,
    weight_travels_with_protein=False,
    seed=2026,
):
    """Calculate observed MI and cell-wise/FWER permutation significance."""
    weights = np.asarray(weights, dtype=float)
    observed = weighted_mi_matrix(dna_matrix, protein_matrix, weights)
    exceedances = np.zeros_like(observed, dtype=int)
    null_maxima = np.empty(n_permutations, dtype=float)
    rng = np.random.default_rng(seed)

    for permutation_number in range(n_permutations):
        order = shuffled_indices(rng, groups)
        permuted_weights = weights[order] if weight_travels_with_protein else weights
        permuted = weighted_mi_matrix(
            dna_matrix,
            protein_matrix[order],
            permuted_weights,
        )
        exceedances += permuted >= observed
        null_maxima[permutation_number] = permuted.max()

    p_values = (exceedances + 1) / (n_permutations + 1)
    max_stat_p = (
        (null_maxima[:, None, None] >= observed[None, :, :]).sum(axis=0) + 1
    ) / (n_permutations + 1)
    return {
        "mi": observed,
        "p_value": p_values,
        "q_value": benjamini_hochberg(p_values),
        "max_stat_p": max_stat_p,
        "null_maxima": null_maxima,
    }


def weighted_conditional_table(dna, amino_acid, weights):
    """Return P(DNA base | amino-acid state) using weighted observations."""
    frame = pd.DataFrame(
        {"dna": dna, "amino_acid": amino_acid, "weight": weights}
    )
    frame = frame.loc[
        frame["dna"].isin(list("ACGT"))
        & ~frame["amino_acid"].isin(["-", ".", "X", "?"])
    ]
    table = frame.pivot_table(
        index="amino_acid",
        columns="dna",
        values="weight",
        aggfunc="sum",
        fill_value=0,
    ).reindex(columns=list("ACGT"), fill_value=0)
    return table.div(table.sum(axis=1), axis=0)
