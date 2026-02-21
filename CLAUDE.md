# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Superconducting circuit quantization tool: symbolically derives Hamiltonians for superconducting qubit circuits using graph theory (spanning trees, fundamental cutsets). Codebase comments and docstrings are in **Chinese (中文)** — maintain this convention.

## Running

```bash
# From the workspace root (parent of this repo), activate the venv
source .venv/Scripts/activate

# Run the main module directly (it has an __main__ demo)
python build_circuit_graph1.py
```

No formal test framework — testing is done by running scripts and inspecting symbolic output.

## Dependencies

`sympy`, `networkx`, `numpy` (all available in the workspace `.venv`).

## Architecture

The entire pipeline lives in `build_circuit_graph1.py` and flows through four stages:

1. **`build_circuit_graph()`** — Takes a list of `Component` objects (and optional `MutualInductance` list), builds a NetworkX `MultiGraph`, computes a minimum spanning tree, and **re-indexes all edges** so tree branches get keys `0..nt-1` and chords get keys `nt..m-1`. All downstream code depends on this ordering.

2. **`fundamental_cut_matrix()`** — Builds the fundamental cutset matrix Q_f (sympy Matrix, size nt×m). For each tree branch, removes it from the tree, finds the connected component, and determines which chords cross the cut (with ±1 orientation).

3. **`build_parameter_matrices()`** — Constructs three diagonal/block matrices from the edge map: capacitance `D_C`, inverse inductance `L_plus` (with mutual inductance support), and Josephson `D_J`.

4. **`calculate_hamiltonian()`** — Derives the full symbolic Hamiltonian: kinetic energy `½ qᵀ M⁻¹ q`, inductive potential `½ Φᵀ L⁺ Φ`, and Josephson terms `-Ej cos(Φk)`. Automatically introduces external flux variables `Phi_ext_k` for each chord (closed loop).

## Key Conventions

- **Edge direction**: Always from lower-numbered node to higher-numbered node, regardless of input order.
- **Spanning tree weights**: C=1, JJ=2, L=3. Capacitors are preferred in the tree so that the charge (kinetic) energy matrix is well-conditioned.
- **Component types**: `'C'` (capacitor), `'L'` (inductor), `'JJ'` (Josephson junction). Values can be symbolic strings (auto-converted via `sp.symbols`) or numeric.
- **Mutual inductance keys**: After re-indexing, stored as `{(smaller_key, larger_key): value}`.
