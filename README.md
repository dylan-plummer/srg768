# Strongly regular graphs with parameters (768, 312, 136, 120)

This repository holds the code, graphs and a draft paper for the result below (Dylan Plummer, draft, November 2026).
The paper is in `paper/`.

## Results

| claim | status |
|---|---|
| srg(768, 312, 136, 120) and its complement srg(768, 455, 262, 280) exist | listed as open (`?`) in Brouwer's table, checked 18 Nov 2026 |
| One of them is a Cayley graph on G = C2^4 × (C2^4 ⋊ C3) = SmallGroup(768, 1085323), so G contains a (768, 312, 136, 120) partial difference set | no abelian group can contain one (Ma: cyclic Sylow 3-subgroup) |
| The construction gives many non-isomorphic graphs | 1000 random members fall into 999 isomorphism classes |

Spectrum: 312^1, 24^182, (−8)^585.

- **Literature check.** We found no earlier construction of any srg on 768 vertices, and no partial difference set of
  order 768. Brouwer's table lists every feasible 768 parameter set as open or ruled out. Section 1 of the paper gives
  the details.
- **Automorphisms.**
  - Cayley graph: Aut = C2^8 ⋊ C15, of order 3840. It is solvable, and its subgroups of order 768 form a single
    conjugacy class, all regular. So the graph is a Cayley graph only of G.
  - Three-orbit graph: |Aut| = 1280, with three orbits of size 256.
  - Random members: typically |Aut| = 256 (the translations).
- **Template, and why prime fields fail.** The same proof applies formally to n fibres F_q^2 on K_n. The parameter
  arithmetic then allows five further open sets with v ≤ 1300: (363,192,106,96), (507,176,70,56),
  (1083,560,298,280), (196,120,74,72) and (676,459,314,306). All have prime q, and Proposition 8.1 of the paper shows
  that over F_p^2 the cross sets must be cosets of lines or their complements, so the template cannot produce any of
  them (`code/families.py`, `code/prime_obstruction.py`, `code/relaxed_template.py`; Section 8 of the paper).

## Construction

Notation: F = F_16 ⊃ K = F_4, N(a) = a^5 is the norm F → K, and V = F^2 is paired with itself by
⟨(x,y),(a,b)⟩ = Tr(ax + by). The Walsh transform is Hf(v) = Σ_β f(β)(−1)^⟨v,β⟩.

1. **Seed.** A function z, ±1 on four punctured spread lines and 0 elsewhere, with Hz ∈ {10, −6}. Then
   B = {Hz = 10} has |B| = 96 and Fourier transform 16z + 96δ.
   - **Norm seeds.** z(a, μa) = ε_μ · (−1)^[N(a) = ν] for μ ∈ K, with ε = −1 at exactly one μ0.
   - **Proof that these are seeds.**
     - Gauss periods: each line contributes ε_μ(−3 + 8δ_μ).
     - Parity: for each (x, y), the number of μ ∈ K with N(x + μy) ∈ {0, ν⁻¹} is even.
   - **All C5-invariant seeds.** There are 4080 of them, forming one GL(2,16)-orbit. They live exactly on the 340
     four-subsets of sublines PG(1,4).
   - **All seeds.** Every seed is z·χ_c for a C5-invariant seed z and a character χ_c, so there are
     340 × 3072 = 1,044,480 seeds, and the sets B form a single AGL(2,16)-orbit (Proposition 3.8; CP-SAT
     infeasibility, `code/newseeds.py`).
2. **Graph.** Take three fibres V and three seeds z_i with pairwise disjoint support lines S_i.
   - Fibre i is joined internally by the punctured annihilators of the eight lines S_{i−1} ∪ S_i.
   - (i, x) ~ (i+1, y) iff x + y ∈ B_i.
3. **Proof.** On the three characters χ_β, one per fibre, A acts by a 3×3 matrix:
   - β = 0: diagonal 120, off-diagonal 96;
   - β on a support line: a 2×2 block [[8, ±16], [±16, 8]] plus a (−8);
   - every other β: −8·I.
4. **Cayley member.** M = [[1,0],[c,ω]] with c ∉ K has order 3 and cycles three additive cosets of K. Set
   z_i = z_0 ∘ M^{−i}. Then (i, x) ↦ (i+1, M^{−T}x), together with the translations, acts regularly.

## Layout

| path | contents |
|---|---|
| `graphs/` | graph6 files, all verified on every pair by `verify_all.sh`; the partial difference set as a text file |
| `code/srg768.py` | the construction: F_16, Walsh transform, seeds (original and norm), assembly, 3×3 spectral certificate, proof lemmas checked by enumeration |
| `code/verify_srg.py` | independent numpy verifier (A² identity on all pairs); reads multi-graph graph6 files |
| `code/invariants768.py` | \|Aut\| and orbits (bliss), four-cliques per fibre, canonical hashes; isomorphism with the original files |
| `code/seed_search.py` | exhaustive search for C5-invariant seeds; `gl` argument: GL(2,16)-orbit check |
| `code/newseeds.py` | classification of all seeds: CP-SAT with the known seeds excluded, one four-set per orbit (needs OR-Tools) |
| `code/cayley768.py`, `code/check_pds.py` | Cayley group law and PDS counts; GAP scripts for IdGroup and Aut |
| `code/prolific768.py` | random members, isomorphism-class count |
| `code/families.py`, `code/template_q3.py` | parameter analysis of the template; the q = 3 instance (Schläfli) |
| `code/export_graphs.py` | regenerates the norm-coordinate graphs, the sample and the PDS file |
| `original/` | the original construction and seed-search scripts (byte-identical rebuild of the original graphs) |
| `checks/` | logs of every check cited in the paper; see `CHECKS.md` |
| `paper/` | draft preprint `srg768.tex` |

### Graph files

| file | graph |
|---|---|
| `srg768_312_136_120_cayley.g6` | first Cayley graph (original coordinates); \|Aut\| = 3840; 445,040 four-cliques per vertex |
| `srg768_312_136_120_threeorbit.g6` | first three-orbit graph (original coordinates); \|Aut\| = 1280; four-cliques per vertex 444,230 / 443,250 / 443,760 by fibre |
| `srg768_455_262_280_threeorbit_complement.g6` | its complement |
| `srg768_312_136_120_norm_cayley.g6` | Cayley graph in norm coordinates (c = t, μ0 = 0, ν = 1); \|Aut\| = 3840; not isomorphic to the first Cayley graph |
| `srg768_312_136_120_norm_threeorbit.g6` | shear graph in norm coordinates (cosets 0, t, t^3 of K); \|Aut\| = 1280 |
| `srg768_312_136_120_family_sample20.g6` | 20 pairwise non-isomorphic random members (one graph per line) |
| `pds768_norm_cayley.txt` | the partial difference set of the norm Cayley graph, with the group law in the header |

Integer encoding: F_16 = F_2[t]/(t^4+t+1), with the bits of an integer 0..15 giving the coefficients of
1, t, t^2, t^3. A vector (x, y) is encoded as 16x + y, and vertex (i, v) has index 256i + v.

## Reproduce

Python 3.10+ with numpy. The invariant scripts also need python-igraph, the seed classification needs OR-Tools, and the group checks need GAP.

```sh
./verify_all.sh                         # all graph files, all pairs
cd code
python srg768.py                        # rebuild both constructions, spectral certificate, all-pairs check
python invariants768.py                 # Aut, four-cliques, isomorphism with original/ files
python seed_search.py; python seed_search.py gl
python newseeds.py 0,1,2,3 4; python newseeds.py 0,1,6,7 4   # about 5 and 40 minutes
python cayley768.py && gap -q ../checks/logs/cayley768.g
python check_pds.py ../graphs/pds768_norm_cayley.txt
python prolific768.py 1000 2026 > prolific.jsonl
python families.py 256
cd ../original && python -B construct768.py --output /tmp/rebuild768
```
