# Checks

Each claim in the paper and the log that supports it (`checks/logs/`). All checks were run on 18 November 2026 with
the srg conda environment (Python 3.12, numpy 2.3, python-igraph 1.0, GAP 4).

The verification code in `code/` is independent of the original construction code in `original/`. It was written
from the mathematical description only and shares no code with the original scripts.

## Strong regularity

| claim | check | log |
|---|---|---|
| The original Cayley and three-orbit graphs are srg(768,312,136,120); the complement is srg(768,455,262,280) | `verify_all.sh` (all 294,528 pairs each) | `checks/verify_all.out` |
| The norm-coordinate Cayley and shear graphs and 20 random members | `verify_all.sh` | `checks/verify_all.out` |
| `code/srg768.py` rebuilds both constructions and their complements, all pairs; the 3×3 block spectrum is 312^1 24^182 (−8)^585 | `srg768.py` | `srg768_main.log` |
| The rebuilt graphs are isomorphic to the original files (canonical hashes; also \|Aut\| and four-cliques) | `invariants768.py` | `invariants768.log` |
| The original script rebuilds the three original files byte for byte | `original/construct768.py`, `cmp` | (interactive; see `original/README.md`) |
| The original scripts in `original/jordan/` rebuild `crossratio.json`, `planes.json`, the Jordan graph file and its certificate byte for byte | `cayley_control.py`, `cmp` | (interactive; see `original/jordan/README.md`) |
| Every random member in the isomorphism count was verified on all pairs before hashing | `prolific768.py` | `prolific768.jsonl` |

## Seeds

| claim | check | log |
|---|---|---|
| Lemma 3.5 (Gauss periods), Lemma 3.6 (parity), all 48 norm seeds | `srg768.parity_lemma`, `gauss_period_lemma`, `is_seed` by enumeration | `norm_lemmas.log` |
| C5-invariant seeds: 4080 in total; 12 on each of the 340 subline four-sets; none on the other 2040 four-sets | `seed_search.py` (all 2380 four-sets × 4096 sign patterns) | `seed_search.log`, `seeds_base.json` |
| The 4080 C5-invariant seeds are exactly the GL(2,16)-images of one norm seed | `seed_search.py gl` | `seed_gl_orbit.log` |
| Proposition 3.8: every seed is z·χ_c with z C5-invariant; 3072 seeds on each subline four-set, none on the other orbit. One four-set per PΓL(2,16)-orbit (orbits by cross-ratio class, `seed_search.log`); CP-SAT (OR-Tools 9.15, 4 threads) with all known seeds excluded reports INFEASIBLE for both. CP-SAT gives no infeasibility certificate | `newseeds.py 0,1,2,3 4`, `newseeds.py 0,1,6,7 4` | `newseeds_0_1_2_3.log`, `newseeds_0_1_6_7.log` |
| The original seed is z_{0,0,1} ∘ M0^{-1} with M0 = ((t²+1, t²), (0, 1)) (Section 7) | `original_seed_check.py` | `original_seed_check.log` |

## Cayley graph and partial difference set

| claim | check | log |
|---|---|---|
| The spectral construction equals Cay(G, D); D is a (768,312,136,120) PDS; G is nonabelian | `cayley768.py` (all 312² quotients) | `cayley768.log` |
| Original Cayley graph: G = SmallGroup(768,1085323) = C2^4 × (C2^4 ⋊ C3); Aut = C2^8 ⋊ C15, of order 3840, solvable; vertex stabilizer C5; a single class of subgroups of order 768, all regular. For the norm-seed Cayley graph (Proposition 6.2) the paper proves this from \|Aut\| = 3840 (bliss, `invariants768.log`) | GAP on `cayley768.g` | `cayley768_gap.log` |
| `graphs/pds768_norm_cayley.txt` is a PDS in the group given by its own header; IdGroup is the same | `check_pds.py`, GAP | `check_pds.log` |

## Invariants and isomorphism classes

| claim | check | log |
|---|---|---|
| \|Aut\| and orbits of the four explicit graphs; four-cliques through each fiber; the norm and original Cayley graphs are not isomorphic | `invariants768.py` | `invariants768.log` |
| 1000 random members, each verified on all pairs: 999 isomorphism classes (trials 308 and 389 coincide); \|Aut\| = 256 for 996, 512 for 2 (orbits 256+512), 1280 for 2 | `prolific768.py 1000 2026` (bliss canonical form; self-test under relabeling) | `prolific768.jsonl` |

## Template

| claim | check | log |
|---|---|---|
| Parameter sets of the n-fiber template (Table 2 of the paper), with their status in the local copy of Brouwer's table | `families.py 256` | `families.log` |
| The q = 3 instance gives srg(27,16,10,8) | `template_q3.py` | `template_q3.log` |
| Proposition 8.1 (prime fields): the only flat line-sum sequences are a full line and its complement (exhaustive, p = 3, 5, 7) | `prime_obstruction.py` | `prime_obstruction.log` |
| Relaxed template, Remark 8.2 (arbitrary two-eigenvalue Cayley graphs in fibers F_p^2; the code takes the fiber eigenvalue off the supports to be the negative eigenvalue, the other choice gives the complementary parameter sets): every open row with v ≤ 1300 that it reaches has fibers of order p² with p an odd prime and \|H1_B\| = p, so Proposition 8.1 excludes it | `relaxed_template.py current_table.csv` (local copy of Brouwer's table as of 18 Nov 2026, with our own constructions marked as known; not in the repository) | `relaxed_template.log` |

## Section 9: generalized seeds and the Jordan Cayley graph

The graph `graphs/srg768_312_136_120_jordan_cayley.g6` and its certificate `data/jordan_seed.json` (60 signed
frequencies, images of the basis under M) come from the search; F_2^8 with the dot product, 8-bit integers.

| claim | check | log |
|---|---|---|
| The certificate seed is supported on four 4-dim subspaces meeting pairwise in 0, of type (c) in all 24 orders; Hz in {10, −6}; M has order 3 and one fixed point; the twelve subspaces meet pairwise in 0; the graph rebuilt from the certificate equals the file; srg(768,312,136,120) and complement on all pairs; τ is an automorphism; D is a PDS in G_J (all 312² quotients); left multiplications are automorphisms | `jordan768.py --gap jordan768.g` | `jordan768.log` |
| \|Aut Γ_J\| = 768, one orbit; 444,352 four-cliques per vertex | `jordan768.py` (bliss) | `jordan768.log` |
| G_J = SmallGroup(768,1085321), center trivial, C2^8 : C3 | GAP on `jordan768.g` | `jordan768_gap.log` |
| Proposition 9.4: each C_i ∪ {0} contains exactly 8 four-dim subspaces, 21 of their 70 four-sets of type (c). Cross-check: every C′ of Construction 4.1 with exactly 8 subspaces has no type-(c) four-set | `jordan768.py` | `jordan768.log` |
| Table 3 and the matrices Q_1, Q_2 (M(x) = Q_1 x + Q_2 x̄; printed as A, B by the script): the certificate transported to X, Y, Δ, J in F_4^4; the graph built from Table 3 is srg on all pairs and has the same canonical hash as the certificate graph | `jordan_normal.py` | `jordan_normal.log`, `data/jordan_normal.json` |
| Proposition 9.5 (groups): proof only. Proposition 9.6: M = I_2 + W + W + W on F_2^8 (fixed space of dimension 2); random sampling (seed 1) found a type-(c) four-set with the twelve subspaces meeting pairwise in 0 after 6,110,218 samples; CP-SAT found a seed; the graph is srg on all pairs, τ is an automorphism, D is a PDS, \|Aut\| = 768 with one orbit, 443,812 four-cliques per vertex | `dim2_search.py 1 20000000` | `dim2_search.log`, `data/dim2_seed.json` |
| Γ_2 rebuilt from `data/dim2_seed.json` and written to `graphs/srg768_312_136_120_dim2_cayley.g6`; support of type (c); each C_i ∪ {0} contains 9 four-dim subspaces (Construction 4.1: 8, 13 or 18, `jordan768.log`); not isomorphic to Γ_J | `dim2_search.py rebuild` | `dim2_rebuild.log` |
| G_2 = SmallGroup(768,1085320) = C2^2 × (C2^6 : C3), center of order 4 | GAP on `dim2.g` | `dim2_gap.log` |
| Proposition 9.2: exactly 12 twist classes (3072 seeds) on a type-(c) four-set; iterated CP-SAT, final round INFEASIBLE. Run on X, Y, Δ, J (normal form) and, as a cross-check, on the certificate support | `jordan_seeds.py 4 normal`, `jordan_seeds.py 4` | `jordan_seeds_normal.log`, `jordan_seeds.log` |

## Bugs found and fixed while checking

- The first version of the norm three-orbit example used coset representatives 0, t, t². But t and t² lie in the same
  additive coset of F_4 = {0, 1, t²+t, t²+t+1}. The disjointness assertion in `assemble` caught this. The
  representatives are now 0, t, t³, t³+t.
