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
| Every random member in the isomorphism count was verified on all pairs before hashing | `prolific768.py` | `prolific768.jsonl` |

## Seeds

| claim | check | log |
|---|---|---|
| Lemma 3.4 (Gauss periods), Lemma 3.5 (parity), all 48 norm seeds | `srg768.parity_lemma`, `gauss_period_lemma`, `is_seed` by enumeration | `norm_lemmas.log` |
| C5-invariant seeds: 4080 in total; 12 on each of the 340 subline four-sets; none on the other 2040 four-sets | `seed_search.py` (all 2380 four-sets × 4096 sign patterns) | `seed_search.log`, `seeds_base.json` |
| The 4080 C5-invariant seeds are exactly the GL(2,16)-images of one norm seed | `seed_search.py gl` | `seed_gl_orbit.log` |

## Cayley graph and partial difference set

| claim | check | log |
|---|---|---|
| The spectral construction equals Cay(G, D); D is a (768,312,136,120) PDS; G is nonabelian | `cayley768.py` (all 312² quotients) | `cayley768.log` |
| G = SmallGroup(768,1085323) = C2^4 × (C2^4 ⋊ C3); Aut = C2^8 ⋊ C15, of order 3840, solvable; vertex stabiliser C5; a single class of subgroups of order 768, all regular | GAP on `cayley768.g` | `cayley768_gap.log` |
| `graphs/pds768_norm_cayley.txt` is a PDS in the group given by its own header; IdGroup is the same | `check_pds.py`, GAP | `check_pds.log` |

## Invariants and isomorphism classes

| claim | check | log |
|---|---|---|
| \|Aut\| and orbits of the four explicit graphs; four-cliques through each fibre; the norm and original Cayley graphs are not isomorphic | `invariants768.py` | `invariants768.log` |
| 1000 random members, each verified on all pairs: 999 isomorphism classes (trials 308 and 389 coincide); \|Aut\| = 256 for 996, 512 for 2 (orbits 256+512), 1280 for 2 | `prolific768.py 1000 2026` (bliss canonical form; self-test under relabelling) | `prolific768.jsonl` |

## Template

| claim | check | log |
|---|---|---|
| Parameter sets of the n-fibre template (Table 2 of the paper), with their status in the local copy of Brouwer's table | `families.py 256` | `families.log` |
| The q = 3 instance gives srg(27,16,10,8) | `template_q3.py` | `template_q3.log` |

## Bugs found and fixed while checking

- The first version of the norm three-orbit example used coset representatives 0, t, t². But t and t² lie in the same
  additive coset of F_4 = {0, 1, t²+t, t²+t+1}. The disjointness assertion in `assemble` caught this. The
  representatives are now 0, t, t³, t³+t.
