# Verification for periodic mixed Cairo–Prismatic tilings

This repository contains the exact-arithmetic verification programs and input data for the eight constructions in *The Full Wallpaper Groups of Doubly Periodic Mixed Cairo–Prismatic Tilings*. The inputs are four periodic words and four binary matrices. The computations verify each specified construction. The paper gives the proof that these groups exhaust the possibilities.

## Reproduce the results

Use Python 3.9 or later. The scripts use the Python standard library.

```sh
python reproduce.py
```

The runner executes all five scripts in order, checks the computed full groups and primitive-cell counts against `witness-inputs.json`, and saves the generated JSON files to `reproduced/`. Use `python reproduce.py --output run-2` to choose a different output directory. The output directory must be new. The `reference-*.json` files record the reference runs.

| Script | Computation |
| --- | --- |
| `check_groups_round32.py` | Checks the word and matrix inputs, including period, metric, and discrete symmetry data. |
| `develop_matrix_round36.py` | Develops the four binary matrices into exact square–triangle cells. |
| `develop_words_round36.py` | Develops the four words into exact square–triangle cells. |
| `dualize_all_round36.py` | Checks the corresponding pentagonal tilings and writes their coordinates. |
| `check_physical_groups_round36.py` | Enumerates the full isometries of the eight specified tilings using exact arithmetic in $\mathbb{Q}(\sqrt 3)$. |

`witness-inputs.json` lists the inputs in the same order as the construction table in the paper. The five original scripts carry their inputs in the source. The runner checks that the input values reported by both the discrete and the geometric pipeline agree with this list.

## Expected primitive cells

| Full wallpaper group | Cairo pentagons | Prismatic pentagons |
| --- | ---: | ---: |
| `p1` | 8 | 6 |
| `p2` | 8 | 2 |
| `pg` | 8 | 4 |
| `pgg` | 4 | 4 |
| `p4` | 20 | 8 |
| `p4g` | 8 | 4 |
| `cmm` | 16 | 10 |
| `cm` | 72 | 36 |

The physical-group report gives counts in the supplied period cells. To obtain the primitive counts, the runner divides by the index of the supplied translation lattice in the full translation lattice. The groups of the square–triangle constructions agree with those of their Cairo–Prismatic duals by the duality proved in the paper.

## Verification methods

For an input matrix with periods $r,n$, the discrete check tests the signed matrix criterion on a square period cell of side $\operatorname{lcm}(r,n)$. This permits candidate maps that interchange the coordinates. If the action on coarse faces is $i\mapsto Mi+2(a,b)$, its action on coarse vertices is $v\mapsto Mv+2(a,b)+(\boldsymbol 1-M\boldsymbol 1)/2$, where $\boldsymbol 1=(1,1)$. The check also compares the four angles and the edge lengths of each coarse face.

The geometric check uses the developed unit faces and the input period lattice $\Lambda=\mathbb ZU+\mathbb ZV$. Connectedness, oppositely oriented paired edges, and a single cycle of two square and three triangular sectors at each vertex give a local isometry to the period torus. Equality of the total face area with $\lvert\det(U,V)\rvert$ makes the covering degree one. The symmetry search considers twelve rotations and twelve reflections as candidate linear parts. For each linear part $R$, possible images of a fixed C-type vertex give candidate offsets $t$ modulo $\Lambda$. It first determines the full translation lattice $\Gamma$, then checks $RU,RV\in\Gamma$ and matches every transformed face under $x\mapsto Rx+t$ modulo $\Lambda$. A reversing offset class contains a reflection exactly when $(I+R)(t+\lambda)=0$ for some $\lambda\in\Lambda$.
