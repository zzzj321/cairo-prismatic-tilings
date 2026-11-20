"""Exact finite checks for four matrix and four word witnesses.

These checks validate specified periodic inputs. The passage from coarse data
to the complete Euclidean symmetry group is an analytic lemma in the paper,
not a theorem proved by this program. No third-party dependencies.
"""
from fractions import Fraction
from itertools import product
from math import lcm
from pathlib import Path
import json
import platform
import sys
import time

MATS = {
    "I": ((1, 0), (0, 1)), "R90": ((0, -1), (1, 0)),
    "R180": ((-1, 0), (0, -1)), "R270": ((0, 1), (-1, 0)),
    "sx": ((-1, 0), (0, 1)), "sy": ((1, 0), (0, -1)),
    "swap": ((0, 1), (1, 0)), "negative_swap": ((0, -1), (-1, 0)),
}
CORNERS = ((0, 0), (1, 0), (1, 1), (0, 1))
INPUTS = {
    "p4": [[int(q == 2*p % 5) for q in range(5)] for p in range(5)],
    "p4g": [[1, 0], [0, 1]],
    "cmm": [[1, 1, 0, 0], [0, 0, 1, 1]],
    "cm": [[1,1,0,1,0,0], [1,1,0,0,0,1], [0,1,1,0,0,1],
           [1,0,0,1,1,0], [0,0,1,1,1,0], [0,0,1,0,1,1]],
}


def mv(m, v):
    return tuple(sum(m[i][j]*v[j] for j in range(2)) for i in range(2))


def det(m):
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]


def check_matrix(name, rows):
    r, n = len(rows), len(rows[0])
    x = lambda p, q: rows[p % r][q % n]
    eta = lambda p, q: 2*x(p, q)-1
    assert len({sum(row) for row in rows}) == 1
    assert len({sum(x(p,q) for p in range(r)) for q in range(n)}) == 1
    d = sum(x(p, 0) for p in range(r))
    assert 0 < d < r
    D = [[sum(x(p-1,j)-x(p,j) for j in range(q)) for q in range(n)] for p in range(r)]
    R = [[sum(x(i,q)-x(i,q-1) for i in range(p)) for q in range(n)] for p in range(r)]
    ap = [-min(row) for row in D]
    bq = [-min(R[p][q] for p in range(r)) for q in range(n)]
    E = [[D[p][q]+ap[p] for q in range(n)] for p in range(r)]
    F = [[R[p][q]+bq[q] for q in range(n)] for p in range(r)]
    A, B = sum(ap), sum(bq)
    energy = sum(E[p][q]*F[p][q] for p,q in product(range(r), range(n)))
    assert energy == A*B-Fraction(d*n*(r-d), r) == 0
    e = lambda p,q: E[p % r][q % n]
    f = lambda p,q: F[p % r][q % n]

    def qcorners(i, j):
        if (i+j) % 2 == 0:
            return (0,)*4
        p, q = i//2, j//2
        if i % 2 == 0:
            return (-eta(p-1,q), eta(p,q), -eta(p,q), eta(p-1,q))
        return (eta(p,q-1), -eta(p,q-1), eta(p,q), -eta(p,q))

    def edge_length(v, w):
        i, j = min(v[0],w[0]), min(v[1],w[1])
        if v[1] == w[1]:
            assert abs(v[0]-w[0]) == 1
            return 1+e(i//2,j//2) if i % 2 == 0 else 1
        assert v[0] == w[0] and abs(v[1]-w[1]) == 1
        return 1+f(i//2,j//2) if j % 2 == 0 else 1

    # For unequal rectangular periods, testing just r*n source points can miss
    # exchanged-axis incompatibility. L*L is a common period for every predicate.
    L = lcm(r,n)
    phases = {}
    corner_comparisons = metric_edge_comparisons = 0
    accepted = []
    for label, m in MATS.items():
        mone = mv(m, (1,1))
        k0 = tuple((v-1)//2 for v in mone)
        delta = tuple((1-v)//2 for v in mone)
        corner_perm = [CORNERS.index(tuple(u+v for u,v in zip(mv(m,c),delta))) for c in CORNERS]
        phases[label] = []
        for a,b in product(range(r), range(n)):
            def odd_image(p,q):
                z = mv(m,(p,q))
                return z[0]+a+k0[0], z[1]+b+k0[1]
            signed_height = all(eta(*odd_image(p,q)) == det(m)*eta(p,q)
                                for p,q in product(range(L),repeat=2))
            # Independent check: transform each actual quadrilateral corner;
            # compare its angle using the separately written matrix-angle formula.
            corner_ok = True
            for i,j in product(range(2*L), repeat=2):
                fi,fj = mv(m,(i,j))
                dest = qcorners(fi+2*a, fj+2*b)
                src = qcorners(i,j)
                for c in range(4):
                    corner_ok &= src[c] == dest[corner_perm[c]]
                    corner_comparisons += 1
            assert signed_height == corner_ok, (name,label,a,b,"angle mismatch")
            if not signed_height:
                continue
            phases[label].append([a,b])
            accepted.append((label,a,b))
            def vertex_image(v):
                z = mv(m,v)
                return z[0]+2*a+delta[0], z[1]+2*b+delta[1]
            for i,j in product(range(2*L), repeat=2):
                v = (i,j)
                for step in ((1,0),(0,1)):
                    w = (i+step[0],j+step[1])
                    assert edge_length(v,w) == edge_length(vertex_image(v),vertex_image(w)), (name,label,a,b,v,w)
                    metric_edge_comparisons += 1
    # Validate closure of the accepted affine actions, including their shifts.
    accepted_set = set(accepted)
    inverse_names = {m: label for label,m in MATS.items()}
    def k(label,a,b):
        one = mv(MATS[label],(1,1))
        return a+(one[0]-1)//2, b+(one[1]-1)//2
    for l1,a1,b1 in accepted:
        for l2,a2,b2 in accepted:
            m1,m2 = MATS[l1],MATS[l2]
            mm = tuple(tuple(sum(m1[i][s]*m2[s][j] for s in range(2)) for j in range(2)) for i in range(2))
            kk1,kk2 = k(l1,a1,b1),mv(m1,k(l2,a2,b2))
            one = mv(mm,(1,1))
            aa = (kk1[0]+kk2[0]-(one[0]-1)//2) % r
            bb = (kk1[1]+kk2[1]-(one[1]-1)//2) % n
            assert (inverse_names[mm],aa,bb) in accepted_set
    index = len(phases["I"])
    assert all(len(v) in (0,index) for v in phases.values())
    assert 4*r*n % index == 2*(n*A+r*B) % index == 0
    expected = {"p4": {"I","R90","R180","R270"},
                "p4g": set(MATS), "cmm": {"I","R180","sx","sy"},
                "cm": {"I","sx"}}[name]
    assert {label for label,v in phases.items() if v} == expected
    return {"input": rows,"r":r,"n":n,"d":d,"a":ap,"b":bq,"A":A,"B":B,
            "energy":energy,"E":E,"F":F,"cover_L":L,"valid_phases":phases,
            "counts":{label:len(v) for label,v in phases.items()},
            "translation_index":index,"primitive_C":4*r*n//index,
            "primitive_P":2*(n*A+r*B)//index,
            "corner_comparisons":corner_comparisons,"metric_edge_comparisons":metric_edge_comparisons,
            "affine_closure_pairs":len(accepted)**2}


def check_word(name, word):
    N = next(k for k in range(1,len(word)+1) if all(word[i] == word[(i+k)%len(word)] for i in range(len(word))))
    inv = [t for t in range(N) if all(word[i] == word[(t-i)%N] for i in range(N))]
    if N % 2:
        group = "pgg" if inv else "pg"
    elif not inv:
        group = "p1"
    else:
        assert len(inv) == 1
        group = "p2" if inv[0] % 2 == 0 else "pg"
    assert group == name
    ell = lcm(N,2)
    return {"word":word,"minimal_period":N,"inversion_phases":inv,
            "ell":ell,"primitive_C":2*ell,
            "primitive_P":2*sum(word[i%N]-1 for i in range(ell)),"group":group}


def main():
    start = time.perf_counter()
    matrices = {name:check_matrix(name,x) for name,x in INPUTS.items()}
    words = {name:check_word(name,word) for name,word in
             {"p1":[1,1,2,3],"p2":[1,1,1,2],"pg":[1,1,2,2],"pgg":[2]}.items()}
    result = {"schema":"round32-witness-actions-v1","python":platform.python_version(),
              "matrices":matrices,"words":words,"elapsed_seconds":time.perf_counter()-start,
              "evidence":"exact-check of specified inputs, conditional on analytic reconstruction/group criterion",
              "scope":"No universal enumeration or formal certification of the classification."}
    path = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name("group-checks-round32.json")
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"output":str(path),"elapsed_seconds":result["elapsed_seconds"],
                      "matrix_counts":{k:v["counts"] for k,v in matrices.items()},"words":words},ensure_ascii=False))


if __name__ == "__main__":
    main()
