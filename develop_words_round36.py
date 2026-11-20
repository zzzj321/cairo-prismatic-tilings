"""Exact coarse-grid, HS and Cairo/Prismatic coordinates for periodic words.

Standard-library Fraction/Q(sqrt(3)); Python fallback: cached Julia unavailable.
One small exact construction, no group search and no floating geometric tests.
Coordinates have unit HS edges; CP incircle radius is 1/2, CP area (2+sqrt(3))/4.
Import build(word) for structured exact data; run this file to save four JSONs.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from collections import Counter, defaultdict
from pathlib import Path
from math import lcm
import json


@dataclass(frozen=True)
class K:
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    def __add__(self, x):
        x = k(x)
        return K(self.a+x.a, self.b+x.b)
    __radd__ = __add__

    def __neg__(self): return K(-self.a, -self.b)
    def __sub__(self, x): return self+-k(x)
    def __rsub__(self, x): return k(x)+-self

    def __mul__(self, x):
        x = k(x)
        return K(self.a*x.a+3*self.b*x.b, self.a*x.b+self.b*x.a)
    __rmul__ = __mul__

    def __truediv__(self, x):
        x = k(x)
        return self*K(x.a, -x.b)*(1/(x.a*x.a-3*x.b*x.b))

    def sign(self):
        if not self.a: return (self.b > 0)-(self.b < 0)
        if not self.b: return (self.a > 0)-(self.a < 0)
        if self.a*self.b > 0: return (self.a > 0)-(self.a < 0)
        d = self.a*self.a-3*self.b*self.b
        return ((self.a > 0)-(self.a < 0))*((d > 0)-(d < 0))

    def floor(self):
        n = self.a.numerator//self.a.denominator
        while (self-n).sign() < 0: n -= 1
        while (self-n-1).sign() >= 0: n += 1
        return n

    def json(self): return [str(self.a), str(self.b)]
    def __float__(self): return float(self.a)+3**.5*float(self.b)


def k(x): return x if isinstance(x, K) else K(x)
def point(x, y): return k(x), k(y)
def add(a, b): return a[0]+b[0], a[1]+b[1]
def sub(a, b): return a[0]-b[0], a[1]-b[1]
def mul(a, t): return a[0]*t, a[1]*t
def dot(a, b): return a[0]*b[0]+a[1]*b[1]
def cross(a, b): return a[0]*b[1]-a[1]*b[0]
def center(f): return mul(tuple(sum((p[j] for p in f), K()) for j in (0, 1)), F(1, len(f)))
def area(f): return sum((cross(p, f[(j+1) % len(f)]) for j, p in enumerate(f)), K())/2
def pkey(p): return tuple((z.a, z.b) for z in p)
def fkey(f): return tuple(sorted(f, key=pkey))
def jp(p): return [z.json() for z in p]
def jf(f): return [jp(p) for p in f]
def halfpoint(p):
    values = [2*p[0].a, 2*p[0].b, 2*p[1].a, 2*p[1].b]
    assert all(z.denominator == 1 for z in values)
    return [int(z) for z in values]
def quad(p, e, f): return p, add(p, e), add(add(p, e), f), add(p, f)


S = K(0, 1)
E = (point(1, 0), point(S/2, F(1, 2)))
V = (point(0, 1), point(F(-1, 2), S/2))
ZERO = point(0, 0)


def build(word):
    """Return JSON-serializable ordered faces and exact period/dual certificate.

    word must be its minimal positive period if periods are called primitive.
    Construction/checks work with any positive word, but no primitive search is
    done here. The four default words are minimal by inspection.
    """
    word = tuple(word)
    assert word and all(type(n) is int and n >= 1 for n in word)
    ell = lcm(len(word), 2)
    aa = [word[i % len(word)] for i in range(ell)]
    xx = [ZERO]
    for i, n in enumerate(aa): xx.append(add(xx[-1], mul(E[i % 2], n)))
    yy = [ZERO, V[0], add(V[0], V[1])]
    U, W = xx[-1], yy[-1]
    determinant = cross(U, W)
    assert determinant.sign() > 0

    def reduce(p):
        a = (cross(p, W)/determinant).floor()
        b = (cross(U, p)/determinant).floor()
        return sub(p, add(mul(U, a), mul(W, b)))

    def norm(f):
        c = center(f)
        t = sub(c, reduce(c))
        return tuple(sub(p, t) for p in f)

    coarse, raw_faces = [], []
    for i, n in enumerate(aa):
        for j in range(2):
            base = add(xx[i], yy[j])
            e, v = E[i % 2], V[j]
            polygon = quad(base, mul(e, n), v)
            black = (i+j) % 2 == 0
            local = []
            for t in range(n):
                p = add(base, mul(e, t))
                a, b, c, d = quad(p, e, v)
                if black:
                    assert dot(e, v) == K()
                    local.append((a, b, c, d))
                elif dot(e, v) == K(F(1, 2)):
                    local.extend(((a, b, d), (b, c, d)))
                else:
                    assert dot(e, v) == K(F(-1, 2))
                    local.extend(((a, b, c), (a, c, d)))
                if not black:
                    assert dot(e, v).sign() == (-1 if i % 2 == 0 else 1)
            coarse.append({'index': [i, j], 'type': 'black_rectangle' if black else 'white_even_triangle_strip',
                           'vertices': jf(polygon), 'horizontal_length': n, 'vertical_length': 1,
                           'unit_faces': [jf(f) for f in local]})
            raw_faces.extend(local)

    canonical = {fkey(norm(f)): norm(f) for f in raw_faces}
    assert len(canonical) == len(raw_faces)
    faces = list(canonical.values())
    edges, corners = defaultdict(list), defaultdict(list)
    for index, f in enumerate(faces):
        assert area(f).sign() > 0
        for j, p in enumerate(f):
            nxt, prev = f[(j+1) % len(f)], f[j-1]
            assert dot(sub(nxt, p), sub(nxt, p)) == K(1)
            midpoint = mul(add(p, nxt), F(1, 2))
            shift = sub(midpoint, reduce(midpoint))
            a, b = sub(p, shift), sub(nxt, shift)
            edges[fkey((a, b))].append((a, b, index))
            corners[reduce(p)].append((sub(nxt, p), sub(prev, p), len(f), sub(center(f), p), index))
    adjacency = defaultdict(set)
    for es in edges.values():
        assert len(es) == 2 and es[0][:2] == es[1][:2][::-1]
        a, b = es[0][2], es[1][2]
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen, pending = {0}, [0]
    while pending:
        p = pending.pop()
        for q in adjacency[p]-seen:
            seen.add(q)
            pending.append(q)
    assert len(seen) == len(faces)

    dual, counts = [], Counter()
    for vertex, cs in corners.items():
        assert Counter(c[2] for c in cs) == {3: 3, 4: 2}
        links = {c[0]: c for c in cs}
        assert len(links) == 5
        ray = min(links, key=pkey)
        start, rays, sequence, polygon = ray, set(), [], []
        for _ in range(5):
            assert ray not in rays
            rays.add(ray)
            entry = links[ray]
            sequence.append(entry[2])
            polygon.append(add(vertex, entry[3]))
            ray = entry[1]
        assert ray == start and len(rays) == 5
        typ = 'P' if any(sequence[j] == sequence[(j+1) % 5] == 4 for j in range(5)) else 'C'
        counts[typ] += 1
        assert area(polygon) == (2+S)/4
        for j, p in enumerate(polygon):
            e = sub(polygon[(j+1) % 5], p)
            distance_numerator = cross(e, sub(vertex, p))
            assert distance_numerator.sign() > 0
            assert 4*distance_numerator*distance_numerator == dot(e, e)
        dual.append({'type': typ, 'incenter': jp(vertex), 'vertices': jf(polygon),
                     'hs_face_sizes_ccw': sequence})

    assert Counter(map(len, faces)) == {4: sum(aa), 3: 2*sum(aa)}
    assert counts['C'] == 2*ell and counts['P'] == 2*sum(n-1 for n in aa)
    assert sum((area(f) for f in faces), K()) == determinant
    assert sum((area(tuple(tuple(K(*z) for z in p) for p in f['vertices'])) for f in dual), K()) == determinant
    return {'coordinate_encoding': 'Each scalar [a,b] means Fraction(a)+Fraction(b)*sqrt(3). All face vertices are CCW.',
            'word': list(word), 'coarse_periods': [[ell, 0], [0, 2]],
            'periods': [halfpoint(U), halfpoint(W)],
            'faces': [[halfpoint(p) for p in f] for f in faces],
            'half_integer_interface': 'periods and faces points [a,b,c,d] mean ((a+b*sqrt(3))/2,(c+d*sqrt(3))/2).',
            'periods_qsqrt3': [jp(U), jp(W)], 'period_area': determinant.json(),
            'coarse_faces': coarse, 'hs_faces': [jf(f) for f in faces], 'pentagons': dual,
            'counts': {'coarse_faces': len(coarse), 'squares': sum(aa), 'triangles': 2*sum(aa),
                       'hs_vertices': len(corners), 'hs_edges': len(edges), 'C': counts['C'], 'P': counts['P']},
            'exact_checks': ['regular_unit_hs_faces', 'opposite_edge_pairs_mod_period', 'connected_face_adjacency',
                             'single_2S3T_vertex_links', 'hs_area_equals_period_area', 'pentagon_area_equals_(2+sqrt3)/4',
                             'pentagon_incircular_radius_1/2', 'pentagon_area_sum_equals_period_area'],
            'scope': 'Exact specified constructions; no full-group search. Primitive periods/full groups use the analytic word theorem.'}


def main():
    directory = Path(__file__).resolve().parent
    for name in ('1123', '1112', '1122', '2'):
        result = build(tuple(map(int, name)))
        target = directory/f'word-{name}-round36.json'
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        print(name, result['counts'], 'periods=', result['periods'])


if __name__ == '__main__': main()
