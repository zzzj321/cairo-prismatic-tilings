"""Exact HS-to-C/P reconstruction for eight specified input face lists.

Python stdlib fallback: Julia unavailable in the current backend inventory.
No floats in the reconstruction or checks. Reuses the Q(sqrt(3)) field class,
but does not use the optional precomputed pentagons in word input files.
"""
from collections import Counter, defaultdict
from pathlib import Path
from fractions import Fraction as F
import json, time
from develop_words_round36 import K, S, add, sub, mul, dot, cross, center, area, pkey, fkey, jp, jf

HERE = Path(__file__).resolve().parent
CASES = [('p1','word-1123'),('p2','word-1112'),('pg','word-1122'),('pgg','word-2'),
         ('p4','matrix-p4'),('p4g','matrix-p4g'),('cmm','matrix-cmm'),('cm','matrix-cm')]
EXPECTED = {'p1':(8,6),'p2':(8,2),'pg':(8,4),'pgg':(4,4),
            'p4':(100,40),'p4g':(16,8),'cmm':(32,20),'cm':(144,72)}

def decode(p): return K(F(p[0],2),F(p[1],2)), K(F(p[2],2),F(p[3],2))

def build(group, input_path):
    data=json.loads(input_path.read_text(encoding='utf-8'))
    U,V=map(decode,data['periods'])
    determinant=cross(U,V)
    assert determinant.sign()>0
    def reduce(p):
        return sub(p,add(mul(U,(cross(p,V)/determinant).floor()),
                         mul(V,(cross(U,p)/determinant).floor())))
    def norm(f):
        c=center(f); t=sub(c,reduce(c))
        return tuple(sub(p,t) for p in f)
    faces=[norm(tuple(map(decode,f))) for f in data['faces']]
    assert len({fkey(f) for f in faces})==len(faces)
    edges,corners=defaultdict(list),defaultdict(list)
    for index,f in enumerate(faces):
        assert len(f) in (3,4)
        assert area(f)==(K(1) if len(f)==4 else S/4)
        for j,p in enumerate(f):
            nxt,prev=f[(j+1)%len(f)],f[j-1]
            a,b=sub(nxt,p),sub(prev,p)
            assert dot(a,a)==K(1)
            assert dot(a,b)==(K() if len(f)==4 else K(F(1,2)))
            assert cross(a,b).sign()>0
            midpoint=mul(add(p,nxt),F(1,2)); shift=sub(midpoint,reduce(midpoint))
            ep=(sub(p,shift),sub(nxt,shift))
            edges[fkey(ep)].append((*ep,index))
            corners[reduce(p)].append((a,b,len(f),sub(center(f),p),index))
    adjacent=defaultdict(set)
    for es in edges.values():
        assert len(es)==2 and es[0][:2]==es[1][:2][::-1]
        a,b=es[0][2],es[1][2]; adjacent[a].add(b); adjacent[b].add(a)
    reached={0}; pending=[0]
    while pending:
        i=pending.pop()
        for j in adjacent[i]-reached: reached.add(j); pending.append(j)
    assert len(reached)==len(faces)
    assert sum((area(f) for f in faces),K())==determinant
    pentagons=[]; counts=Counter(); dualedges=defaultdict(list)
    for vertex,cs in corners.items():
        assert Counter(c[2] for c in cs)=={3:3,4:2}
        links={c[0]:c for c in cs}; assert len(links)==5
        start=min(links,key=pkey); ray=start; seen=set(); polygon=[]; sizes=[]
        for _ in range(5):
            assert ray not in seen; seen.add(ray)
            c=links[ray]; sizes.append(c[2]); polygon.append(add(vertex,c[3])); ray=c[1]
        assert ray==start and len(seen)==5
        typ='P' if any(sizes[j]==sizes[(j+1)%5]==4 for j in range(5)) else 'C'
        assert area(polygon)==(2+S)/4
        lengths=Counter()
        for j,p in enumerate(polygon):
            nxt=polygon[(j+1)%5]; e=sub(nxt,p); length2=dot(e,e)
            distance_numerator=cross(e,sub(vertex,p))
            assert distance_numerator.sign()>0
            assert 4*distance_numerator*distance_numerator==length2
            lengths[length2]+=1
            midpoint=mul(add(p,nxt),F(1,2)); shift=sub(midpoint,reduce(midpoint))
            edge=(sub(p,shift),sub(nxt,shift)); dualedges[fkey(edge)].append(edge)
        a2=K(F(1,3)); b2=(2+S)/6; c2=K(1)
        assert lengths==({a2:2,b2:2,c2:1} if typ=='P' else {a2:1,b2:4})
        counts[typ]+=1
        pentagons.append({'type':typ,'incenter':jp(vertex),'vertices':jf(polygon)})
    for es in dualedges.values(): assert len(es)==2 and es[0]==es[1][::-1]
    assert (counts['C'],counts['P'])==EXPECTED[group]
    assert len(corners)*(2+S)/4==determinant
    assert len(corners)==2*sum(len(f)==4 for f in faces)==sum(len(f)==3 for f in faces)
    return {'name':group,'source':input_path.name,'periods':[jp(U),jp(V)],
            'coordinate_encoding':'Each scalar [a,b] is Fraction(a)+Fraction(b)*sqrt(3).',
            'period_area':determinant.json(),'pentagons':pentagons,
            'counts':{'S':sum(len(f)==4 for f in faces),'T':sum(len(f)==3 for f in faces),
                      'C':counts['C'],'P':counts['P'],'HS_edges':len(edges),'CP_edges':len(dualedges)},
            'checks':['exact_regular_unit_faces','opposite_HS_edges','connected_HS_faces',
                      'single_2S3T_links','HS_area_equals_period_determinant',
                      'exact_pentagon_area_and_side_multisets','inradius_one_half',
                      'opposite_CP_edges','CP_area_equals_period_determinant'],
            'scope':'Eight specified coordinate constructions; group names identify intended inputs, not results of this checker.'}

def main():
    started=time.perf_counter(); results=[]
    for group,stem in CASES:
        result=build(group,HERE/f'{stem}-round36.json')
        (HERE/f'cp-{group}-round36.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        results.append({'name':group,'counts':result['counts'],'checks':len(result['checks'])})
    summary={'results':results,'elapsed_seconds':time.perf_counter()-started,'evidence':'exact-check'}
    (HERE/'dualization-round36-results.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False))

if __name__=='__main__': main()
