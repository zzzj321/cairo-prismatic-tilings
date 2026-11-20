"""Exact physical Euclidean symmetries of eight specified HS face lists.

Independent of the array/word symmetry predicates: only periods and faces are
read. Python Fraction fallback because cached/live Julia is unavailable.
One finite primary search, not an unrestricted classification certificate.
"""
from collections import defaultdict, Counter
from fractions import Fraction as F
from pathlib import Path
from math import gcd, lcm
import json, time, hashlib, sys
from develop_words_round36 import K, point, add, sub, mul, dot, cross, center, pkey, fkey, jp

ZERO=point(0,0); SQRT3=K(0,1)
I=(K(1),K(),K(),K(1))
J=(K(1),K(),K(),K(-1))
R30=(SQRT3/2,K(F(-1,2)),K(F(1,2)),SQRT3/2)

def apply(M,p):
    a,b,c,d=M; x,y=p
    return a*x+b*y,c*x+d*y

def compose(M,N):
    a,b,c,d=M; e,f,g,h=N
    return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h

def decode(p):
    a,b,c,d=p
    return K(F(a,2),F(b,2)),K(F(c,2),F(d,2))

def egcd(a,b):
    if b==0: return abs(a), 1 if a>=0 else -1, 0
    g,x,y=egcd(b,a%b)
    return g,y,x-(a//b)*y

def integer_pair(a,b,c):
    """Solve a*m+b*n=c in Q(sqrt3), m,n integers, or return None."""
    det=a.a*b.b-a.b*b.a
    if det:
        m=(c.a*b.b-c.b*b.a)/det
        n=(a.a*c.b-a.b*c.a)/det
        return (int(m),int(n)) if m.denominator==n.denominator==1 else None
    if a.a*c.b!=a.b*c.a or b.a*c.b!=b.b*c.a: return None
    if a==K() and b==K(): return (0,0) if c==K() else None
    aa,bb,cc=(a.a,b.a,c.a) if a.a or b.a else (a.b,b.b,c.b)
    den=lcm(aa.denominator,bb.denominator,cc.denominator)
    aa,bb,cc=int(aa*den),int(bb*den),int(cc*den)
    g,m,n=egcd(aa,bb)
    if cc%g: return None
    m,n=m*(cc//g),n*(cc//g)
    assert a*m+b*n==c
    return m,n

def analyze(path,expected):
    started=time.monotonic()
    data=json.loads(path.read_text(encoding='utf-8'))
    U,V=map(decode,data['periods']); det=cross(U,V)
    assert det!=K()
    def red(p):
        m=(cross(p,V)/det).floor(); n=(cross(U,p)/det).floor()
        return sub(p,add(mul(U,m),mul(V,n)))
    def norm(f,c=None):
        c=center(f) if c is None else c
        shift=sub(c,red(c))
        return fkey(tuple(sub(p,shift) for p in f))
    faces=[tuple(map(decode,f)) for f in data['faces']]
    keys={norm(f) for f in faces}
    assert len(keys)==len(faces)
    vertices=defaultdict(list)
    for f in faces:
        for j,p in enumerate(f):
            nxt=sub(f[(j+1)%len(f)],p); prev=sub(f[j-1],p)
            assert dot(nxt,nxt)==K(1)
            vertices[red(p)].append((nxt,prev,len(f)))
    stars={}; cpoints=[]; types=Counter()
    for p,entries in vertices.items():
        assert Counter(e[2] for e in entries)=={3:3,4:2}
        links={e[0]:(e[1],e[2]) for e in entries}
        assert len(links)==5
        ray=next(iter(links)); start=ray; visited=set(); sizes=[]
        for _ in range(5):
            assert ray not in visited
            visited.add(ray); ray,size=links[ray]; sizes.append(size)
        assert ray==start
        typ='P' if any(sizes[i]==sizes[(i+1)%5]==4 for i in range(5)) else 'C'
        types[typ]+=1
        stars[p]=frozenset(links)
        if typ=='C': cpoints.append(p)
    assert cpoints
    origin=min(cpoints,key=pkey); star0=stars[origin]
    vertex_set=set(vertices)
    stats={'linear_candidates':24,'star_phase_candidates':0,'face_tests':0}
    transformations=[]
    R=I
    for k in range(12):
        transformations.extend([(k,False,R),(k,True,compose(R,J))])
        R=compose(R30,R)
    assert R==I and len({M for k,s,M in transformations})==24

    def phase_search(M):
        star=frozenset(apply(M,e) for e in star0)
        mo=apply(M,origin)
        phases={red(sub(q,mo)) for q in cpoints if stars[q]==star}
        stats['star_phase_candidates']+=len(phases)
        tf=[(tuple(apply(M,p) for p in f),apply(M,center(f))) for f in faces]
        tv=[apply(M,p) for p in vertices]
        good=[]
        for t in sorted(phases,key=pkey):
            # Cheap necessary filter before the complete exact face test.
            if any(red(add(p,t)) not in vertex_set for p in tv): continue
            stats['face_tests']+=1
            if all(norm(tuple(add(p,t) for p in f),add(c,t)) in keys for f,c in tf):
                good.append(t)
        return good

    translations=phase_search(I)
    translation_set=set(translations)
    assert ZERO in translation_set
    def istranslation(p): return red(p) in translation_set
    # Cheap group invariant: exact closure of discovered finite translations.
    assert all(istranslation(add(a,b)) for a in translations for b in translations)
    elements=[]; linears={}; mirror_directions=set()
    for k,reversing,M in transformations:
        # The displayed period lattice need not be preserved by M: its images
        # need only be translations of the full physical tiling.
        if not all(istranslation(apply(M,p)) for p in (U,V)): continue
        phases=translations if M==I else phase_search(M)
        if phases: linears[(k,reversing)]=M
        for t in phases:
            row={'rotation_step_30deg':k,'orientation_reversing':reversing,
                 'linear_matrix':[z.json() for z in M], 'phase_mod_input_lattice':jp(t)}
            if reversing:
                # An eigenvector with eigenvalue +1, obtained without trig.
                axis=add(point(1,0),apply(M,point(1,0)))
                if axis==ZERO: axis=add(point(0,1),apply(M,point(0,1)))
                assert axis!=ZERO and apply(M,axis)==axis
                witness=integer_pair(dot(axis,U),dot(axis,V),-dot(axis,t))
                row['true_mirror_in_this_input_lattice_coset']=witness is not None
                row['glide_square_translation']=jp(add(t,apply(M,t)))
                if witness is not None:
                    m,n=witness; tt=add(t,add(mul(U,m),mul(V,n)))
                    assert add(tt,apply(M,tt))==ZERO
                    row['mirror_lattice_shift']=[m,n]
                    row['mirror_phase']=jp(tt)
                    mirror_directions.add(k)
            elements.append(row)
    point_order=len(linears)
    nrot=sum(not b for k,b in linears)
    refs=[(k,M) for (k,b),M in linears.items() if b]
    centered=None
    if not refs:
        group={1:'p1',2:'p2',4:'p4'}.get(point_order,'unclassified_rotation_group')
    elif point_order==8:
        group='p4m' if len(mirror_directions)==4 else 'p4g' if len(mirror_directions)==2 else 'unclassified_D4'
    elif not mirror_directions:
        group={2:'pg',4:'pgg'}.get(point_order,'unclassified_glide_group')
    else:
        # Gamma splits into its reflection-axis/perpendicular sublattices iff
        # (R-I)Gamma/2 is contained in Gamma. U,V and all phases generate Gamma.
        M=refs[0][1]
        rectangular=all(istranslation(mul(sub(apply(M,g),g),F(1,2)))
                        for g in [U,V]+translations)
        centered=not rectangular
        if point_order==2: group='cm' if centered else 'pm'
        elif point_order==4:
            if len(mirror_directions)==1: group='pmg'
            elif len(mirror_directions)==2: group='cmm' if centered else 'pmm'
            else: group='unclassified_D2'
        else: group='unclassified_reflection_group'
    assert all(sum(e['rotation_step_30deg']==k and e['orientation_reversing']==b for e in elements)==len(translations)
               for k,b in linears), 'Every present linear coset must have |Gamma/Lambda| phases'
    result={'input':path.name,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'expected':expected,'physical_group':group,'matches_expected':group==expected,
            'input_periods':[jp(U),jp(V)],'unit_faces':dict(Counter(map(len,faces))),
            'vertex_types':dict(types),'translation_coset_count':len(translations),
            'translation_phases':[jp(t) for t in translations],
            'full_translation_covolume':((det if det.sign()>0 else -det)/len(translations)).json(),
            'point_group_order':point_order,'rotation_subgroup_order':nrot,
            'true_mirror_linear_steps':sorted(mirror_directions),'centered_reflection_lattice':centered,
            'elements_mod_input_lattice':elements,'search_statistics':stats,
            'elapsed_seconds':time.monotonic()-started}
    return result

def main():
    base=Path(__file__).resolve().parent
    cases=[('matrix-p4','p4'),('matrix-p4g','p4g'),('matrix-cmm','cmm'),('matrix-cm','cm'),
           ('word-1123','p1'),('word-1112','p2'),('word-1122','pg'),('word-2','pgg')]
    results=[]; started=time.monotonic()
    out=base/'physical-groups-round36.json'
    for name,expected in cases:
        result=analyze(base/f'{name}-round36.json',expected)
        results.append(result)
        print(name,result['physical_group'],'expected',expected,'translations',result['translation_coset_count'],
              'point',result['point_group_order'],'seconds',round(result['elapsed_seconds'],3),flush=True)
        payload={'evidence':'exact-check','scope':'Complete physical isometries of eight specified periodic HS constructions; not unrestricted classification.',
                 'backend':sys.version,'arithmetic':'fractions.Fraction in Q(sqrt(3)); no floating decisions',
                 'algorithm':'24 linear parts; a C vertex image fixes phase; exact rays, vertices, faces and image-period translations; integer mirror-phase equation',
                 'results':results,'elapsed_seconds':time.monotonic()-started}
        out.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    assert all(r['matches_expected'] for r in results), 'Physical output differs from supplied expectation: inspect saved exact elements.'

if __name__=='__main__': main()
