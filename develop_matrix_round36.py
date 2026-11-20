"""Exact developing map for FOUR specified regular binary matrix models.

Python stdlib fallback: Julia unavailable in the retained task environment.
No symmetry classification and no unbounded classification claim.
Points encode as [a,b,c,d] = ((a+b*sqrt(3))/2,(c+d*sqrt(3))/2).
Run directly to write develop-matrix-round36.json next to this module.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from collections import deque, Counter
from pathlib import Path
import json
from math import isqrt

@dataclass(frozen=True)
class Q:
    a: F = F(0)
    b: F = F(0)
    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))
    def __add__(self, z):
        z = z if isinstance(z,Q) else Q(z)
        return Q(self.a+z.a,self.b+z.b)
    __radd__=__add__
    def __neg__(self): return Q(-self.a,-self.b)
    def __sub__(self,z): return self+-asq(z)
    def __rsub__(self,z): return asq(z)+-self
    def __mul__(self,z):
        z=asq(z)
        return Q(self.a*z.a+3*self.b*z.b,self.a*z.b+self.b*z.a)
    __rmul__=__mul__
    def __truediv__(self,z):
        z=asq(z); norm=z.a*z.a-3*z.b*z.b
        assert norm
        return self*Q(z.a/norm,-z.b/norm)
    def __eq__(self,z):
        z=asq(z); return self.a==z.a and self.b==z.b
    def sign(self):
        if not self.b: return (self.a>0)-(self.a<0)
        if not self.a: return (self.b>0)-(self.b<0)
        if (self.a>0)==(self.b>0): return 1 if self.a>0 else -1
        delta=self.a*self.a-3*self.b*self.b
        return ((delta>0)-(delta<0))*(1 if self.a>0 else -1)
    def json(self): return [self.a.numerator,self.a.denominator,self.b.numerator,self.b.denominator]

def asq(z): return z if isinstance(z,Q) else Q(z)
def add(p,q): return (p[0]+q[0],p[1]+q[1])
def sub(p,q): return (p[0]-q[0],p[1]-q[1])
def scale(a,p): return (a*p[0],a*p[1])
def cross(p,q): return p[0]*q[1]-p[1]*q[0]
def dot(p,q): return p[0]*q[0]+p[1]*q[1]
O=(Q(),Q())
DIRECTIONS=[(Q(1),Q()),(Q(0,F(1,2)),Q(F(1,2))),
 (Q(F(1,2)),Q(0,F(1,2))),(Q(),Q(1)),
 (Q(-F(1,2)),Q(0,F(1,2))),(Q(0,-F(1,2)),Q(F(1,2)))]
DIRECTIONS += [scale(-1,p) for p in DIRECTIONS]

def encode(p):
    vals=[2*p[0].a,2*p[0].b,2*p[1].a,2*p[1].b]
    assert all(z.denominator==1 for z in vals), vals
    return [int(z) for z in vals]
def decode(p): return (Q(F(p[0],2),F(p[1],2)),Q(F(p[2],2),F(p[3],2)))
def area(poly):
    return sum((cross(p,q) for p,q in zip(poly,poly[1:]+poly[:1])),Q())/2
def keyface(face): return tuple(sorted(tuple(encode(p)) for p in face))

EXAMPLES={
 'p4': [[int(q==(2*p)%5) for q in range(5)] for p in range(5)],
 'p4g': [[1,0],[0,1]],
 'cmm': [[1,1,0,0],[0,0,1,1]],
 'cm': [[1,1,0,1,0,0],[1,1,0,0,0,1],[0,1,1,0,0,1],
        [1,0,0,1,1,0],[0,0,1,1,1,0],[0,0,1,0,1,1]],
}

def build_model(matrix):
    r,n=len(matrix),len(matrix[0]); I,J=2*r,2*n
    assert all(len(row)==n and set(row)<={0,1} for row in matrix)
    assert len({sum(row) for row in matrix})==1
    assert len({sum(matrix[p][q] for p in range(r)) for q in range(n)})==1
    eta=lambda p,q:2*matrix[p%r][q%n]-1
    q={}; dh={}; dv={}
    for i in range(I):
        for j in range(J):
            if (i+j)%2:
                p,s=i//2,j//2
                if i%2==0: z=(-eta(p-1,s),eta(p,s),-eta(p,s),eta(p-1,s))
                else: z=(eta(p,s-1),-eta(p,s-1),eta(p,s),-eta(p,s))
                assert sum(z)==0
                q[i,j]=z; dh[i,j]=(z[2]+z[3])//2; dv[i,j]=(z[1]+z[2])//2
            else: dh[i,j]=dv[i,j]=0
    H={}; V={}
    for i in range(I):
        pre=[0]
        for j in range(J): pre.append(pre[-1]+dh[i,j])
        assert pre[-1]==0
        for j in range(J): H[i,j]=1+pre[j]-min(pre)
    for j in range(J):
        pre=[0]
        for i in range(I): pre.append(pre[-1]+dv[i,j])
        assert pre[-1]==0
        for i in range(I): V[i,j]=1+pre[i]-min(pre)
    for i in range(I):
        for j in range(J):
            hs=H[i,j],H[i,(j+1)%J]; vs=V[i,j],V[(i+1)%I,j]
            if dh[i,j]: assert vs==(1,1)
            elif dv[i,j]: assert hs==(1,1)
            else: assert min(hs[0],vs[0])==1
    return {'I':I,'J':J,'q':q,'H':H,'V':V,'matrix':matrix}

def face_data(model,i,j):
    I,J=model['I'],model['J']; i%=I; j%=J
    H,V=model['H'],model['V']
    lengths=(H[i,j],V[(i+1)%I,j],H[i,(j+1)%J],V[i,j])
    angles=(3,3,3,3) if (i+j)%2==0 else tuple(3-z for z in model['q'][i,j])
    offsets=[0]
    for k in (1,2,3): offsets.append(offsets[-1]+6-angles[k])
    assert offsets[-1]+6-angles[0]==12
    return lengths,angles,offsets

def develop(model):
    """BFS on a fundamental rectangle PLUS one row/column collar.

    Every encounter with an already developed edge/vertex must agree exactly;
    the collar checks both fundamental periods and their corner composition.
    """
    I,J=model['I'],model['J']
    domain={(i,j) for i in range(I+1) for j in range(J+1)}
    states={(0,0):(O,0)}; vertices={}; faces={}; queue=deque([(0,0)])
    steps=((0,-1),(1,0),(0,1),(-1,0)); cs=((0,0),(1,0),(1,1),(0,1))
    checks=0
    while queue:
        i,j=queue.popleft(); origin,theta=states[i,j]
        lens,angles,off=face_data(model,i,j)
        vecs=[scale(L,DIRECTIONS[(theta+d)%12]) for L,d in zip(lens,off)]
        poly=[origin]
        for v in vecs: poly.append(add(poly[-1],v))
        assert poly[-1]==origin, ('face nonclosure',i,j)
        poly=poly[:-1]; assert area(poly).sign()>0
        for k in range(4):
            e=sub(poly[(k+1)%4],poly[k]); f=sub(poly[k-1],poly[k])
            assert dot(e,e)==lens[k]**2
            cos={2:Q(F(1,2)),3:Q(),4:Q(-F(1,2))}[angles[k]]
            assert dot(e,f)==lens[k]*lens[k-1]*cos
            ij=i+cs[k][0],j+cs[k][1]
            if ij in vertices: assert vertices[ij]==poly[k], ('vertex monodromy',ij)
            else: vertices[ij]=poly[k]
        faces[i,j]=poly
        for k,(di,dj) in enumerate(steps):
            ij=i+di,j+dj
            if ij not in domain: continue
            nk=(k+2)%4
            nl,na,no=face_data(model,*ij)
            nt=(theta+off[k]+6-no[nk])%12
            shift=O
            for t in range(nk): shift=add(shift,scale(nl[t],DIRECTIONS[(nt+no[t])%12]))
            no_origin=sub(poly[(k+1)%4],shift)
            candidate=(no_origin,nt)
            if ij in states:
                assert states[ij]==candidate, ('path-dependent face',ij)
                checks+=1
            else: states[ij]=candidate; queue.append(ij)
    U=sub(vertices[I,0],vertices[0,0]); W=sub(vertices[0,J],vertices[0,0])
    assert cross(U,W).sign()>0
    for i,j in list(vertices):
        if (i+I,j) in vertices: assert sub(vertices[i+I,j],vertices[i,j])==U
        if (i,j+J) in vertices: assert sub(vertices[i,j+J],vertices[i,j])==W
    for i,j in states:
        if (i+I,j) in states: assert states[i+I,j][1]==states[i,j][1]
        if (i,j+J) in states: assert states[i,j+J][1]==states[i,j][1]
    fundamental={(i,j):faces[i,j] for i in range(I) for j in range(J)}
    assert sum((area(z) for z in fundamental.values()),Q())==cross(U,W)
    return {'model':model,'vertices':vertices,'faces':fundamental,'periods':(U,W),
            'directions':states,'path_checks':checks}

def integer(z):
    assert z.b==0 and z.a.denominator==1
    return int(z.a)

def edge_length(p,q):
    squared=integer(dot(sub(q,p),sub(q,p)))
    length=isqrt(squared)
    assert length*length==squared and length>0
    return length

def refine_face(poly,is_square):
    """Tile one convex coarse polygon in its exact square/triangle lattice."""
    L=edge_length(poly[0],poly[1])
    e=scale(F(1,L),sub(poly[1],poly[0]))
    if is_square: f=(-e[1],e[0])
    else: f=(e[0]/2-Q(0,F(1,2))*e[1],Q(0,F(1,2))*e[0]+e[1]/2)
    det=cross(e,f)
    lp=[]
    for p in poly:
        v=sub(p,poly[0]); lp.append((integer(cross(v,f)/det),integer(cross(e,v)/det)))
    def inside(p):
        return all((b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])>=0 for a,b in zip(lp,lp[1:]+lp[:1]))
    cells=[]
    for a in range(min(x for x,y in lp),max(x for x,y in lp)+1):
        for b in range(min(y for x,y in lp),max(y for x,y in lp)+1):
            if is_square: candidates=[[(a,b),(a+1,b),(a+1,b+1),(a,b+1)]]
            else: candidates=[[(a,b),(a+1,b),(a,b+1)],[(a+1,b),(a+1,b+1),(a,b+1)]]
            for local in candidates:
                if all(inside(p) for p in local):
                    cells.append([add(poly[0],add(scale(x,e),scale(y,f))) for x,y in local])
    edges=Counter()
    for cell in cells:
        assert area(cell)==(Q(1) if is_square else Q(0,F(1,4)))
        for p,q in zip(cell,cell[1:]+cell[:1]):
            assert dot(sub(q,p),sub(q,p))==1
            edges[tuple(encode(p)),tuple(encode(q))]+=1
    boundary=Counter()
    for p,q in zip(poly,poly[1:]+poly[:1]):
        length=edge_length(p,q)
        for k in range(length):
            p1=add(p,scale(F(k,length),sub(q,p))); p2=add(p,scale(F(k+1,length),sub(q,p)))
            boundary[tuple(encode(p1)),tuple(encode(p2))]+=1
    for edge,n in edges.items():
        assert n==1
        assert n-edges[edge[::-1]]==boundary[edge]-boundary[edge[::-1]], ('boundary mismatch',edge)
    assert sum((area(z) for z in cells),Q())==area(poly)
    assert all(edges[e]==1 for e in boundary)
    return cells

def refine(developed):
    out=[]
    for ij,poly in developed['faces'].items():
        typ='S' if sum(ij)%2==0 else 'T'
        for cell in refine_face(poly,typ=='S'): out.append({'type':typ,'coarse':ij,'vertices':cell})
    assert len({keyface(z['vertices']) for z in out})==len(out)
    return out

def run_example(name,matrix):
    d=develop(build_model(matrix)); fine=refine(d); I,J=d['model']['I'],d['model']['J']
    counts=Counter(z['type'] for z in fine)
    assert counts['T']==2*counts['S']
    return {'name':name,'matrix':matrix,'coarse_periods':[I,J],
      'coordinate_encoding':'[a,b,c,d] means ((a+b sqrt(3))/2,(c+d sqrt(3))/2)',
      'period_vectors':[encode(p) for p in d['periods']],
      'cell_area':cross(*d['periods']).json(),
      'coarse_vertices':[[i,j,*encode(p)] for (i,j),p in sorted(d['vertices'].items()) if i<=I and j<=J],
      'coarse_faces':[{'index':list(ij),'type':'S' if sum(ij)%2==0 else 'T',
                       'vertices':[encode(p) for p in poly]} for ij,poly in d['faces'].items()],
      'fine_faces':[{'type':z['type'],'coarse':list(z['coarse']),'vertices':[encode(p) for p in z['vertices']]} for z in fine],
      'checks':{'path_consistency_encounters':d['path_checks'],'period_and_corner_closure':True,
                'coarse_lengths_angles':True,'fine_unit_edges_area':True,'refinement_oriented_boundary':True,
                'coarse_area_equals_period_determinant':True},
      'display_cell_counts':dict(counts)}

if __name__=='__main__':
    data={'evidence':'exact-check','scope':'four specified matrix inputs only; no full-group claim',
          'examples':{name:run_example(name,x) for name,x in EXAMPLES.items()}}
    target=Path(__file__).resolve().with_name('develop-matrix-round36.json')
    target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    for name,example in data['examples'].items():
        payload={**example,'periods':example['period_vectors'],
                 'faces':[face['vertices'] for face in example['fine_faces']]}
        target.with_name('matrix-'+name+'-round36.json').write_text(
            json.dumps(payload,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:{'periods':v['period_vectors'],'counts':v['display_cell_counts'],'checks':v['checks']}
                      for k,v in data['examples'].items()},ensure_ascii=False))
