import json,pickle,sys,numpy as np,random,collections
sys.path.insert(0,'/home/sandbox/work/repos/mega27-20-drug-target-prediction/src/targetscan')
from collision_floor import collision_mse_floor
L=list(json.load(open('ligands_can.txt')).values());P=list(json.load(open('proteins.txt')).values())
Y=pickle.load(open('Y','rb'),encoding='latin1');pk=-np.log10(Y/1e9)
tr=json.load(open('train_fold_setting1.txt'));te=json.load(open('test_fold_setting1.txt'))
print('train structure',type(tr),len(tr),[len(x) for x in tr],'test',len(te))
pool=sorted(set(i for f in tr for i in f));print('pool',len(pool),'overlap pool/test',len(set(pool)&set(te)))
N=68*442
li=[i//442 for i in range(N)];pj=[i%442 for i in range(N)]
y=np.array([pk[li[i],pj[i]] for i in range(N)])
key=lambda i:(L[li[i]]+'\x00'+P[pj[i]])
def h1(trn,tst):
    d=collections.defaultdict(list)
    for i in trn:d[key(i)].append(y[i])
    hit=[i for i in tst if key(i) in d]
    se=[(y[i]-np.mean(d[key(i)]))**2 for i in hit]
    return len(hit),(float(np.mean(se)) if se else None)
print('H1 fold',h1(pool,te))
print('dup protein entries',442-len(set(P)))
rs=[];
for s in range(100):
    r=random.Random(s);idx=list(range(N));r.shuffle(idx);t=idx[:len(te)];t0=idx[len(te):]
    rs.append(h1(t0,t)[0])
print('H1 random 100 splits hit count min/med/max',min(rs),np.median(rs),max(rs))
# H3
ytr=y[pool];mu=ytr.mean()
a=np.zeros(68);b=np.zeros(442);li_a=np.array([li[i] for i in pool]);pj_a=np.array([pj[i] for i in pool])
lam=1.0
for _ in range(20):
    r=ytr-mu-b[pj_a]
    for k in range(68):
        m=li_a==k;a[k]=r[m].sum()/(m.sum()+lam) if m.any() else 0
    r=ytr-mu-a[li_a]
    for k in range(442):
        m=pj_a==k;b[k]=r[m].sum()/(m.sum()+lam) if m.any() else 0
seenL=set(li_a);seenP=set(pj_a)
def pred(i,use_a=True,use_b=True):
    return mu+(a[li[i]] if use_a and li[i] in seenL else 0)+(b[pj[i]] if use_b and pj[i] in seenP else 0)
yt=y[te]
def ci(yt,p):
    yt=np.asarray(yt);p=np.asarray(p);s=0;n=0
    o=np.argsort(yt)
    for ii in range(len(o)):
        d=yt[o[ii+1:]]>yt[o[ii]];s+=((p[o[ii+1:]][d]>p[o[ii]])+0.5*(p[o[ii+1:]][d]==p[o[ii]])).sum();n+=d.sum()
    return s/n
for nm,fa,fb in [('additive',1,1),('ligand-only',1,0),('protein-only',0,1)]:
    p=[pred(i,fa,fb) for i in te];print('H3',nm,'MSE',float(np.mean((yt-p)**2)),'CI',ci(yt,p))
print('H3 global-mean MSE',float(np.mean((yt-mu)**2)),'CI 0.5 by construction')
# H2
mk=lambda ids,f:[f(i).encode() for i in ids]
full=lambda i:key(i)
def fl(ids,f):return collision_mse_floor(mk(ids,f),[float(y[i]) for i in ids])['empirical_mse_floor']
print('H2 floor pool',fl(pool,full),'test',fl(te,full))
print('H2 ligand-only key',fl(range(N),lambda i:L[li[i]]),'protein-seq-only key',fl(range(N),lambda i:P[pj[i]]))
sh=[]
for s in range(20):
    r=np.random.RandomState(s);yy=r.permutation(y)
    sh.append(collision_mse_floor(mk(range(N),full),[float(v) for v in yy])['empirical_mse_floor'])
print('H2 shuffled median',float(np.median(sh)),'label var',float(y.var()))
# H5
g=collections.defaultdict(list)
for j in range(442):g[P[j]].append(j)
dis=[]
for s,js in g.items():
    if len(js)>1:
        for i in range(68):
            v=[pk[i,j] for j in js]
            if max(v)!=min(v):dis.append(max(v)-min(v))
print('H5 dup-seq groups',sum(len(v)>1 for v in g.values()),'ligand-groups differing',len(dis),'median diff',float(np.median(dis)) if dis else None,'max',max(dis) if dis else None)
