import re,json,pickle,collections,numpy as np,random
D='/tmp/davis-src/'
h=open(D+'kinhub.html').read()
rows=re.findall(r'<tr>\s*((?:<td>.*?</td>\s*){8})</tr>',h,re.S)
m={}
for r in rows:
    c=[re.sub(r'<.*?>','',x).strip() for x in re.findall(r'<td>(.*?)</td>',r,re.S)]
    for k in (c[0],c[1],c[2]): m.setdefault(k.upper(),c)
L=list(json.load(open(D+'ligands_can.txt')).values());Pd=json.load(open(D+'proteins.txt'));names=list(Pd);P=list(Pd.values())
Y=pickle.load(open(D+'Y','rb'),encoding='latin1');y=-np.log10(Y/1e9)
def base(n):
    n=re.sub(r'\(.*?\)','',n);n=re.sub(r'-(phosphorylated|nonphosphorylated|autoinhibited|cyclic).*','',n,flags=re.I)
    if n.upper() not in m and n.endswith('p') and n[:-1].upper() in m:n=n[:-1]
    return n.upper()
fam={}
for j,n in enumerate(names):
    b=base(n)
    if b in m: fam[j]=m[b][5]
print('mapped',len(fam),'families',len(set(fam.values())),'unmapped',442-len(fam))
cnt=collections.Counter(fam.values())
folds=[[] for _ in range(5)];load=[0]*5
for f,c in sorted(cnt.items(),key=lambda x:(-x[1],x[0])):
    k=load.index(min(load));folds[k].append(f);load[k]+=c
print('proteins per fold',load)
# 3-mer
idx={};
def kmers(s):
    v=collections.Counter(s[i:i+3] for i in range(len(s)-2));return v
allk=sorted(set(k for s in P for k in kmers(s)));ki={k:i for i,k in enumerate(allk)}
X=np.zeros((442,len(allk)),dtype=np.float32)
for j,s in enumerate(P):
    for k,c in kmers(s).items():X[j,ki[k]]=c
X/=np.linalg.norm(X,axis=1,keepdims=True)
S=X@X.T
def ci(yt,p):
    yt=np.asarray(yt);p=np.asarray(p);o=np.argsort(yt);s=0;n=0
    for i in range(len(o)):
        d=yt[o[i+1:]]>yt[o[i]];q=p[o[i+1:]][d];s+=((q>p[o[i]])+0.5*(q==p[o[i]])).sum();n+=d.sum()
    return s/n
res=[];tp=[]
for k in range(5):
    tf=set(folds[k]);test=[j for j,f in fam.items() if f in tf];tr=[j for j in range(442) if j not in set(test)]
    mu=y[:,tr].mean();lm=y[:,tr].mean(axis=1)
    for j in test:
        sim=S[j,tr];nn=[tr[i] for i in np.argsort(-sim)[:5]]
        p2=y[:,nn].mean(axis=1)
        for i in range(68):
            res.append((k,j,i,y[i,j],mu,lm[i],p2[i]))
R=np.array(res);yt=R[:,3]
for nm,c in [('M0 global mean',4),('M1 ligand mean',5),('M2 seq-kNN',6)]:
    pf=[float(np.mean((R[R[:,0]==k,3]-R[R[:,0]==k,c])**2)) for k in range(5)]
    print(nm,'pooled MSE',float(np.mean((yt-R[:,c])**2)),'CI',ci(yt,R[:,c]),'per-fold',[round(x,3) for x in pf])
rng=np.random.RandomState(0);prots=sorted(set(R[:,1].astype(int)));byp={p:R[R[:,1]==p] for p in prots}
dl=[]
for _ in range(2000):
    s=rng.choice(prots,len(prots));A=np.vstack([byp[p] for p in s])
    dl.append(np.mean((A[:,3]-A[:,6])**2)-np.mean((A[:,3]-A[:,5])**2))
print('M2-M1 MSE diff',float(np.mean((yt-R[:,6])**2)-np.mean((yt-R[:,5])**2)),'bootstrap 95% CI',np.percentile(dl,[2.5,97.5]).tolist(),'n scored rows',len(R),'scored proteins',len(prots))
