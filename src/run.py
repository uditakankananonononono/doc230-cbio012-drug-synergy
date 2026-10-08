import pandas as pd,numpy as np,json,itertools
from sklearn.ensemble import HistGradientBoostingClassifier as HGC,HistGradientBoostingRegressor as HGR
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold,GroupKFold
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import TruncatedSVD
l=pd.read_csv('/tmp/syn/labels.csv',index_col=0)
n=pd.read_csv('/tmp/het/nodes.tsv',sep='\t'); c=n[n.kind=='Compound']; nm={x.lower():i for x,i in zip(c.name,c.id)}
al={'5-FU':'fluorouracil','MITOMYCINE':'mitomycin','ZOLINZA':'vorinostat','SN-38':'irinotecan'}
E=pd.read_csv('/tmp/het/edges.sif.gz',sep='\t'); drugs=sorted(set(l.drug_a_name)|set(l.drug_b_name)); cells=sorted(l.cell_line.unique())
cb=E[E.metaedge=='CbG']; T={}
for d in drugs:
    k=al.get(d,d.lower()).lower()
    if k in nm: T[d]=set(cb[cb.source==nm[k]].target)
mapped=[d for d in drugs if d in T and len(T[d])>0]
gg=E[E.metaedge.isin(['GiG','Gr>G'])]; adj={}
for a,b in zip(gg.source,gg.target): adj.setdefault(a,set()).add(b); adj.setdefault(b,set()).add(a)
pw=E[E.metaedge=='GpPW']; g2p={}
for a,b in zip(pw.source,pw.target): g2p.setdefault(a,set()).add(b)
P={d:set().union(*[g2p.get(g,set()) for g in T[d]]) for d in mapped}
allg=sorted(set().union(*T.values())); gi={g:i for i,g in enumerate(allg)}; V=np.zeros((len(mapped),len(allg)))
for i,d in enumerate(mapped):
    for g in T[d]: V[i,gi[g]]=1
Z=TruncatedSVD(8,random_state=0).fit_transform(V); Zd={d:Z[i] for i,d in enumerate(mapped)}
def links(a,b): return sum(len(adj.get(g,set())&T[b]) for g in T[a])
def gfe(a,b):
    ta,tb=T[a],T[b]; pa,pb=P[a],P[b]
    ov=[len(ta)+len(tb),len(ta)*len(tb),len(ta&tb)/len(ta|tb)]; nl=[links(a,b),float(links(a,b)>0)]; pwj=[len(pa&pb)/max(1,len(pa|pb))]; sv=list(Zd[a]+Zd[b])+list(np.abs(Zd[a]-Zd[b]))
    return ov,nl,pwj,sv
m=l.drug_a_name.isin(mapped)&l.drug_b_name.isin(mapped); L=l[m].reset_index(drop=True)
di={d:i for i,d in enumerate(drugs)}; ci={x:i for i,x in enumerate(cells)}
ID=np.zeros((len(L),len(drugs)+len(cells)))
for r,(a,b,cl) in enumerate(zip(L.drug_a_name,L.drug_b_name,L.cell_line)): ID[r,di[a]]+=1; ID[r,di[b]]+=1; ID[r,len(drugs)+ci[cl]]=1
G={'overlap':[],'links':[],'pathway':[],'svd':[]}; cache={}
for a,b in zip(L.drug_a_name,L.drug_b_name):
    k=tuple(sorted((a,b)))
    if k not in cache: cache[k]=gfe(*k)
    o,nl,pwj,sv=cache[k]; G['overlap'].append(o); G['links'].append(nl); G['pathway'].append(pwj); G['svd'].append(sv)
G={k:np.array(v,float) for k,v in G.items()}; GALL=np.hstack(list(G.values()))
pair=np.array(['|'.join(sorted((a,b))) for a,b in zip(L.drug_a_name,L.drug_b_name)]); score=L.synergy.values
R={'d1':dict(n_rows_all=len(l),n_drugs=len(drugs),n_mapped=len(mapped),unmapped=[d for d in drugs if d not in mapped],n_rows_both_mapped=len(L),n_pairs=len(set(pair)),n_cells=len(cells),pos_rate=round(float((score>=10).mean()),3),SN38_proxy='irinotecan')}
def mdl(): return HGC(max_iter=150,learning_rate=0.1,max_depth=4,random_state=0)
def folds(kind,seed):
    if kind=='random': return list(KFold(5,shuffle=True,random_state=seed).split(ID))
    grp={'cold_pair':pair,'cold_cell':L.cell_line.values}[kind]; u=sorted(set(grp)); rng=np.random.default_rng(seed); mp=dict(zip(u,rng.permutation(len(u)))); return list(GroupKFold(5).split(ID,groups=np.array([mp[g] for g in grp])))
def evalsplit(kind,feat,y,seeds=(0,1,2),F=None):
    res=[]; oof=np.zeros(len(y)); 
    for s in seeds:
        o=np.zeros(len(y))
        for tr,te in folds(kind,s):
            X=ID if F is None else np.hstack([ID,F]); o[te]=mdl().fit(X[tr],y[tr]).predict_proba(X[te])[:,1]
        res.append(roc_auc_score(y,o)); oof=o
    return float(np.mean(res)),oof
def gain(kind,y,F,seeds=(0,1,2)):
    a,_=evalsplit(kind,None,y,seeds,None); b,oo=evalsplit(kind,None,y,seeds,F); return dict(auroc_identity=round(a,4),auroc_with_gene=round(b,4),gain=round(b-a,4)),oo
y=(score>=10).astype(int)
def save(): json.dump(R,open('results/results.json','w'),indent=1,default=float)
lr=LogisticRegression(C=1,max_iter=500); b1=[]
for tr,te in folds('cold_pair',0): b1.append(roc_auc_score(y[te],lr.fit(ID[tr],y[tr]).decision_function(ID[te])))
R['d2_B1_LR_identity_cold_pair_auroc']=round(float(np.mean(b1)),4)
R['d2_random'],_=gain('random',y,GALL); print('d2',R['d2_random'],flush=True); save()
R['d3_cold_pair'],oo3=gain('cold_pair',y,GALL); R['d3_gate_gain_ge_0.01']=bool(R['d3_cold_pair']['gain']>=.01); print('d3',R['d3_cold_pair'],flush=True); save()
R['d4_cold_cell'],_=gain('cold_cell',y,GALL); print('d4',R['d4_cold_cell'],flush=True); save()
dg=[];
for d in mapped:
    te=np.where((L.drug_a_name==d)|(L.drug_b_name==d))[0]; tr=np.setdiff1d(np.arange(len(L)),te)
    if len(set(y[te]))<2: continue
    ai=mdl().fit(ID[tr],y[tr]).predict_proba(ID[te])[:,1]; X=np.hstack([ID,GALL]); ag=mdl().fit(X[tr],y[tr]).predict_proba(X[te])[:,1]; dg.append((d,roc_auc_score(y[te],ai),roc_auc_score(y[te],ag)))
R['d5_cold_drug']=dict(n_drugs=len(dg),mean_auroc_identity=round(float(np.mean([x[1] for x in dg])),4),mean_auroc_with_gene=round(float(np.mean([x[2] for x in dg])),4),mean_gain=round(float(np.mean([x[2]-x[1] for x in dg])),4),drugs_with_gain=int(sum(x[2]>x[1] for x in dg)),per_drug={x[0]:[round(x[1],3),round(x[2],3)] for x in dg}); print('d5',{k:v for k,v in R['d5_cold_drug'].items() if k!='per_drug'},flush=True); save()
R['d6_ablation_cold_pair']={k:gain('cold_pair',y,v,seeds=(0,1))[0]['gain'] for k,v in G.items()}; print('d6',R['d6_ablation_cold_pair'],flush=True); save()
up=sorted(set(pair)); rng=np.random.default_rng(7); perm=dict(zip(up,rng.permutation(up))); first={p:np.where(pair==p)[0][0] for p in up}; Fsh=np.array([GALL[first[perm[p]]] for p in pair])
R['d7_shuffled_gene_features_cold_pair'],_=gain('cold_pair',y,Fsh); R['d7_gate_gain_le_0.01']=bool(R['d7_shuffled_gene_features_cold_pair']['gain']<=.01); print('d7',R['d7_shuffled_gene_features_cold_pair'],flush=True); save()
cl=L.cell_line.values; R['d8_per_cell_auroc_cold_pair']={k:round(float(roc_auc_score(y[cl==k],oo3[cl==k])),3) for k in cells if (cl==k).sum()>=50 and len(set(y[cl==k]))>1}
v=list(R['d8_per_cell_auroc_cold_pair'].values()); R['d8_summary']=dict(min=min(v),max=max(v),mean=round(float(np.mean(v)),3))
R['d9_threshold_sensitivity']={}
for th in (5,10,20): R['d9_threshold_sensitivity'][str(th)]=gain('cold_pair',(score>=th).astype(int),GALL,seeds=(0,1))[0]
save()
rr=[]
for F in (None,GALL):
    X=ID if F is None else np.hstack([ID,F]); o=np.zeros(len(L))
    for tr,te in folds('cold_pair',0): o[te]=HGR(max_iter=150,max_depth=4,random_state=0).fit(X[tr],score[tr]).predict(X[te])
    rr.append(float(np.corrcoef(o,score)[0,1]))
R['d10_regression_pearson_cold_pair']=dict(identity=round(rr[0],4),with_gene=round(rr[1],4)); save(); print(json.dumps(R,indent=1,default=float)[:2500])
