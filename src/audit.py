import json
SRC=open('src/run.py').read(); exec(SRC.split("def save()")[0].replace("R={'d1'","R0={'d1'"))
from sklearn.preprocessing import StandardScaler
y=(score>=10).astype(int); ups=sorted(set(pair)); pidx=np.array([ups.index(p) for p in pair])
Gs=StandardScaler().fit_transform(np.hstack([G['overlap'],G['links'],G['pathway']]))
def oof_avg(make,X,seeds=(0,1,2)):
    acc=np.zeros(len(y))
    for s in seeds:
        for tr,te in folds('cold_pair',s): acc[te]+=make().fit(X[tr],y[tr]).predict_proba(X[te])[:,1] if True else 0
    return acc/len(seeds)
def boot(pa,pb,B=200):
    rng=np.random.default_rng(0); bypair=[np.where(pidx==i)[0] for i in range(len(ups))]; out=[]; g0=roc_auc_score(y,pb)-roc_auc_score(y,pa)
    for _ in range(B):
        pick=rng.integers(0,len(ups),len(ups)); ii=np.concatenate([bypair[i] for i in pick])
        if len(set(y[ii]))<2: continue
        out.append(roc_auc_score(y[ii],pb[ii])-roc_auc_score(y[ii],pa[ii]))
    return dict(point=float(g0),lo=float(np.percentile(out,2.5)),hi=float(np.percentile(out,97.5)))
R={}
a=oof_avg(mdl,ID); b=oof_avg(mdl,np.hstack([ID,GALL])); R['Q1_HGB']=boot(a,b); R['Q1_HGB'].update(auroc_id=float(roc_auc_score(y,a)),auroc_gene=float(roc_auc_score(y,b))); print('Q1',R['Q1_HGB'],flush=True)
lrm=lambda:LogisticRegression(C=1,max_iter=500)
a2=oof_avg(lrm,ID); b2=oof_avg(lrm,np.hstack([ID,Gs])); R['Q2_LR']=boot(a2,b2); R['Q2_LR'].update(auroc_id=float(roc_auc_score(y,a2)),auroc_gene=float(roc_auc_score(y,b2))); print('Q2',R['Q2_LR'],flush=True)
g3=oof_avg(mdl,GALL); rng=np.random.default_rng(1); bypair=[np.where(pidx==i)[0] for i in range(len(ups))]; au=[]
for _ in range(200):
    pick=rng.integers(0,len(ups),len(ups)); ii=np.concatenate([bypair[i] for i in pick])
    if len(set(y[ii]))>1: au.append(roc_auc_score(y[ii],g3[ii]))
R['Q3_gene_only']=dict(auroc=float(roc_auc_score(y,g3)),lo=float(np.percentile(au,2.5)),hi=float(np.percentile(au,97.5))); print('Q3',R['Q3_gene_only'])
R['reading']=dict(Q1_below_0.01=bool(R['Q1_HGB']['hi']<0.01),Q2_below_0.01=bool(R['Q2_LR']['hi']<0.01),Q3_signal_alone=bool(R['Q3_gene_only']['lo']>0.5)) if False else None
R['reading']={'Q1_hi_below_0.01':bool(R['Q1_HGB']['hi']<0.01),'Q2_hi_below_0.01':bool(R['Q2_LR']['hi']<0.01),'Q3_lo_above_0.5':bool(R['Q3_gene_only']['lo']>0.5)}
json.dump(R,open('results/audit.json','w'),indent=1); print(R['reading'],'done')
