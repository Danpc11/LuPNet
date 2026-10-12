#!/usr/bin/env python3
"""Every number quoted in the text of the manuscript that is not in a figure or table script, recomputed from the stored
cohorts; also supplementary table S10 (PH by FVC range). Writes results/figures/manuscript_numbers.tsv.

    python3 src/manuscript_numbers.py
"""
import sys, numpy as np, pandas as pd
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_io import cases
from ph_stats import point, bootstrap, prevalence_ci
from scipy.stats import pearsonr, spearmanr
from sklearn.mixture import GaussianMixture
from sklearn.metrics import roc_auc_score
ROOT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
F=os.path.join(ROOT,'results','ph_cohorts'); OUT=os.path.join(ROOT,'results','figures'); os.makedirs(OUT,exist_ok=True)
def last(tag):
    R=[]
    for i,d in cases(F,tag):
        r0,x=d.iloc[0],d.iloc[-1]
        R.append(dict(mPAP=x.mPAP,PVR=x.PVR,FVC=100*x.FVC/r0.FVC,DLCO=100*x.DLCO/r0.DLCO,A=float(x.get('vasc_A',0)),eta=x.get('eta_p',np.nan),hpv=x.get('hpv_max',np.nan)))
    return pd.DataFrame(R)
C,A,AF,H,HT=last('big2'),last('adv2'),last('altF'),last('hom'),last('het')
M=C[(C.FVC>=50)&(C.FVC<=90)].reset_index(drop=True)
o=[]
def add(k,v): o.append((k,str(v))); print(f'{k:55s} {v}')
add('FVC range big2/adv2/DLCO big2',f'{C.FVC.min():.0f}-{C.FVC.max():.0f} / {A.FVC.min():.0f}-{A.FVC.max():.0f} / {C.DLCO.min():.0f}-{C.DLCO.max():.0f}')
add('vasc prevalence big2/adv2/altF',f'{100*(C.A>0).mean():.0f}/{100*(A.A>0).mean():.0f}/{100*(AF.A>0).mean():.0f}')
nv=M[M.A==0]; v=M[M.A>0]
add('mild n / no-vasc n / vasc n',f'{len(M)}/{len(nv)}/{len(v)}')
add('no vasc mPAP',f'{nv.mPAP.mean():.1f}±{nv.mPAP.std():.1f} max {nv.mPAP.max():.1f}')
add('vasc mPAP, %PH',f'{v.mPAP.mean():.1f}±{v.mPAP.std():.1f}, {100*(v.mPAP>=25).mean():.0f}%')
add('PH without vasc (all 200)',int(((C.mPAP>=25)&(C.A==0)).sum()))
ph=M.mPAP>=21; add('DLCO ratio at 21',f'{M.DLCO[ph].mean()/M.DLCO[~ph].mean():.2f}')
add('precap 2022 no vasc (mild)',f'{int(((nv.mPAP>20)&(nv.PVR>2)).sum())}/{len(nv)} maxPVR {nv.PVR.max():.2f} max mPAP {nv.mPAP.max():.1f}')
for lab,f in (('mPAP>20 & PVR>2',lambda D:(D.mPAP>20)&(D.PVR>2)),('mPAP>20',lambda D:D.mPAP>20),('PVR>5',lambda D:D.PVR>5)):
    a,b,c=prevalence_ci(M,f); add(f'mild {lab}',f'{a:.0f}% ({b:.0f}-{c:.0f})')
add('PVR>5 vasc/novasc (mild)',f'{int(((M.PVR>5)&(M.A>0)).sum())}/{int(((M.PVR>5)&(M.A==0)).sum())}')
P=pd.concat([C,A]).reset_index(drop=True)
for lab,D in (('big2 FVC<50',C[C.FVC<50]),('adv2 FVC<50',A[A.FVC<50]),('pooled FVC<50',P[P.FVC<50]),('pooled FVC 34-66',P[(P.FVC>=34)&(P.FVC<=66)])):
    D=D.reset_index(drop=True); a,b,c=prevalence_ci(D,lambda D:D.mPAP>=25); s,sl,sh=prevalence_ci(D,lambda D:D.mPAP>40); r5,_,_=prevalence_ci(D,lambda D:D.PVR>5)
    add(lab,f'n={len(D)} PH {a:.0f}% ({b:.0f}-{c:.0f}); >40 {s:.0f}% ({sl:.0f}-{sh:.0f}); PVR>5 {r5:.0f}%; FVC PH/noPH {D.FVC[D.mPAP>=25].mean():.1f}/{D.FVC[D.mPAP<25].mean():.1f}; nPH {int((D.mPAP>=25).sum())} nNo {int((D.mPAP<25).sum())}')
nov=C[C.A==0]; dl=np.polyfit(nov.FVC,np.log(nov.DLCO),2); res=(np.log(C.DLCO)-np.polyval(dl,C.FVC)).values.reshape(-1,1)
g2=GaussianMixture(2,random_state=0).fit(res); g1=GaussianMixture(1,random_state=0).fit(res); L=g2.predict(res)==int(np.argmin(g2.means_.ravel())); PH=(C.mPAP>=25).values
add('phenotype',f'dBIC {g1.bic(res)-g2.bic(res):.1f}; n={L.sum()}; vasc {100*(C.A[L]>0).mean():.0f}%; PH {100*PH[L].mean():.0f}%; rest vasc {100*(C.A[~L]>0).mean():.0f}% PH {100*PH[~L].mean():.0f}%; sens {100*(L&PH).sum()/PH.sum():.0f}% ({(L&PH).sum()}/{PH.sum()}) spec {100*((~L)&(~PH)).sum()/(~PH).sum():.0f}%; vasc captured {100*(L&(C.A.values>0)).sum()/(C.A>0).sum():.0f}%')
Hm=H[(H.FVC>=50)&(H.FVC<=90)]; HTm=HT[(HT.FVC>=50)&(HT.FVC<=90)]
add('hom',f'n={len(Hm)} {Hm.mPAP.mean():.1f}±{Hm.mPAP.std():.1f} max {Hm.mPAP.max():.1f} PH {int((Hm.mPAP>=25).sum())} rFVC {pearsonr(Hm.mPAP,Hm.FVC)[0]:.2f} rDLCO {pearsonr(Hm.mPAP,Hm.DLCO)[0]:.2f}')
add('het',f'n={len(HTm)} sd {HTm.mPAP.std():.2f} max {HTm.mPAP.max():.1f} rho eta {spearmanr(HTm.mPAP,HTm.eta)[0]:.2f} hpv {spearmanr(HTm.mPAP,HTm.hpv)[0]:.2f}')
AFm=AF[(AF.FVC>=50)&(AF.FVC<=90)].reset_index(drop=True); pa=point(AFm); ca=bootstrap(AFm)
add('altF',f'n={len(AFm)} vasc {int((AFm.A>0).sum())}; '+'; '.join(f'{k} {pa[k]:.2f} ({ca[k][0]:.2f}-{ca[k][1]:.2f})' for k in pa))
# Furukawa threshold (mPAP >= 21 mmHg) DLCO ratio with bootstrap interval
g=np.random.default_rng(1); b=[]
for _ in range(1000):
    S=M.iloc[g.integers(0,len(M),len(M))]; ph=S.mPAP>=21; b.append(S.DLCO[ph].mean()/S.DLCO[~ph].mean())
ph=M.mPAP>=21; add('DLCO ratio at 21 mmHg (95% CI)', f'{M.DLCO[ph].mean()/M.DLCO[~ph].mean():.2f} ({np.percentile(b,2.5):.2f}-{np.percentile(b,97.5):.2f})')
# table S10: PH by FVC range, pooled calibrated and advanced cohorts
P=pd.concat([C,A]).reset_index(drop=True); rows=[]
for lab,lo,hi in (('< 35%',0,35),('35-40%',35,40),('40-45%',40,45),('45-50%',45,50),('50-60%',50,60),('60-70%',60,70),('70-80%',70,81),('< 50% (main text)',0,50),('34-66% (range of transplant candidates)',34,66.0001)):
    D=P[(P.FVC>=lo)&(P.FVC<hi)].reset_index(drop=True); v,l,h=prevalence_ci(D,lambda D:D.mPAP>=25); s_,sl,sh=prevalence_ci(D,lambda D:D.mPAP>40)
    rows.append([lab,len(D),f'{v:.0f} ({l:.0f}-{h:.0f})',f'{s_:.0f} ({sl:.0f}-{sh:.0f})'])
pd.DataFrame(rows,columns=['FVC','n','PH, % (95% CI)','mPAP > 40 mmHg, % (95% CI)']).to_csv(os.path.join(OUT,'TableS10_PH_by_FVC.tsv'),sep='\t',index=False)
# partial association of PaO2 with lesion severity at given FVC and DLCO (trajectory points in established disease)
R=[]
for tag in ('big2','adv2'):
    for i,d in cases(F,tag):
        A_=float(d.vasc_A.iloc[0]); t0=float(d.t_vasc0.iloc[0]); r0=d.iloc[0]
        for _,x in d.iterrows():
            if x.converged and x.healthy<=0.85:
                R.append(dict(FVC=100*x.FVC/r0.FVC,DLCO=100*x.DLCO/r0.DLCO,PaO2=x.PaO2,Aeff=A_*(1-np.exp(-(x.years+t0)/2.0)) if A_>0 else 0.0))
V=pd.DataFrame(R); Xc=np.c_[np.ones(len(V)),V.FVC,V.DLCO,V.FVC**2,V.DLCO**2,V.FVC*V.DLCO]
res=lambda y: y-Xc@np.linalg.lstsq(Xc,y,rcond=None)[0]
add('partial r(PaO2, effective lesion | FVC, DLCO)', f'{pearsonr(res(V.PaO2.values),res(V.Aeff.values))[0]:.2f}')
pd.DataFrame(o,columns=['quantity','value']).to_csv(os.path.join(OUT,'manuscript_numbers.tsv'),sep='\t',index=False)
print('written results/figures/manuscript_numbers.tsv and TableS10_PH_by_FVC.tsv')
