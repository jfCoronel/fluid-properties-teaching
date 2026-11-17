"""Accuracy of manual fluid-property methods, emulated with CoolProp 6.4.1.

Study A: double linear interpolation in superheated-vapour tables built on
textbook-like grids (water and R-134a), compared with the equation of state.
Study B: a standard R-134a vapour-compression cycle solved with the equation of
state, by table interpolation, and with Monte Carlo chart-reading errors. The
steps of the table interpolation are stored under 'worked' (Table 2).

Writes analysis.json and the intermediate *.npy files to the current folder.
The random seed is fixed, so the results are reproducible.
"""
import numpy as np, json
from CoolProp.CoolProp import PropsSI
rng=np.random.default_rng(2026)
def P(o,a,av,b,bv,f): return PropsSI(o,a,av,b,bv,f)
# ---------- Study A: double interpolation in superheated tables ----------
def table_study(fluid, pgrid_MPa, tgrid_C, n=4000):
    pg=np.array(pgrid_MPa)*1e6; tg=np.array(tgrid_C)+273.15
    # tabulate (NaN where not superheated)
    tab={}
    for k in ['H','S','D']:
        M=np.full((len(pg),len(tg)),np.nan)
        for i,p in enumerate(pg):
            ts=P('T','P',p,'Q',1,fluid)
            for j,t in enumerate(tg):
                if t>ts+1e-6: M[i,j]=P(k,'P',p,'T',t,fluid)
        tab[k]=M
    tab['V']=1/tab['D']
    res=[]
    while len(res)<n:
        p=np.exp(rng.uniform(np.log(pg[0]),np.log(pg[-1]))); t=rng.uniform(tg[0],tg[-1])
        i=np.searchsorted(pg,p)-1; j=np.searchsorted(tg,t)-1
        if i<0 or j<0 or i>=len(pg)-1 or j>=len(tg)-1: continue
        cell={k:tab[k][i:i+2,j:j+2] for k in ['H','S','V']}
        if np.isnan(cell['H']).any(): continue   # needs saturation row: skip
        fp=(p-pg[i])/(pg[i+1]-pg[i]); ft=(t-tg[j])/(tg[j+1]-tg[j])
        r={'p':p,'t':t,'tsat':P('T','P',p,'Q',1,fluid)}
        for k in ['H','S','V']:
            c=cell[k]; est=(1-fp)*((1-ft)*c[0,0]+ft*c[0,1])+fp*((1-ft)*c[1,0]+ft*c[1,1])
            ex=P(k if k!='V' else 'D','P',p,'T',t,fluid); ex=1/ex if k=='V' else ex
            r[k]=(est-ex)/ex*100 if k=='V' else (est-ex)
        res.append(r)
    return res
steam_p=[0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0]
steam_t=list(range(50,601,50))+[700,800]
r134_p=[0.06,0.10,0.14,0.18,0.20,0.24,0.28,0.32,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60]
r134_t=list(range(-20,161,10))
out={}
for name,f,pg,tg in [('steam','Water',steam_p,steam_t),('R134a','R134a',r134_p,r134_t)]:
    res=table_study(f,pg,tg)
    sup=np.array([r['t']-r['tsat'] for r in res])
    dh=np.abs([r['H']/1e3 for r in res]); dv=np.abs([r['V'] for r in res]); ds=np.abs([r['S']/1e3 for r in res])
    near=sup<20
    out[name]={'n':len(res),
      'dh_mean':dh.mean(),'dh_p95':np.percentile(dh,95),'dh_max':dh.max(),
      'dv_mean':dv.mean(),'dv_p95':np.percentile(dv,95),'dv_max':dv.max(),
      'ds_mean':ds.mean(),'ds_max':ds.max(),
      # v varies as ~1/p: report the first (widest, 0.01-0.05 MPa) pressure interval apart
      'dv_firstp_mean':dv[np.array([r['p'] for r in res])<pg[1]*1e6].mean(),'dv_restp_mean':dv[np.array([r['p'] for r in res])>=pg[1]*1e6].mean(),
      'dh_near_mean':dh[near].mean(),'dh_far_mean':dh[~near].mean(),'dv_near_mean':dv[near].mean(),'dv_far_mean':dv[~near].mean(),'frac_near':near.mean()}
    np.save(f'{name}_res.npy',np.array([[r['p'],r['t'],r['tsat'],r['H'],r['S'],r['V']] for r in res]))
# ---------- Study B: refrigeration cycle COP ----------
f='R134a'; Te=263.15; Tc=313.15; eta=0.75
pe=P('P','T',Te,'Q',1,f); pc=P('P','T',Tc,'Q',0,f)
h1=P('H','T',Te,'Q',1,f); s1=P('S','T',Te,'Q',1,f)
h2s=P('H','P',pc,'S',s1,f); h2=h1+(h2s-h1)/eta; h3=P('H','T',Tc,'Q',0,f); h4=h3
cop=(h1-h4)/(h2-h1); T2=P('T','P',pc,'H',h2,f)-273.15
# table route: superheated table at 1.0 and 1.2 MPa, interpolate in s then in p
def h_at_ps(p,s,pg,tg,pasos=None):
    # pasos, if given, collects the rows read and each interpolation (Table 2 of the article)
    i=np.searchsorted(pg,p/1e6)-1; hs=[]
    for pp in (pg[i],pg[i+1]):
        T=np.array([t+273.15 for t in tg if t+273.15>P('T','P',pp*1e6,'Q',1,f)])
        S=np.array([P('S','P',pp*1e6,'Q',1,f)]+[P('S','P',pp*1e6,'T',t,f) for t in T]); H=np.array([P('H','P',pp*1e6,'Q',1,f)]+[P('H','P',pp*1e6,'T',t,f) for t in T])
        k=np.searchsorted(S,s)-1; hs.append(H[k]+(s-S[k])*(H[k+1]-H[k])/(S[k+1]-S[k]))
        if pasos is not None:
            Tfila=lambda j: 'sat' if j==0 else T[j-1]-273.15
            pasos.append({'p_MPa':pp,'rows':[[Tfila(j),S[j]/1e3,H[j]/1e3] for j in (k,k+1)],'h':hs[-1]/1e3})
    return hs[0]+(p/1e6-pg[i])*(hs[1]-hs[0])/(pg[i+1]-pg[i])
pasos=[]; h2s_tab=h_at_ps(pc,s1,r134_p,r134_t,pasos); h2_tab=h1+(h2s_tab-h1)/eta; cop_tab=(h1-h4)/(h2_tab-h1)
# chart reading Monte Carlo: each student reads h1,h3 and h2s (along isentrope) with sigma
N=20000
res_mc={}
for sig in [1.0,2.0,3.0]:
    e=rng.normal(0,sig*1e3,(N,3))
    H1=h1+e[:,0]; H3=h3+e[:,1]; H2s=h2s+e[:,2]
    C=(H1-H3)/((H2s-H1)/eta)
    res_mc[sig]={'mean':C.mean(),'sd':C.std(),'p5':np.percentile(C,5),'p95':np.percentile(C,95)}
    if sig==2.0: np.save('cop_mc.npy',C)
out['cycle']={'pe_kPa':pe/1e3,'pc_kPa':pc/1e3,'h1':h1/1e3,'s1':s1/1e3,'h2s':h2s/1e3,'h2':h2/1e3,'T2':T2,'h3':h3/1e3,'cop':cop,
  'qe':(h1-h4)/1e3,'w':(h2-h1)/1e3,'h2s_tab':h2s_tab/1e3,'cop_tab':cop_tab,'cop_tab_err_pct':(cop_tab-cop)/cop*100,'mc':res_mc,
  'T2_student42_s':P('S','P',pc,'T',315.15,f)/1e3,'h2_student':P('H','P',pc,'T',315.15,f)/1e3,
  'h4_student':P('H','T',Te,'Q',0.2,f)/1e3,'x4':P('Q','P',pe,'H',h4,f)}
out['cycle']['eta_student']=(h2s-h1)/(out['cycle']['h2_student']*1e3-h1)
out['worked']={'blocks':pasos,'p_MPa':pc/1e6,'s1':s1/1e3,'h2s_tab':h2s_tab/1e3,'h2s_eos':h2s/1e3,'T2s':P('T','P',pc,'S',s1,f)-273.15}
json.dump(out,open('analysis.json','w'),indent=1,default=float)
print(json.dumps(out,indent=1,default=float))
