"""Figures 1-3 of the article. Run analysis.py first, in the same folder.

Usage: python figs.py <output folder>
"""
import numpy as np, json, os, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D=sys.argv[1]; os.makedirs(D, exist_ok=True); A=json.load(open("analysis.json"))
BLUE='#2a78d6'; ORANGE='#eb6834'; INK='#2b2b2a'; MUTED='#77766f'
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.edgecolor':MUTED,'axes.labelcolor':INK,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'pdf.fonttype':42})
# ---- Fig 1 timeline ----
ev=[(1614,'Napier: tables of logarithms'),(1904,'Mollier h–s chart'),(1911,'Carrier psychrometric\nformulae'),
(1936,'Keenan & Keyes steam tables'),(1972,'Hand-held scientific calculator'),(1979,'Electronic spreadsheet'),
(1993,'EES: equations +\nproperty functions'),(1995,'IAPWS-95 equation\nof state for water'),(2014,'CoolProp\n(open source)'),
(2019,'WebAssembly\nweb standard')]
fig,ax=plt.subplots(figsize=(7.2,3.2))
eras=[(1600,1970,'Printed tables & charts'),(1970,1990,'Calculators &\nspreadsheets'),(1990,2010,'Equation\nsolvers'),(2010,2030,'Open &\nweb-based')]
xs=lambda y: np.interp(y,[1600,1900,1950,2030],[0,1.2,3.2,10])
shades=['#eef3fb','#dbe7f8','#c4d8f4','#a9c7ef']
for (a,b,l),c in zip(eras,shades):
    ax.axvspan(xs(a),xs(b),ymin=0.42,ymax=0.58,color=c,lw=0); ax.text((xs(a)+xs(b))/2,0.5,l,ha='center',va='center',fontsize=8,color=INK)
for k,(y,l) in enumerate(ev):
    up,lvl=[(1,0),(0,0),(1,1),(0,1),(1,0),(0,0),(1,1),(0,1),(1,0),(0,0)][k]; yy=0.62 if up else 0.38; yt=(0.78+0.2*lvl) if up else (0.22-0.2*lvl)
    ax.plot([xs(y)]*2,[yy,yt+(-0.04 if up else 0.04)],color=MUTED,lw=0.7)
    ax.plot(xs(y),yy,'o',ms=4,color=BLUE)
    ax.text(xs(y),yt,f'{y}\n{l}',ha='center',va='bottom' if up else 'top',fontsize=6.8,color=INK)
ax.set_ylim(-0.2,1.2); ax.set_xlim(-0.6,10.2); ax.axis('off')
fig.tight_layout(); [fig.savefig(f'{D}/Figure_1_timeline.{e}') for e in ('png','tif','pdf')]; plt.close()
# ---- Fig 2 interpolation error vs superheat ----
fig,axs=plt.subplots(1,2,figsize=(7.2,2.9),sharey=False)
for ax,(name,lab) in zip(axs,[('steam','Water (steam table grid)'),('R134a','R-134a (refrigerant table grid)')]):
    r=np.load(f'{name}_res.npy'); sup=r[:,1]-r[:,2]; dh=np.abs(r[:,3])/1e3
    ax.scatter(sup,np.maximum(dh,1e-4),s=3,color=BLUE,alpha=0.35,lw=0)
    ax.set_yscale('log'); ax.set_xscale('log'); ax.set_xlabel('Degree of superheat (K)'); ax.set_ylabel('|h interpolated − h EOS| (kJ/kg)')
    ax.set_title(lab,fontsize=9,color=INK,loc='left')
    ax.axhline(1,color=ORANGE,lw=1,ls='--'); ax.text(sup.min()*1.1,1.25,'1 kJ/kg',ha='left',fontsize=7,color=INK)
    ax.grid(True,which='major',color='#e6e5e0',lw=0.5)
fig.tight_layout(); [fig.savefig(f'{D}/Figure_2_interpolation_error.{e}') for e in ('png','tif','pdf')]; plt.close()
# ---- Fig 3 COP spread ----
C=np.load('cop_mc.npy'); c=A['cycle']
fig,ax=plt.subplots(figsize=(4.6,2.8))
ax.hist(C,bins=60,range=(2.0,4.2),color=BLUE,edgecolor='white',lw=0.4)
ax.axvline(c['cop'],color=INK,lw=1.2); ax.text(3.62,ax.get_ylim()[1]*0.97,f"EOS / software: {c['cop']:.3f}\nTable interpolation: {c['cop_tab']:.3f}",fontsize=7.5,color=INK,va='top')
p5,p95=np.percentile(C,[5,95])
for v in (p5,p95): ax.axvline(v,color=ORANGE,lw=1,ls='--')
ax.text(p95+0.04,ax.get_ylim()[1]*0.62,f'90% of readings:\n{p5:.2f} – {p95:.2f}',fontsize=7.5,color=INK)
ax.set_xlabel('Refrigeration COP obtained'); ax.set_ylabel('Simulated readings'); 
fig.tight_layout(); [fig.savefig(f'{D}/Figure_3_chart_reading_COP.{e}') for e in ('png','tif','pdf')]; plt.close()
print('ok')
