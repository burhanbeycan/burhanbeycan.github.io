"""Plot three independent physical checks from the released model and recorded measurements."""
import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

root = Path(__file__).resolve().parent.parent
web = root / 'dist' if (root / 'dist').is_dir() else root
data = json.loads((web/'dense-results.js').read_text().split('=',1)[1].strip().rstrip(';'))
fold = next(x for x in data['folds'] if x['id'] == 'hold44')
mid = next(x for x in fold['rows'] if x['T_K']==317.15 and x['c_g_L']==178)
read = """global.window=globalThis;for(const f of ['phase-engine','osmotic-engine','dense-engine','dense-results','activity-engine'])require('./'+f+'.js');process.stdout.write(JSON.stringify(SolventActivity.compare()));"""
activity = json.loads(subprocess.check_output(['node','-e',read],cwd=web))
last = activity[-1]
colors = ['#233e38','#517382','#138065','#ac7854']
fig, axes = plt.subplots(1, 3, figsize=(12.8,5.35))
fig.subplots_adjust(left=.065,right=.986,top=.72,bottom=.23,wspace=.34)
fig.patch.set_facecolor('#f7f9f5')
fig.text(.058,.917,'Bir basınç düzeltmesi üç ayrı ölçütü karşılıyor mu?',size=22,fontfamily='DejaVu Serif',color='#163a32')
fig.text(.060,.845,'Polistiren / siklohekzan    •    Ölçüm, önceki model ve basınca uydurulan model',size=11,color='#476558')
fmt=lambda x: f'{x:.2f}'.replace('.',',')
labels=['Ölçüm','Önceki','D serbest']
for ax in axes:
    ax.set_facecolor('#f7f9f5')
    ax.spines[['left','top','right']].set_visible(False)
    ax.spines['bottom'].set_color('#9cab9f')
    ax.tick_params(axis='both',length=0,labelsize=10)
    ax.grid(axis='y',color='#dae4dc',zorder=0)
    ax.set_axisbelow(True)
    ax.set_xticks(range(3),labels)
    ax.set_xlim(-.55,2.55)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value,_:f'{value:g}'.replace('.',',')))
ax=axes[0]
pred=[mid['observed_kPa'],mid['base_kPa'],mid['extended_kPa']]
ax.bar(range(3),pred,color=colors[:3],width=.57)
ax.set_ylim(0,32);ax.set_title('Ozmotik basınç',loc='left',fontweight='bold');ax.set_ylabel('Π (kPa)')
for x,y in enumerate(pred):ax.text(x,y+.7,fmt(y),ha='center',size=10)
ax=axes[1]; a3=[6.5,data['base']['A3']*1e4,fold['A3']*1e4]
ax.axhline(0,color='#74847b',lw=1)
ax.bar(range(3),a3,color=colors[:3],width=.57)
ax.errorbar(0,6.5,yerr=2.5,fmt='none',color='#233e38',capsize=4)
ax.set_title('Üçüncü virial katsayısı',loc='left',fontweight='bold')
ax.set_ylabel('A₃ (10⁻⁴ mL² mol g⁻³)');ax.set_ylim(-9.5,11.1)
for x,y in enumerate(a3):ax.text(x,9.28 if x==0 else y+.43 if y>0 else y-.54,fmt(y),ha='center',va='bottom' if y>0 else 'top',size=10)
ax=axes[2];a=[last[k]['activity'] for k in ('observed','baseline','extended')]
ax.scatter(range(3),a,c=colors[:3],s=85,zorder=3)
ax.set_ylim(.850,.948);ax.set_title('Çözücü aktivitesi',loc='left',fontweight='bold');ax.set_ylabel('aₛ')
for x,y in enumerate(a):ax.text(x,y+.006,f'{y:.4f}'.replace('.',','),ha='center',size=10)
for ax,cap in zip(axes,['44 °C  ·  178 g/L  ·  Mₙ = 72 kg/mol','θ ≈ 34,5 °C  ·  ışık saçılması¹','44 °C  ·  φ = 0,660  ·  Mₙ = 72 kg/mol']):
    ax.text(0,-.23,cap,transform=ax.transAxes,fontsize=9.5,color='#365449')
fig.text(.058,.052,'Kaynaklar: Krigbaum ve Geymer (1959), Tablolar I ve III.  ¹ Nakamura ve ark. (1991), ayrı PS fraksiyonları; çizgi yazarların bildirdiği ±2,5 aralığıdır.',size=8.7,color='#405a4d')
fig.savefig(web/'PS_Basinc_Aktivite.png',dpi=160,facecolor=fig.get_facecolor())
