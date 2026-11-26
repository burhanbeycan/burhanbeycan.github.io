/* Dilute-solution viscometry. Concentrations are g/L; intrinsic, reduced and
   inherent viscosities are mL/g. Percentage errors use each observed response. */
(function (root) {
  'use strict';
  function positive(x, name) {
    if (!Number.isFinite(x) || x<=0) throw new RangeError(name+' sonlu ve pozitif olmalı.');
    return x;
  }
  function concentration(c) {
    if (!Number.isFinite(c) || c<0) throw new RangeError('Derişim sonlu ve negatif olmayan bir değer olmalı.');
    return c/1000;
  }
  function line(xs,ys) {
    const mean=a=>a.reduce((s,x)=>s+x,0)/a.length;
    const mx=mean(xs), my=mean(ys);
    const xx=xs.reduce((s,x)=>s+(x-mx)**2,0);
    if (!(xx>0)) throw new RangeError('Doğrusal uyum için iki farklı apsis gerekir.');
    const slope=xs.reduce((s,x,i)=>s+(x-mx)*(ys[i]-my),0)/xx;
    return {slope,intercept:my-slope*mx,mx,my};
  }
  function distinctRows(rows) {
    if (!Array.isArray(rows) || rows.length<1) throw new RangeError('Veri dizisi boş.');
    const ids=new Set();
    for (const r of rows) {
      if(typeof r.id!=='string'||!r.id||ids.has(r.id)) throw new RangeError('Her ölçüm benzersiz kimlik taşımalı.');
      ids.add(r.id);
    }
  }
  function compatible(rows) {
    const vals = k=>new Set(rows.map(r=>r[k]).filter(v=>v!=null));
    if(vals('system').size>1||vals('solvent').size>1||vals('T_K').size>1)
      throw new RangeError('Mol kütlesi uyumu aynı polimer, çözücü ve sıcaklık gerektirir.');
  }
  function fitMH(rows) {
    distinctRows(rows);compatible(rows);
    const xs=rows.map(r=>Math.log(positive(r.Mw_g_mol,'M𝑤')));
    const ys=rows.map(r=>Math.log(positive(r.eta_mL_g,'[η]')));
    const {slope:a,intercept,mx,my}=line(xs,ys);
    return {a,K_mL_g:Math.exp(intercept),Mref_g_mol:Math.exp(mx),etaRef_mL_g:Math.exp(my),
      // K is the numeric coefficient for Mw expressed numerically in g/mol.
      minMw:Math.min(...rows.map(r=>r.Mw_g_mol)), maxMw:Math.max(...rows.map(r=>r.Mw_g_mol)), n:rows.length};
  }
  function predictMH(Mw,fit) {
    positive(Mw,'M𝑤');positive(fit.Mref_g_mol,'Referans mol kütlesi');positive(fit.etaRef_mL_g,'Referans intrinsik viskozitesi');
    if(!Number.isFinite(fit.a))throw new RangeError('Mark–Houwink üssü geçersiz.');
    return fit.etaRef_mL_g*Math.exp(fit.a*Math.log(Mw/fit.Mref_g_mol));
  }
  function errors(rows, observed, predicted) {
    if (!rows.length)return {n:0,mae:null,rmse:null,rmse_pct:null};
    let abs=0, sq=0, relSq=0;
    for(const r of rows) {
      const y=positive(r[observed],'Ölçülen yanıt'), p=r[predicted];
      if(!Number.isFinite(p))throw new RangeError('Sonlu model yanıtı gerekir.');
      const d=p-y;
      abs+=Math.abs(d);sq+=d*d;relSq+=(d/y)**2;
    }
    return {n:rows.length,mae:abs/rows.length,rmse:Math.sqrt(sq/rows.length),rmse_pct:100*Math.sqrt(relSq/rows.length)};
  }
  function mhSet(train,test) {
    distinctRows(train);distinctRows(test);compatible([...train,...test]);
    if(train.some(r=>test.some(v=>v.id===r.id)))throw new RangeError('Eğitim ve sınama fraksiyonları çakışıyor.');
    const model=fitMH(train);
    return test.map(r=>({
      id:r.id,study:r.study??null,system:r.system??null,Mw_g_mol:r.Mw_g_mol,
      measured_mL_g:positive(r.eta_mL_g,'[η]'),predicted_mL_g:predictMH(r.Mw_g_mol,model),a:model.a,
      domain:r.Mw_g_mol<model.minMw||r.Mw_g_mol>model.maxMw?'extrapolation':'interpolation',
      training_min_Mw:model.minMw,training_max_Mw:model.maxMw,train_ids:train.map(x=>x.id)
    }));
  }
  function mhReport(folds) {
    const e=errors(folds,'measured_mL_g','predicted_mL_g');
    return {folds,metrics:{...e,mae_mL_g:e.mae,rmse_mL_g:e.rmse}};
  }
  function studyBenchmark(train,out=[]) {
    if(train.length<3)throw new RangeError('Fraksiyon sınaması için en az üç eğitim numunesi gerekir.');
    const sources = new Set(train.map(row => row.study).filter(Boolean));
    if (out.some(row => sources.has(row.study))) throw new RangeError('Bağımsız çalışma sınamasında kaynaklar çakışıyor.');
    const loo=mhReport(train.flatMap(r=>mhSet(train.filter(x=>x.id!==r.id),[r])));
    return {loo,external:out.length?mhReport(mhSet(train,out)):null};
  }
  function scale(eta,c) {
    positive(eta,'[η]');
    return {u:eta*concentration(c),c_eta_g_L:1000/eta};
  }
  function dilute(eta,kH,c,eta0=null) {
    const scales=scale(eta,c);
    if(eta0!=null)positive(eta0,'Saf çözücünün dinamik viskozitesi');
    if(c>0&&!Number.isFinite(kH))throw new RangeError('Sonlu derişimde Huggins katsayısı gerekir.');
    const {u}=scales;
    const sp=u+(c===0?0:kH*u*u);
    if(!Number.isFinite(sp)||sp<=-1)throw new RangeError('Bu derişimde model fiziksel bir bağıl viskozite vermiyor.');
    const rel=1+sp, C=concentration(c);
    return {...scales,eta_rel:rel,eta_sp:sp,reduced_mL_g:c===0?eta:sp/C,
      inherent_mL_g:c===0?eta:Math.log1p(sp)/C,eta_mPa_s:eta0==null?null:eta0*rel};
  }
  function fitDilute(rows) {
    distinctRows(rows);compatible(rows);
    if(rows.length<2)throw new RangeError('Viskozite uyumu için en az iki farklı derişim gerekir.');
    const cs=rows.map(r=>positive(r.c_g_L,'Uyum derişimi')/1000);
    const rel=rows.map(r=>positive(r.eta_rel,'Bağıl viskozite'));
    const h=line(cs,rel.map((v,i)=>(v-1)/cs[i]));
    const k=line(cs,rel.map((v,i)=>Math.log(v)/cs[i]));
    positive(h.intercept,'Huggins sıfır-derişim kesişimi');
    positive(k.intercept,'Kraemer sıfır-derişim kesişimi');
    return {eta_mL_g:h.intercept,kH:h.slope/h.intercept**2,kK:k.slope/k.intercept**2,
      kraemerIntercept_mL_g:k.intercept,hugginsSlope_mL2_g2:h.slope,kraemerSlope_mL2_g2:k.slope,
      maxC_g_L:Math.max(...cs)*1000,minC_g_L:Math.min(...cs)*1000,n:rows.length};
  }
  function finiteSplit(rows,holdCount,middleLow=false) {
    let sorted=[...rows].sort((a,b)=>a.c_g_L-b.c_g_L);
    const levels=[...new Set(sorted.map(r=>r.c_g_L))];
    if(levels.length<holdCount+2)throw new RangeError('Derişim dış sınaması için yeterli bağımsız derişim düzeyi yok.');
    // The interior low-c check uses only the first three distinct concentration levels.
    // Its two-point fit provides an interpolation challenge with no residual degrees of freedom.
    const testLevels=new Set(middleLow ? [levels[1]] : levels.slice(-holdCount));
    if(middleLow)sorted=sorted.filter(r=>r.c_g_L<=levels[2]);
    const train=sorted.filter(r=>!testLevels.has(r.c_g_L));
    const fit=fitDilute(train);
    const folds=sorted.map(r=>{
      const c=concentration(r.c_g_L);
      const h=dilute(fit.eta_mL_g,fit.kH,r.c_g_L);
      const logK=fit.kraemerIntercept_mL_g*c+fit.kraemerSlope_mL2_g2*c*c;
      return {id:r.id,c_g_L:r.c_g_L,eta_rel:r.eta_rel,observed_eta_rel:r.eta_rel,
        huggins_eta_rel:h.eta_rel,kraemer_eta_rel:Math.exp(logK),
        observed_reduced_mL_g:(r.eta_rel-1)/c,observed_inherent_mL_g:Math.log(r.eta_rel)/c,
        u:h.u,role:testLevels.has(r.c_g_L)?'test':'train'};
    });
    const test=folds.filter(r=>r.role==='test');
    return {fit,train_ids:train.map(r=>r.id),test_ids:test.map(r=>r.id),folds,
      metrics:{huggins:errors(test,'observed_eta_rel','huggins_eta_rel'),kraemer:errors(test,'observed_eta_rel','kraemer_eta_rel')}};
  }
  function finiteBenchmark(rows) {
    distinctRows(rows);compatible(rows);
    const out=finiteSplit(rows,1);
    out.upperTwo=new Set(rows.map(r=>r.c_g_L)).size>=5?finiteSplit(rows,2):null;
    out.lowInterpolation=finiteSplit(rows,1,true);
    return out;
  }
  const api={fitMH,predictMH,studyBenchmark,mhSet,mhReport,scale,dilute,fitDilute,finiteBenchmark};
  root.Viscosity=api;
  if(typeof module!=='undefined' && module.exports)module.exports=api;
})(typeof window==='undefined'?globalThis:window);
