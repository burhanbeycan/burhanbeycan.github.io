// Kato, Katsuki & Takahashi, Macromolecules 17 (1984), Table I, p. 1728.
// Each row pairs three experimental measurements on the identical pullulan fraction.
(function (root) {
  'use strict';
  const source = {
    id:'kato1984', title:'Static and Dynamic Solution Properties of Pullulan in a Dilute Solution',
    authors:'T. Kato, T. Katsuki, A. Takahashi', year:1984, doi:'10.1021/ma00139a016',
    url:'https://doi.org/10.1021/ma00139a016', table:'I', page:1728,
    solvent:'water-NaN3-0.2wt', T_K:298.15, solvent_label:'Suda kütlece %0,2 sodyum azit',
    Rg_method:'Statik ışık saçılması, z-ortalama kare yarıçapın karekökü.',
    Mw_method:'Statik ışık saçılması, (Mw)LS sütunu.',
    viscosity_method:'25 °C Ubbelohde kılcal viskozimetri.',
    dispersity_method:'PEO standartlı evrensel kalibrasyon ile sulu GPC.'
  };
  // Published numeric columns: (Mw)LS / 10⁴, Rg (nm), [η] / 10² in cm³/g, Mw/Mn.
  // For P400 the first three significant digits of [η] are legible, a trailing mark is ambiguous.
  const records=[
    ['P2000A',236,75.2,3.38,1.27,'primary'],
    ['F81',210,66.0,3.10,1.28,'primary'],
    ['P2000B',149,56.0,2.50,1.21,'primary'],
    ['F83',105,45.0,2.25,1.25,'primary'],
    ['P800',81.1,37.5,1.75,1.19,'primary'],
    ['F85',57.2,31.5,1.45,1.20,'primary'],
    ['P400',38.4,25.1,1.02,1.11,'eta-truncated'],
    ['P200',18.8,17.0,.700,1.11,'primary'],
    ['P100',10.7,15.0,.442,1.09,'radius-parenthesized']
  ];
  const flags={
    primary:'',
    'eta-truncated':'[η] sütununda 1,02 sonrasındaki işaret/basamak belirsizdir. 102 mL/g yaklaşık değeri yalnız duyarlılık hesabı içindir.',
    'radius-parenthesized':'Kaynak Rg değerini parantez içinde vermiştir. Duyarlılık satırıdır.'
  };
  const rows=records.map(([id,mass,radius,intrinsic,dispersity,quality])=>({
    id,sample:id,row_key:'kato1984/pullulan/'+id,family:'pullulan',study:source.id,
    solvent:source.solvent,T_K:source.T_K,Mw_g_mol:Math.round(mass*1e4),Rg_nm:radius,
    eta_mL_g:Number((intrinsic*100).toFixed(3)),dispersity,quality,use_primary:quality==='primary',
    notes:flags[quality],source_doi:source.doi
  }));
  const api={source,rows};
  if (root) root.FAMILY_EXTRA=api;
  if (typeof module!=='undefined'&&module.exports) module.exports=api;
})(typeof window==='undefined'?null:window);
