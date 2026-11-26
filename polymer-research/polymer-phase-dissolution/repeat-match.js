/* Match two labelled, single-bond attachment sites on a connected hydrocarbon repeat.
   OCL's interpretation of plain '*' is ambiguous; neutral monovalent sentinels
   preserve the two ports during canonicalisation. Only input without Cl is admitted. */
(function(root){
  'use strict';
  function graph(input,ocl){
    const raw=String(input??'').trim();
    if(!raw||raw.length>400)throw new RangeError('İki bağlantı noktası içeren, en fazla 400 karakterlik bir tekrar birimi gir.');
    if(/[@/\\]/.test(raw))throw new RangeError('Stereokimya belirtildi. Bu kohort stereokimyaya özgü bir eşleşme sağlamıyor.');
    if(/\s|\./.test(raw))throw new RangeError('Tek bir bağlı tekrar birimi gir.');
    let count=0;
    const maps=new Set();
    const marked=raw.replace(/\[\*(?::([1-9][0-9]*))?\]|\*/g,(_,map)=>{
      if(map&&maps.has(map))throw new RangeError('Bağlantı noktalarının etiketleri farklı olmalı.');
      if(map)maps.add(map);count++;return '~';
    });
    if(count!==2)throw new RangeError('İki uç için * veya [*:1] ve [*:2] kullan.');
    if(!/^[Cc0-9()%=#:\-~]+$/.test(marked)||raw.includes('~'))
      throw new RangeError('Bu eşleştirme nötr hidrokarbon tekrarları için C/c atom yazımını ve iki bağlantı noktasını destekler. Girilen yapı kapsam dışında.');
    let depth=0;
    for(const c of marked){
      if(c==='(')depth++;
      else if(c===')'&&--depth<0)throw new RangeError('Açılmamış dal parantezi var.');
    }
    if(depth!==0||/\(\)|^[=(#:\-)0-9%]|[=(#:\-]$|[-=#:]{2}/.test(marked))
      throw new RangeError('Dal veya bağ gösterimi tamamlanmamış.');
    let mol;
    try{mol=ocl.Molecule.fromSmiles(marked.replaceAll('~','Cl'));mol.ensureHelperArrays(ocl.Molecule.cHelperCIP);mol.validate();}
    catch{throw new RangeError('Tekrar SMILES grafiği okunamadı. Atom değerliklerini ve halka kapanışlarını kontrol et.');}
    if(mol.getFragments().length!==1)throw new RangeError('Tek bir bağlı tekrar birimi gerekli.');
    const ports=[];
    for(let i=0;i<mol.getAllAtoms();i++){
      const z=mol.getAtomicNo(i);
      if(mol.getAtomCharge(i)!==0||mol.getAtomRadical(i)!==0||mol.getOccupiedValence(i)>mol.getMaxValence(i))
        throw new RangeError('Atom yükü, radikali veya değerliği bu modelin kapsamı dışında.');
      if(z===17){
        if(mol.getConnAtoms(i)!==1||mol.getOccupiedValence(i)!==1||mol.getAtomicNo(mol.getConnAtom(i,0))!==6)
          throw new RangeError('Her bağlantı noktası bir karbon atomuna tek bağla bağlanmalı.');
        ports.push(i);
      }else if(z!==6)throw new RangeError('Bu kohort yalnızca karbon ve hidrojen içeren tekrarları kapsar.');
    }
    if(ports.length!==2)throw new RangeError('İki ayrı bağlantı noktası tanınmalı.');
    return {key:mol.getIDCode(),mol,ports};
  }
  function registry(catalog,ocl){
    const keys=new Map();
    for(const [id,r] of Object.entries(catalog.chemistries)){
      const entries=[r.repeat_smiles];
      // The same 1,4-diene can be cut at either side of an allylic single bond.
      if(id==='PBD')entries.push('*C=CCC*');
      if(id==='PI')entries.push('*C(C)=CCC*','*C=C(C)CC*');
      for(const smiles of entries){
        const {key}=graph(smiles,ocl);
        if(keys.has(key)&&keys.get(key)!==id)throw Error('Tekrar grafiği birden fazla kayıtla eşleşiyor.');
        keys.set(key,id);
      }
    }
    return keys;
  }
  function match(input,ocl,keys){
    const item=graph(input,ocl);
    const chemistry=keys.get(item.key);
    if(!chemistry)throw new RangeError('Bu bağlantı düzeni kayıtlı beş tekrar birimiyle eşleşmiyor. Karşılaştırma için kayıtlı bir kimya seçebilirsin.');
    return {chemistry, ...item};
  }
  function depiction(item){
    const m=item.mol.getCompactCopy();
    for(const i of item.ports){m.setAtomicNo(i,0);m.setAtomCustomLabel(i,'*');}
    return m.toSVG(340,140);
  }
  const api={registry,match,depiction};
  root.RepeatViscosity=api;if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof window==='undefined'?globalThis:window);
