(() => {
  "use strict";

  const valences = {
    B: {0: [3], "-1": [4]}, C: {0: [4], "-1": [3], 1: [3]},
    N: {0: [3], "-1": [2], 1: [4]}, O: {0: [2], "-1": [1], 1: [3]},
    F: {0: [1], "-1": [0]}, P: {0: [3,5], 1: [4]}, S: {0: [2,4,6], "-1": [1], 1: [3]},
    Cl: {0: [1], "-1": [0]}, Br: {0: [1], "-1": [0]}, I: {0: [1], "-1": [0]},
    H: {0: [1]}, Si: {0:[4]}, Se: {0:[2,4,6]}, As: {0:[3,5]}, "*": {0: [1]}
  };
  const empty = () => ({atoms: [], bonds: []});
  const copy = obj => JSON.parse(JSON.stringify(obj));
  const atomName = a => a.port ? 'P' + a.port : a.element + (a.charge > 0 ? '⁺' : a.charge < 0 ? '⁻' : '');
  const attached = (g,id) => g.bonds.filter(b => b.a === id || b.b === id);
  const other = (b,id) => b.a === id ? b.b : b.a;
  const byId = (g,id) => g.atoms.find(a => a.id === id);
  const allowed = a => valences[a.element]?.[a.charge];
  const used = (g,a) => attached(g,a.id).reduce((s,b) => s+b.order, a.h || 0);
  const hCount = (g,a) => {
    if (a.h !== null) return a.h;
    if (a.port || a.element === 'H') return 0;
    const count = used(g,a);
    return Math.max(0, (allowed(a)?.find(n => n >= count) ?? count) - count);
  };
  const fail = msg => { throw new Error(msg); };

  function check(g) {
    const ports = new Set();
    for (const a of g.atoms) {
      const v = allowed(a);
      if (!v) fail(atomName(a) + ' için bu prototipte yük ve değerlik kuralı tanımlı değil.');
      const count = used(g,a);
      if (count > v[v.length-1]) fail(atomName(a) + ' (' + a.id + ') üzerinde bağ dereceleri ve sabit H toplamı ' + count + ' olur. İzin verilen üst sınır ' + v[v.length-1] + '.');
      if (a.port) {
        if (ports.has(a.port)) fail('P' + a.port + ' bağlantı noktası zaten mevcut.');
        ports.add(a.port);
        for (const b of attached(g,a.id)) {
          if (b.order !== 1 || byId(g,other(b,a.id))?.port) fail('Her bağlantı noktası tek bağla gerçek bir atoma bağlanmalı.');
        }
      }
    }
  }

  function readiness(g) {
    if (!g.atoms.length) return ['Tekrar birimini çizin veya SMILES girin.'];
    const problems = [];
    const hasPort=g.atoms.some(a=>a.port);
    for (const p of (hasPort ? [1,2] : [])) {
      const atom = g.atoms.find(a => a.port === p);
      if (!atom || attached(g,atom.id).length !== 1) problems.push('P' + p + ' için bir zincir bağlantısı tanımlayın.');
    }
    if (!g.atoms.some(a => !a.port)) problems.push('Birime gerçek atom ekleyin.');
    const seen = new Set([g.atoms[0].id]), queue = [g.atoms[0].id];
    for (let i=0; i<queue.length; i++) {
      for (const b of attached(g,queue[i])) {
        const next = other(b,queue[i]);
        if (!seen.has(next)) {seen.add(next); queue.push(next);}
      }
    }
    if (seen.size !== g.atoms.length) problems.push('Tekrar birimi tek bir bağlı yapı olmalı.');
    if (!g.readOnly && g.atoms.some(a => !a.port && a.h !== null && !allowed(a).includes(used(g,a)))) {
      problems.push('Köşeli parantezdeki H sayısını ve yükü denetleyin. Açık kabuklu atomlar kapsam dışında.');
    }
    return problems;
  }

  function atomCode(a) {
    if (a.port) return '[*:' + a.port + ']';
    if (a.h === null && a.charge === 0 && a.element !== 'H') return a.element;
    const hydrogen = a.h > 0 ? 'H' + (a.h > 1 ? a.h : '') : '';
    return '[' + a.element + hydrogen + (a.charge > 0 ? '+' : a.charge < 0 ? '-' : '') + ']';
  }
  const bondCode = order => ({1:'',2:'=',3:'#'})[order];
  function serialize(g) {
    if (!g.atoms.length) return '';
    const seen = new Set(), visitedBonds = new Set(), children = new Map(), rings = new Map(), roots = [];
    let ringId = 0;
    function explore(id) {
      seen.add(id); children.set(id,[]);
      for (const b of attached(g,id)) {
        if (visitedBonds.has(b)) continue;
        visitedBonds.add(b);
        const to = other(b,id);
        if (!seen.has(to)) { children.get(id).push([b,to]); explore(to); }
        else {
          ringId++;
          const tag = ringId < 10 ? String(ringId) : '%' + ringId;
          rings.set(id,(rings.get(id)||'') + bondCode(b.order) + tag);
          rings.set(to,(rings.get(to)||'') + tag);
        }
      }
    }
    for (const a of g.atoms) if (!seen.has(a.id)) {roots.push(a.id); explore(a.id);}
    function part(id) {
      const kids = children.get(id);
      const out = kids.map(([b,to],i) => {
        const s = bondCode(b.order) + part(to);
        return i+1 === kids.length ? s : '(' + s + ')';
      });
      return atomCode(byId(g,id)) + (rings.get(id)||'') + out.join('');
    }
    return roots.map(part).join('.');
  }

  function place(g,parentId) {
    const parent = byId(g,parentId);
    if (!parent) return {x:260, y:110};
    let best = -Infinity, dest;
    const list = attached(g,parentId).map(b => byId(g,other(b,parentId)));
    for (const deg of [-30,30,-150,150,90,-90,0,180,60,-60,120,-120]) {
      const rad = deg*Math.PI/180, x = parent.x + 53*Math.cos(rad), y=parent.y+53*Math.sin(rad);
      const sep = g.atoms.filter(a => a !== parent).reduce((n,a) => Math.min(n, Math.hypot(a.x-x,a.y-y)), 108);
      const spacing = list.reduce((s,a) => {
        const ang = Math.atan2(a.y-parent.y,a.x-parent.x);
        const delta = Math.acos(Math.min(1,Math.max(-1,Math.cos(rad-ang))));
        return Math.min(s,delta);
      },Math.PI);
      const score = Math.min(sep,85) - 19*Math.abs(spacing-2*Math.PI/3) - Math.abs(y-110)*0.07 + x*0.014;
      if (score > best) {best=score;dest={x,y};}
    }
    return dest;
  }

  // Relax a parsed graph for legibility in a 2D drawing. No chemical energy or 3D coordinates are inferred here.
  function layout(g) {
    if (g.atoms.length < 3 || g.atoms.length > 100 || g.bonds.length < g.atoms.length) return;
    for (let iter=0; iter<170; iter++) {
      const force = new Map(g.atoms.map(a => [a.id,{x:0,y:0}]));
      function shift(a,b,magnitude) {
        let dx=a.x-b.x, dy=a.y-b.y, d=Math.hypot(dx,dy);
        if (d < 0.01) { dx=0.6;dy=0.8;d=1; }
        const fx=magnitude*dx/d, fy=magnitude*dy/d;
        force.get(a.id).x+=fx;force.get(a.id).y+=fy;force.get(b.id).x-=fx;force.get(b.id).y-=fy;
      }
      for (let i=0;i<g.atoms.length;i++) {
        for (let j=i+1;j<g.atoms.length;j++) {
          const a=g.atoms[i], b=g.atoms[j], d=Math.max(4,Math.hypot(a.x-b.x,a.y-b.y));
          if (d<120) shift(a,b,Math.min(8,8000/(d*d)));
        }
      }
      for (const b of g.bonds) {
        const a=byId(g,b.a), z=byId(g,b.b);
        shift(a,z,(52-Math.hypot(a.x-z.x,a.y-z.y))*0.26);
      }
      for (const a of g.atoms) {
        const f = force.get(a.id);
        a.x+=Math.max(-4,Math.min(4,f.x));a.y+=Math.max(-4,Math.min(4,f.y));
      }
    }
  }

  function parseSmiles(input, syntaxOnly=false) {
    const s = String(input).trim();
    if (!s) fail('Önce tekrar biriminin SMILES gösterimini girin.');
    if (s.length > 1600) fail('Çizim önizlemesi için girdi sınırı 1.600 karakter.');
    if (/[\\/@]/.test(s)) fail('Stereokimya gösterimini denetleyin. @ ve @@ köşeli parantez içinde, / ve ters eğik çizgi yönlü bağlarda kullanılır.');
    if (s.includes('.')) fail('Tek bir bağlı tekrar birimi girin. Noktayla ayrılmış bileşenler ayrı yapılar olarak değerlendirilir.');
    const g=empty(), stack=[], rings=new Map();
    let current=null, order=null, prev='start', index=0;
    function connect(a,b,n) {
      if (a===b || g.bonds.some(k => (k.a===a && k.b===b) || (k.a===b && k.b===a))) fail('Aynı atom veya mevcut bağ üzerinden ikinci bir halka kapatılamaz.');
      g.bonds.push({a,b,order:n ?? 1});
    }
    while (index<s.length) {
      const c=s[index];
      if (c==='(') {
        if (current===null || order!==null || !['atom','ring','close'].includes(prev)) fail('Dal başlangıcını bir atomdan sonra yazın.');
        stack.push(current);prev='open';index++;continue;
      }
      if (c===')') {
        if (!stack.length || !['atom','ring','close'].includes(prev) || order!==null) fail('Dal parantezlerini ve dal içindeki atomları denetleyin.');
        current=stack.pop();prev='close';index++;continue;
      }
      if (['-','=','#'].includes(c)) {
        if (current===null || order!==null) fail('Bağ işaretini iki atom arasında kullanın.');
        order={'-':1,'=':2,'#':3}[c];prev='bond';index++;continue;
      }
      if (/\d/.test(c) || c==='%') {
        const match=s.slice(index).match(/^(?:%(\d{2})|(\d))/);
        if (current===null || !match || prev==='open') fail('Halka numarasını bir atomun hemen ardından girin.');
        const key=Number(match[1]||match[2]);
        if (rings.has(key)) {
          const start=rings.get(key);
          if (order!==null && start.order!==null && order!==start.order) fail('Halka kapama bağ dereceleri eşleşmiyor.');
          connect(start.id,current,order ?? start.order);
          rings.delete(key);
        } else rings.set(key,{id:current,order});
        index+=match[0].length;order=null;prev='ring';continue;
      }
      let atom, size;
      if (c==='[') {
        const stop=s.indexOf(']',index);
        if (stop<0) fail('Köşeli parantez kapanmadan girdi bitti.');
        size=stop-index+1;
        const v=s.slice(index+1,stop), port=v.match(/^\*:([12])$/);
        if (port) atom={element:'*',port:Number(port[1]),charge:0,h:0};
        else {
          const a=v.match(/^(Cl|Br|[BCNOFPSIH])(H(\d+)?)?([+-])?$/);
          if (!a) fail('Köşeli parantezli atom bu önizlemenin kapsamı dışında. P1 ve P2 için [*:1] ve [*:2] kullanın.');
          atom={element:a[1],port:null,charge:a[4]==='+'?1:a[4]==='-'?-1:0,h:a[2] ? Number(a[3]||1) : 0};
        }
      } else {
        const a=s.slice(index).match(/^(Cl|Br|[BCNOFPSI])/);
        if (!a) {
          if (/[bcnosp]/.test(c)) fail('Aromatik gösterim bu biçimde okunamadı. Halka kapanışlarını denetleyin veya açık tek ve çift bağlarla Kekulé gösterimini kullanın.');
          fail('Desteklenmeyen SMILES öğesi (' + (index+1) + '. karakter). Bağlantı noktaları [*:1] ve [*:2] biçimindedir.');
        }
        size=a[0].length;atom={element:a[0],port:null,charge:0,h:null};
      }
      if (g.atoms.length>=150) fail('Bu etkileşimli çizimin sınırı 150 atom ve bağlantı noktasıdır.');
      const id=g.atoms.length+1;
      g.atoms.push({...atom,id,...place(g,current)});
      if (current!==null) connect(current,id,order);
      current=id;order=null;index+=size;prev='atom';
    }
    if (stack.length || rings.size || order!==null || prev==='open') fail('Açık halka, parantez veya atom bekleyen bağ var.');
    if (!syntaxOnly) {check(g);layout(g);}
    return g;
  }

  // Preserve the imported string. Coordinates and atom numbering may be normalized by OpenChemLib.
  function advancedSmiles(raw) {
    if (!window.PolyOcl || !raw || raw.length>1600 || /\s/.test(raw) || raw.includes('.')) return null;
    const lex=/\[([^\]]+)\]|Cl|Br|[BCNOFPSIbcnops]|%\d\d|[0-9()=#:\/\\-]/gy;
    let i=0, normalized='', aromatics=0, advanced=false;
    const ports=[];
    while (i<raw.length) {
      lex.lastIndex=i;
      const t=lex.exec(raw);
      if (!t) return null;
      const token=t[0], bracket=t[1];
      if (bracket!==undefined) {
        if (/^\*:[12]$/.test(bracket)) {ports.push(Number(bracket[2])); normalized+=token;}
        else {
          // No atom maps on real atoms, no arbitrary query/SMARTS syntax.
          const b=bracket.match(/^(\d+)?(Cl|Br|Si|Se|se|As|as|[BCNOFPSIHbcnops])(@@?)?(H\d*)?([+-](?:[1-3])?)?$/);
          if (!b) return null;
          if (/^[a-z]/.test(b[2])) aromatics++;
          normalized+='C';
          advanced ||= !!(b[1] || b[3] || /^[a-z]/.test(b[2]) || ['Si','Se','As'].includes(b[2]) || (b[5] && b[5].length>1));
        }
      } else if (/^[bcnops]$/.test(token)) {aromatics++;normalized+='C';advanced=true;}
      else if (token==='/' || token==='\\' || token===':') {normalized+='-';advanced=true;}
      else normalized+=token;
      i=lex.lastIndex;
    }
    if (!advanced) return null;
    parseSmiles(normalized,true);
    if (new Set(ports).size!==ports.length) fail('P1 ve P2 işaretlerinin her biri en çok bir kez kullanılabilir.');
    try {
      const ocl=window.PolyOcl(), m=ocl.Molecule.fromSmiles(raw);
      // OCL treats '*' as a query fragment, suppressing implicit H on every atom.
      // Restore normal implicit H handling and keep the two source-mapped ports outside atom chemistry.
      m.setFragment(false);
      m.ensureHelperArrays(ocl.Molecule.cHelperCIP);
      const totalCharge=Array.from({length:m.getAllAtoms()},(_,k)=>m.getAtomCharge(k)).reduce((sum,v)=>sum+v,0);
      try {m.validate();} catch(err) {
        if (!(totalCharge && /(?:^|: )unbalanced atom charge$/.test(String(err?.message)))) throw err;
      }
      if (m.getAllAtoms()>150) fail('Bu etkileşimli çizimin sınırı 150 atom ve bağlantı noktasıdır.');
      if (!m.getAllAtoms() || m.getFragments().length!==1) fail('Tek bir bağlı yapı girin.');
      const g=empty();
      let aromaticCount=0;
      for (let k=0;k<m.getAllAtoms();k++) {
        const port=m.getAtomMapNo(k)||null;
        if (port && !ports.includes(port)) fail('Tanımlanamayan atom eşlemesi.');
        if (m.getAtomRadical(k)) fail('Açık kabuklu atomlar bu modelin kapsamı dışındadır.');
        if (!port) {
          const spec={element:m.getAtomLabel(k),charge:m.getAtomCharge(k)};
          const v=allowed(spec), actual=m.getOccupiedValence(k)+m.getImplicitHydrogens(k);
          if (!v) fail(atomName(spec)+' için bu modelde yük ve değerlik kapsamı tanımlı değil.');
          if (!v.includes(actual)) fail(atomName(spec)+' ('+(k+1)+') için bağ dereceleri ve H toplamı '+actual+'. Desteklenen değerlik '+v.join(' / ')+'.');
        }
        if (m.isAromaticAtom(k)) aromaticCount++;
        g.atoms.push({id:k+1,element:port?'*':m.getAtomLabel(k),port,charge:m.getAtomCharge(k),h:port?0:m.getImplicitHydrogens(k),x:260+44*m.getAtomX(k),y:110+44*m.getAtomY(k)});
      }
      for (let k=0;k<m.getAllBonds();k++) g.bonds.push({a:1+m.getBondAtom(0,k),b:1+m.getBondAtom(1,k),order:m.getBondOrder(k)});
      if (aromatics>aromaticCount) fail('Aromatik olarak yazılan atomlar için uygun aromatik halka doğrulanamadı.');
      for (const a of g.atoms.filter(a=>a.port)) {
        const b=attached(g,a.id);
        if (b.length!==1 || b[0].order!==1 || byId(g,other(b[0],a.id)).port) fail('Her zincir bağlantısı gerçek bir atoma tek bağla bağlanmalı.');
      }
      for (const a of g.atoms.filter(a=>a.port)) m.setAtomCustomLabel(a.id-1,'P'+a.port);
      g.depiction=m.toSVG(495,200,'ps-imported-structure',{autoCrop:true,autoCropMargin:12});
      g.readOnly=true;
      return g;
    } catch (e) {
      if (!e?.message || /Class\$|atom valence exceeded|Exception|Error:/i.test(e.message)) fail('OpenChemLib girdiyi doğrulayamadı. Atom değerliklerini, aromatikliği ve stereokimyayı denetleyin.');
      throw e;
    }
  }

  function mount(root) {
    if (!root || typeof root.querySelector!=='function') throw new Error('Çizim kapsayıcısı gerekli.');
    if (root._polySketch) return root._polySketch;
    let graph=empty(), focusId=null, focusBond=null, pendingLink=false, element='C', mode='draw', history=[];
    let box={x:0,y:0,w:520,h:220};
    root.classList.add('ps-editor');
    root.innerHTML=`
      <div class="ps-switch" aria-label="Yapı giriş yöntemi">
        <button type="button" data-mode="draw" aria-pressed="true">Atom atom çiz</button>
        <button type="button" data-mode="smiles" aria-pressed="false">SMILES gir</button>
      </div>
      <div class="ps-input" hidden>
        <label>Tekrar birimi veya sonlu yapı için SMILES
          <textarea rows="2" data-smiles spellcheck="false" autocapitalize="off" autocomplete="off" placeholder="Örn. [*:1]CCO[*:2]"></textarea>
        </label>
        <div class="ps-input-foot"><span>Doğrusal birim uçları <code>[*:1]</code> ve <code>[*:2]</code></span>
        <button class="ps-main" type="button" data-action="load">Yapıyı aç</button></div>
      </div>
      <div class="ps-worktop">
        <div class="ps-palette" aria-label="Atom veya zincir bağlantısı" data-palette></div>
        <div class="ps-bondtype"><label>Bağ <select data-order aria-label="Bağ derecesi"><option value="1">Tek</option><option value="2">Çift</option><option value="3">Üçlü</option></select></label></div>
      </div>
      <div class="ps-board">
        <svg class="ps-svg" viewBox="0 0 520 220" aria-label="Etkileşimli iki boyutlu kimyasal yapı çizimi"></svg>
        <div class="ps-rich" hidden></div>
        <div class="ps-blank"><div class="ps-empty-circle">+</div><b>İlk atomu yerleştir</b><span>Bir element seçip tuvale dokun.</span></div>
        <div class="ps-board-badge">2D</div>
      </div>
      <div class="ps-actions">
        <button type="button" data-action="add">+ Atom ekle</button>
        <button type="button" data-action="link" aria-pressed="false">İki atomu bağla</button>
        <button type="button" data-action="remove">Sil</button>
        <button type="button" data-action="undo">Geri al</button>
        <button type="button" data-action="clear">Yeni çizim</button>
        <button type="button" class="ps-example" data-action="sample">Örnek aç</button>
      </div>
      <div class="ps-readout" role="status" aria-live="polite"></div>
      <details class="ps-inspect"><summary>Seçili atomun yükü ve hidrojenleri</summary>
        <div class="ps-props">
          <label>Formal yük <select data-charge><option value="-1">−1</option><option value="0" selected>0</option><option value="1">+1</option></select></label>
          <label>H sayısı <select data-h><option value="auto">Örtük</option><option>0</option><option>1</option><option>2</option><option>3</option><option>4</option></select></label>
          <button type="button" data-action="apply">Atoma uygula</button>
        </div>
      </details>
      <div class="ps-smiles-line"><span>SMILES</span><output data-code aria-label="Çizimden üretilen SMILES">—</output></div>
      <div class="ps-meta" data-meta></div>
      <details class="ps-scope"><summary>Girdi kapsamı</summary><p>B, C, N, O, F, P, S, Cl, Br ve I ile kurulan tek bileşenli yapılar, belirli formal yükler, dallar, halka kapanışları ve tek, çift, üçlü bağlar okunur. H atomu [H] olarak girilebilir. OpenChemLib, desteklenen aromatik, stereokimyasal, izotoplu ve Si içeren SMILES yapılarını özgün girdiyle korur. Bu yapılarda değişikliği SMILES alanından yapabilirsiniz. Radikaller ve koordinasyon kompleksleri bu prototipin kapsamı dışındadır. Yerel kontrol bağ dereceleri ile tanımlı H sayısını denetler. Sentezlenebilirlik, kararlılık ve kuvvet alanı kapsamı ayrı değerlendirme gerektirir.</p></details>`;
    const $=selector=>root.querySelector(selector);
    const $$=selector=>Array.from(root.querySelectorAll(selector));
    const svg=$('.ps-svg'), log=$('.ps-readout'), input=$('[data-smiles]');
    const selectPalette=$('[data-palette]');
    const elements=['C','N','O','S','F','Cl','P1','P2'];
    const extras=['B','P','Br','I'];
    const names={C:'Karbon',N:'Azot',O:'Oksijen',S:'Kükürt',F:'Flor',Cl:'Klor',P1:'Baş bağlantı noktası',P2:'Son bağlantı noktası',P:'Fosfor',B:'Bor',Br:'Brom',I:'İyot'};
    selectPalette.innerHTML=elements.map(e=>'<button type="button" data-element="'+e+'" aria-label="'+names[e]+' ('+e+')" title="'+names[e]+'" aria-pressed="'+(e==='C')+'">'+e+'</button>').join('')
      + '<label class="ps-other"><span class="ps-sr">Diğer elementler</span><select data-other aria-label="Diğer elementler"><option value="">Diğer</option>'+extras.map(e=>'<option value="'+e+'">'+e+'</option>').join('')+'</select></label>';
    function hint(s,error=false) {log.textContent=s;log.classList.toggle('ps-error',error);}
    function payload() {
      const g=copy(graph), errors=readiness(g);
      const {depiction,rawSmiles,...data}=g;
      const hasPort=g.atoms.some(a=>a.port);
      return {...data,smiles:rawSmiles || serialize(g),canonical:false,ready:hasPort && !errors.length,finiteReady:!hasPort && !errors.length,errors};
    }
    function emit() {window.dispatchEvent(new CustomEvent('polymer-graph-change',{detail:payload()}));}
    function commit(next,msg,selected=null,imported=false) {
      if (!imported) delete next.rawSmiles;
      try {if (!next.readOnly) check(next);} catch(e) {hint(e.message,true);return false;}
      history.push(copy(graph));
      if (history.length>60) history.shift();
      graph=next;focusId=selected;focusBond=null;pendingLink=false;render();hint(msg);emit();return true;
    }
    function load(smiles) {
      const raw=String(smiles).trim();
      let next;
      try {
        try {next=parseSmiles(raw);} catch(simpleError) {
          next=advancedSmiles(raw);
          if (!next) throw simpleError;
        }
        next.rawSmiles=raw;
      } catch(e) {
        hint(e.message,true);
        window.dispatchEvent(new CustomEvent('polymer-smiles-review',{detail:{smiles:raw,reason:e.message,acceptedInEditor:false}}));
        return false;
      }
      const msg=next.readOnly?'OpenChemLib yapıyı okudu. Özgün gösterimi korumak için bu yapıda düzenlemeyi SMILES alanından yapın.':'Girdi okundu. Çizimde bir atoma dokunarak düzenleyebilirsiniz.';
      const ok=commit(next,msg,null,true);
      if (ok) input.value=raw;
      return ok;
    }
    function activeElement(e) {
      element=e;
      $$('[data-element]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.element===e)));
      $('[data-other]').value=extras.includes(e)?e:'';
    }
    function modeSet(m) {
      mode=m;$$('[data-mode]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.mode===m)));
      $('.ps-input').hidden=m!=='smiles';
      if (m==='smiles') {input.value=graph.rawSmiles || serialize(graph) || input.value; input.focus();}
    }
    function bondSelected(g) {
      return focusBond && g.bonds.find(b=>b.a===focusBond.a && b.b===focusBond.b);
    }
    function add(pos) {
      if (graph.readOnly) return hint('Bu yapıyı SMILES alanından düzenleyebilirsiniz.');
      if (graph.atoms.length>=150) return hint('Çizimin sınırı 150 atom ve bağlantı noktasıdır.',true);
      if (!focusId && graph.atoms.length) return hint('Yeni atomun bağlanacağı atomu önce seçin.',true);
      const n=copy(graph), port=['P1','P2'].includes(element)?Number(element[1]):null;
      const id=Math.max(0,...n.atoms.map(a=>a.id))+1;
      n.atoms.push({id,element:port?'*':element,port,charge:0,h:port?0:null,...(pos||place(n,focusId))});
      if (focusId!==null) n.bonds.push({a:focusId,b:id,order:Number($('[data-order]').value)});
      commit(n,port ? 'P'+port+' eklendi. Bu işaret zincirin bağlanacağı yeri tanımlar.' : names[element]+' eklendi.', id);
    }
    function chooseAtom(id) {
      if (pendingLink && focusId!==null) {
        if (id===focusId) return hint('Bağ için ikinci, farklı bir atom seçin.',true);
        const n=copy(graph);
        if (attached(n,id).some(b=>other(b,id)===focusId)) return hint('Bu iki atom zaten bağlı. Bağ derecesi için aradaki çizgiyi seçin.',true);
        n.bonds.push({a:focusId,b:id,order:Number($('[data-order]').value)});
        commit(n,'İki atom arasında bağ kuruldu.',id);return;
      }
      focusId=id;focusBond=null;
      const a=byId(graph,id);
      render();hint(atomName(a)+' ('+id+') seçildi. '+(pendingLink?'Şimdi ikinci atomu seçin.':'Yeni atom seçili atoma bağlanır.'));
    }
    function chooseBond(a,b) {
      focusId=null;pendingLink=false;focusBond={a,b};
      $('[data-order]').value=String(bondSelected(graph).order);
      render();hint('Bağ seçildi. Bağ menüsünden derecesini değiştirebilirsiniz.');
    }
    function actions(key) {
      if (key==='add') return add(null);
      if (key==='sample') {
        const s='[*:1]CCO[*:2]';
        if (load(s)) {if (mode==='smiles') input.value=s;hint('Poli(etilen oksit) tekrar birimi açıldı. Çizimi değiştirebilir veya geri alabilirsiniz.');}
      }
      if (key==='load') load(input.value);
      if (key==='clear' && graph.atoms.length) {commit(empty(),'Yeni çizim açıldı. Önceki yapı için Geri al seçeneğini kullanabilirsiniz.');input.value='';}
      if (key==='undo' && history.length) {
        graph=history.pop();focusId=null;focusBond=null;pendingLink=false;render();hint('Son işlem geri alındı.');emit();
      }
      if (key==='link') {
        pendingLink=!pendingLink;focusBond=null;render();
        hint(pendingLink?(focusId===null?'Sırayla bağlanacak iki atomu seçin.':'Bağlanacak ikinci atomu seçin.'):'Bağ ekleme seçimi kapatıldı.');
      }
      if (key==='remove') {
        const n=copy(graph);
        if (focusId!==null) {n.atoms=n.atoms.filter(a=>a.id!==focusId);n.bonds=n.bonds.filter(b=>b.a!==focusId && b.b!==focusId);}
        else if (focusBond) n.bonds=n.bonds.filter(b=>b.a!==focusBond.a || b.b!==focusBond.b);
        else return;
        commit(n,'Seçim silindi. Geri al ile önceki çizime dönebilirsiniz.');
      }
      if (key==='apply') {
        if (focusId===null) return;
        const n=copy(graph), a=byId(n,focusId);
        if (a.port) return;
        const ch=Number($('[data-charge]').value), hydrogen=$('[data-h]').value;
        if (ch!==0 && hydrogen==='auto') return hint('Yüklü atom için H sayısını açıkça seçin. Örneğin bağlı H içermeyen bir atomda 0 kullanılır.',true);
        a.charge=ch;a.h=hydrogen==='auto'?null:Number(hydrogen);
        commit(n,'Atom ayarları uygulandı.',focusId);
      }
    }
    function svgPaint() {
      const rich=$('.ps-rich');
      rich.hidden=!graph.readOnly;
      if (graph.readOnly) {
        svg.replaceChildren();rich.innerHTML=graph.depiction;$('.ps-blank').hidden=true;return;
      }
      rich.replaceChildren();
      const ns='http://www.w3.org/2000/svg';
      const node=(tag,attrs,label)=>{
        const el=document.createElementNS(ns,tag);
        for (const [k,v] of Object.entries(attrs)) el.setAttribute(k,v);
        if (label!==undefined) el.textContent=label;
        return el;
      };
      const xs=graph.atoms.map(a=>a.x), ys=graph.atoms.map(a=>a.y);
      const left=Math.min(0,...xs.map(x=>x-30)),right=Math.max(520,...xs.map(x=>x+30));
      const top=Math.min(0,...ys.map(y=>y-32)),bottom=Math.max(220,...ys.map(y=>y+32));
      box={x:left,y:top,w:right-left,h:bottom-top};
      svg.setAttribute('viewBox',[box.x,box.y,box.w,box.h].join(' '));svg.replaceChildren();
      for (const b of graph.bonds) {
        const a=byId(graph,b.a), z=byId(graph,b.b),dx=z.x-a.x,dy=z.y-a.y,len=Math.max(.1,Math.hypot(dx,dy));
        const selected=focusBond?.a===b.a && focusBond?.b===b.b;
        const group=node('g',{'data-bonda':b.a,'data-bondb':b.b,tabindex:'0',role:'button','aria-label':'Bağ '+b.a+'–'+b.b+', derece '+b.order,class:'ps-bond'+(selected?' is-selected':'')});
        group.append(node('line',{x1:a.x,y1:a.y,x2:z.x,y2:z.y,class:'ps-hit'}));
        const offsets=b.order===1?[0]:b.order===2?[-2.2,2.2]:[-4,0,4];
        offsets.forEach(off=>group.append(node('line',{x1:a.x-off*dy/len,y1:a.y+off*dx/len,x2:z.x-off*dy/len,y2:z.y+off*dx/len,class:'ps-stick'})));
        svg.append(group);
      }
      for (const a of graph.atoms) {
        const sel=a.id===focusId, isPort=!!a.port;
        const group=node('g',{transform:'translate('+a.x+' '+a.y+')','data-atom':a.id,tabindex:'0',role:'button','aria-label':atomName(a)+', atom '+a.id,class:'ps-atom ps-atom-'+(a.port?'port':a.element)+(sel?' is-selected':'')});
        group.append(node('circle',{r:17,class:'ps-atom-target'}));
        group.append(node('circle',{r:isPort?13:12,class:'ps-disc'}));
        group.append(node('text',{'text-anchor':'middle',dy:'.34em',class:'ps-symbol'},atomName(a)));
        if (!isPort && a.h>0) group.append(node('text',{x:10,y:-8,class:'ps-small-h'},'H'+(a.h>1?a.h:'')));
        svg.append(group);
      }
      $('.ps-blank').hidden=!!graph.atoms.length;
    }
    function render() {
      svgPaint();
      const chosen=byId(graph,focusId), g=payload();
      $('[data-action="undo"]').disabled=!history.length;
      $('[data-action="clear"]').disabled=!graph.atoms.length;
      $('[data-action="add"]').disabled=!!graph.readOnly;
      $('[data-order]').disabled=!!graph.readOnly;
      $$('[data-element],[data-other]').forEach(b=>b.disabled=!!graph.readOnly);
      $('[data-action="remove"]').disabled=!chosen && !focusBond;
      $('[data-action="link"]').disabled=graph.atoms.length<2 || !!graph.readOnly;
      $('[data-action="link"]').setAttribute('aria-pressed',String(pendingLink));
      $('[data-action="apply"]').disabled=!chosen || !!chosen.port;
      $('[data-charge]').disabled=!chosen || !!chosen.port;
      $('[data-h]').disabled=!chosen || !!chosen.port;
      $('[data-charge]').value=String(chosen?.charge||0);
      $('[data-h]').value=chosen?.h!=null ? String(chosen.h) : 'auto';
      $('[data-code]').textContent=g.smiles||'—';
      const p=graph.atoms.filter(a=>a.port).length, n=graph.atoms.length-p;
      $('[data-meta]').textContent=graph.atoms.length?n+' atom · '+p+'/2 bağlantı noktası · '+(g.ready?'Doğrusal birim tanımlı':g.finiteReady?'Sonlu yapı tanımlı. Tekrar birimi oluşturmak için P1 ve P2 ekleyin.':g.errors.join(' ')):'';
    }
    function pointOnCanvas(ev) {
      const matrix=svg.getScreenCTM?.();
      if (matrix && typeof svg.createSVGPoint==='function') {
        const pt=svg.createSVGPoint();pt.x=ev.clientX;pt.y=ev.clientY;
        const p=pt.matrixTransform(matrix.inverse());
        return {x:Math.max(box.x+15,Math.min(p.x,box.x+box.w-15)),y:Math.max(box.y+15,Math.min(p.y,box.y+box.h-15))};
      }
      return place(graph,focusId);
    }
    root.addEventListener('click',e=>{
      const target=e.target.closest('button,[data-atom],[data-bonda],svg');
      if (!target || !root.contains(target)) return;
      if (target.hasAttribute('data-element')) activeElement(target.dataset.element);
      else if (target.hasAttribute('data-mode')) modeSet(target.dataset.mode);
      else if (target.hasAttribute('data-action')) actions(target.dataset.action);
      else if (target.hasAttribute('data-atom')) chooseAtom(Number(target.dataset.atom));
      else if (target.hasAttribute('data-bonda')) chooseBond(Number(target.dataset.bonda),Number(target.dataset.bondb));
      else if (target===svg) {
        if (pendingLink) return hint('Çizimde bağlanacak atomu seçin.',true);
        add(pointOnCanvas(e));
      }
    });
    svg.addEventListener('keydown',e=>{
      if (!['Enter',' '].includes(e.key)) return;
      const t=e.target.closest('[data-atom],[data-bonda]');if (!t) return;
      e.preventDefault();
      if (t.hasAttribute('data-atom')) {const id=Number(t.dataset.atom);chooseAtom(id);svg.querySelector('[data-atom="'+id+'"]')?.focus();}
      else {const {bonda,bondb}=t.dataset;chooseBond(Number(bonda),Number(bondb));svg.querySelector('[data-bonda="'+bonda+'"][data-bondb="'+bondb+'"]')?.focus();}
    });
    $('[data-order]').addEventListener('change',()=>{
      if (!focusBond) return;
      const n=copy(graph), b=bondSelected(n), val=Number($('[data-order]').value);
      const prior=b.order;b.order=val;
      if (!commit(n,'Bağ derecesi güncellendi.')) $('[data-order]').value=String(prior);
    });
    $('[data-other]').addEventListener('change',e=>{if (e.target.value) activeElement(e.target.value);});
    input.addEventListener('keydown',e=>{if (e.key==='Enter' && (e.ctrlKey || e.metaKey)) {e.preventDefault();load(input.value);}});
    render();hint('Yeni bir tekrar birimi çizebilir veya SMILES girebilirsiniz.');
    const api={getGraph:payload,loadSmiles:load};
    root._polySketch=api;emit();return api;
  }
  window.PolySketch={mount,parseSmiles,advancedSmiles,serialize,readiness,hydrogenCount:hCount};
})();
