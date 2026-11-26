// Experimental records for bulk PEtOx in water. Nulls retain unreported measurements.
(() => {
  const optical = {
    id: "hacioglu2024", authors: "Hacioglu ve ark.", year: 2024,
    title: "Thionation of Poly(2-ethyl-2-oxazoline) and Its Effect on Solubility and Cloud Point",
    doi: "10.1021/acs.macromol.4c01341", url: "https://doi.org/10.1021/acs.macromol.4c01341",
    supplement: "https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.4c01341/suppl_file/ma4c01341_si_001.pdf", dataUrl: "https://zenodo.org/records/13919956",
    method: "600 nm türbidimetri, ısıtma 1 °C/dk. Geçişin başlangıç ve bitiş sıcaklıklarının aritmetik ortalaması."
  };
  const thermal = {
    id: "gauche2019", authors: "Gauche ve Felisberti", year: 2019,
    title: "Colloidal Behavior of Cellulose Nanocrystals Grafted with Poly(2-alkyl-2-oxazoline)s",
    doi: "10.1021/acsomega.9b01269", url: "https://doi.org/10.1021/acsomega.9b01269",
    supplement: "https://pstorage-acs-6854636.s3.amazonaws.com/16191263/ao9b01269_si_001.pdf",
    method: "Yüksek duyarlılıklı diferansiyel taramalı kalorimetri (HSDSC), geçiş eğrisinin büküm noktası."
  };
  const commonOptical = {
    sourceId: optical.id, sample: "S0 (Aquazol 50)", polymer: "PEtOx", solvent: "H₂O",
    observable: "T_cloud", method: "turbidimetry", criterion: optical.method,
    wavelengthNm: 600, heatingRateCmin: 1, preparation: "4 °C’de en az bir gece bekletilmiş çözelti",
    concentrationBasis: "g polymer per L solution", volumeFraction: null, endGroups: null,
    mwGmol: 50000, molarMassBasis: "Tedarikçinin nominal Mw değeri", mnGmol: null,
    dispersity: null, dispersityRange: [3,4], table: "S2"
  };
  const records = [
    [0.1,88,"anchor"], [0.25,80,"interiorCheck"], [0.5,76,"anchor"],
    [1,72,"interiorCheck"], [2,68,"interiorCheck"], [5,65,"anchor"]
  ].map(([concentrationGL,tcpC,role]) => ({...commonOptical,concentrationGL,tcpC,role}));
  const repeat = {...commonOptical, concentrationGL:5, tcpC:66.0, sampleSdC:0.8, n:3,
    replicateC:[65.1,66.4,66.5], role:"within-study replicate series", table:"S3"};
  const thermalRecords = [
    ["PEtOx106",1,5000,10604,12831,1.20,87,"Yazarların doğrusal referansı"],
    ["PEtOx105",3,7000,10545,11916,1.13,85,"Dallanma miktarı belirtilmemiş"],
    ["PEtOx123",16,6100,12376,14480,1.17,84,"GPC'de yüksek mol kütleli omuz"],
    ["PEtOx133",24,5900,13379,14984,1.12,82,"GPC'de yüksek mol kütleli omuz"]
  ].map(([sample,reactionHours,mnNmrGmol,mnGmol,mwGmol,dispersity,transitionC,architecture]) =>
    ({sourceId:thermal.id,sample,reactionHours,mnNmrGmol,mnGmol,mwGmol,dispersity,transitionC,
      architecture,observable:"HSDSC inflection",criterion:thermal.method,concentrationGL:20,
      polymer:"PEtOx",solvent:"H₂O",concentrationBasis:"g polymer per L solution",volumeFraction:null,
      endGroups:null,molarMassBasis:"GPC, kloroform, polistiren kalibrasyonu",table:"S2"}));
  const external = {
    id:"zhang2014", authors:"Qilu Zhang", year:2014, title:"Smart Materials Based on Thermoresponsive Polymers",
    url:"https://biblio.ugent.be/publication/5878271", supplement:"https://backoffice.biblio.ugent.be/download/5878271/5878291",
    relatedArticle:"https://doi.org/10.1039/C7MH00016B",
    method:"Isıtmada %50 geçirgenlik, 5 K/dk, Tablo 2-2."
  };
  const externalRecords = [[1,71.5],[2.5,68.3],[5,66.6],[10,65.9]].map(([concentrationGL,tcpC]) =>
    ({sourceId:external.id,sample:"Aquazol 50",polymer:"PEtOx",solvent:"H₂O",observable:"T_cloud",
      concentrationGL,tcpC,heatingRateCmin:5,criterion:external.method,role:"external protocol",
      concentrationBasis:"g polymer per L solution",volumeFraction:null,mwGmol:50000,mnGmol:null,
      dispersity:null,dispersityRange:[3,4],molarMassBasis:"Nominal Mw",table:"2-2"}));
  window.PETOX_DATA = {sources:[optical,thermal,external],records,repeat,thermalRecords,externalRecords};
})();
