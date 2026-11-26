"""Structure, sample and solution descriptors. Masses use average atomic weights."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

HERE = Path(__file__).resolve().parent
CATALOG = json.loads((HERE / 'repeat_catalog.json').read_text())
SOLVENTS = dict(water='O', dmso='CS(=O)C', dmf='CN(C)C=O', dmac='CC(=O)N(C)C',
                nmp='CN1CCCC1=O', acetone='CC(=O)C')
UNIT_CONCENTRATION = {'g/mL': 1., 'g/dL': .01, 'g/L': .001, 'mg/mL': .001}
UNIT_INTRINSIC = {'mL/g': 1., 'dL/g': 100., 'L/g': 1000.}
OH = Chem.MolFromSmarts('[O;H1;+0;!$([O]-[C,S,P]=O)]')
ACETATE = Chem.MolFromSmarts('[O;X2]-[C;X3](=[O;X1])-[CH3]')
SODIUM_MASS = Chem.GetPeriodicTable().GetAtomicWeight(11)


def finite(x, field, *, minimum=None, maximum=None):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        raise ValueError(f'{field} must be a finite number')
    if minimum is not None and x < minimum or maximum is not None and x > maximum:
        raise ValueError(f'{field} outside [{minimum}, {maximum}]')
    return float(x)


def positive(x, field):
    if finite(x, field) <= 0:
        raise ValueError(f'{field} must be positive')
    return float(x)


def repeat_descriptors(psmiles: str) -> dict:
    mol = Chem.MolFromSmiles(psmiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('Use one connected repeat with exactly two attachment atoms')
    ends = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(ends) != 2 or mol.GetNumAtoms() <= 2:
        raise ValueError('Expected two terminal attachment atoms and at least one real atom')
    for atom in ends:
        if (atom.GetDegree() != 1 or atom.GetIsotope() or atom.GetFormalCharge() or
            atom.GetBonds()[0].GetBondType() != Chem.BondType.SINGLE or
            atom.GetNeighbors()[0].GetAtomicNum() == 0):
            raise ValueError('Use neutral, singly bonded terminal attachment atoms')
    counts = {
        'CH2': sum(a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 2 for a in mol.GetAtoms()),
        'alcoholic_or_phenolic_OH': len(mol.GetSubstructMatches(OH)),
        'acetate_ester': len(mol.GetSubstructMatches(ACETATE)),
        'HBD': rdMolDescriptors.CalcNumHBD(mol),
        'HBA': rdMolDescriptors.CalcNumHBA(mol),
        'formal_charge': Chem.GetFormalCharge(mol)
    }
    return {'canonical_graph': Chem.MolToSmiles(mol),
            'covalent_repeat_mass_g_mol': Descriptors.MolWt(mol),
            'sites_per_repeat': counts}


def statistical_repeat(components: list[dict], *, sodium_per_repeat: float = 0) -> dict:
    """Number-fraction compositional average per repeat. No fabricated sequence."""
    if not components:
        raise ValueError('At least one component is required')
    x = [finite(i['mole_fraction'], 'mole_fraction', minimum=0, maximum=1) for i in components]
    if not math.isclose(sum(x), 1, abs_tol=1e-8, rel_tol=0):
        raise ValueError('Repeat mole fractions must add to one')
    sodium = finite(sodium_per_repeat, 'sodium_per_repeat', minimum=0)
    repeats = [repeat_descriptors(i['psmiles']) for i in components]
    coval_mass = sum(f*m['covalent_repeat_mass_g_mol'] for f,m in zip(x,repeats))
    counts = {k: sum(f*m['sites_per_repeat'][k] for f,m in zip(x,repeats))
              for k in repeats[0]['sites_per_repeat']}
    if sodium > max(-counts['formal_charge'], 0)+1e-8:
        raise ValueError('Dry sodium counterion fraction exceeds negative charge equivalents')
    mass = coval_mass + sodium*SODIUM_MASS
    return {'components': [{'mole_fraction':f,'canonical_graph':m['canonical_graph']}
                           for f,m in zip(x,repeats)],
            'covalent_repeat_mass_g_mol': coval_mass, 'included_sodium_per_repeat':sodium,
            'repeat_mass_g_mol_on_selected_basis':mass, 'sites_per_repeat':counts,
            'sites_mmol_per_g_on_selected_basis':{k:1000*v/mass for k,v in counts.items()}}


def pva_composition(hydrolysis_mol_percent):
    f = finite(hydrolysis_mol_percent,'hydrolysis_mol_percent',minimum=0,maximum=100)/100
    return statistical_repeat([{'psmiles':CATALOG['PVA']['psmiles'],'mole_fraction':f},
                               {'psmiles':CATALOG['PVAc']['psmiles'],'mole_fraction':1-f}])


def molar_masses(mn=None, mw=None, dispersity=None, mv=None):
    """Leave missing moments explicit. Retain measured values subject to rounding tolerance."""
    a = positive(mn,'Mn') if mn is not None else None
    b = positive(mw,'Mw') if mw is not None else None
    v = positive(mv,'Mv') if mv is not None else None
    ratio = finite(dispersity,'dispersity',minimum=1) if dispersity is not None else None
    origin = {}
    if a is not None and b is not None:
        if b < a or ratio is not None and not math.isclose(b/a,ratio,rel_tol=.01):
            raise ValueError('Inconsistent Mn, Mw or dispersity (1% ratio rounding tolerance)')
        if ratio is None:
            ratio=b/a
            origin['dispersity']='Mw/Mn'
    elif ratio is not None and a is not None:
        b=ratio*a
        origin['Mw_g_mol']='dispersity × Mn'
    elif ratio is not None and b is not None:
        a=b/ratio
        origin['Mn_g_mol']='Mw/dispersity'
    return {'Mn_g_mol':a,'Mw_g_mol':b,'dispersity':ratio,'Mv_g_mol':v,'computed_fields':origin}


def mass_concentration(value, unit, density_g_mL=None):
    c=finite(value,'concentration',minimum=0)
    if unit in UNIT_CONCENTRATION:
        return c*UNIT_CONCENTRATION[unit]
    if unit not in ('wt%','mass_fraction'):
        raise ValueError(f'Unsupported concentration unit: {unit}')
    upper=100 if unit=='wt%' else 1
    if c>upper:
        raise ValueError('Mass fraction exceeds one')
    if density_g_mL is None:
        return None
    return c/upper*positive(density_g_mL,'density_g_mL')


def solvent_info(name):
    key=name.lower()
    if key not in SOLVENTS:
        raise ValueError('Solvent not registered. Use a documented neat solvent from SOLVENTS.')
    mol=Chem.MolFromSmiles(SOLVENTS[key])
    # Water is an acceptor and a donor site. The small-molecule Lipinski SMARTS
    # do not count its isolated oxygen and must not define its physical role.
    return {'id':key,'smiles':SOLVENTS[key],'molar_mass_g_mol':Descriptors.MolWt(mol),
            'donor_sites_per_molecule':1 if key=='water' else rdMolDescriptors.CalcNumHBD(mol),
            'acceptor_sites_per_molecule':1 if key=='water' else rdMolDescriptors.CalcNumHBA(mol)}


def describe(sample: dict[str,Any]) -> dict:
    """Accept either nominal PVA DH, any registered motif, or a known repeat mixture."""
    family=sample.get('family')
    form='exact repeat graph' if not sample.get('components') else 'composition average'
    if sample.get('components'):
        struct=statistical_repeat(sample['components'],sodium_per_repeat=sample.get('sodium_per_repeat',0))
    elif family=='PVA':
        if sample.get('hydrolysis_mol_percent') is None:
            struct=None
            form='parent chemistry known; quantitative vinyl-alcohol/acetate composition unknown'
        else:
            struct=pva_composition(sample['hydrolysis_mol_percent'])
            form='nominal vinyl-alcohol/vinyl-acetate composition'
    elif family in CATALOG:
        struct=statistical_repeat([dict(psmiles=CATALOG[family]['psmiles'],mole_fraction=1)],
                                  sodium_per_repeat=sample.get('sodium_per_repeat',0))
        form=CATALOG[family]['semantic']
    else:
        raise ValueError('Supply a registered family or explicit repeat components')
    moments=molar_masses(sample.get('Mn_g_mol'),sample.get('Mw_g_mol'),
                         sample.get('dispersity'),sample.get('Mv_g_mol'))
    if moments['Mn_g_mol'] is not None and not sample.get('Mn_origin') and 'Mn_g_mol' not in moments['computed_fields']:
        raise ValueError('Specify Mn_origin such as hypothetical target, measured SEC, or supplier')
    t=positive(sample['temperature_K'],'temperature_K')
    solvent=solvent_info(sample['solvent'])
    conc=sample.get('concentration')
    c=None
    if conc is not None:
        c=mass_concentration(conc['value'],conc['unit'],conc.get('solution_density_g_mL'))
    u=None
    iv=sample.get('intrinsic_viscosity')
    if iv is not None:
        if not iv.get('provenance') or iv.get('sample_id')!=sample.get('sample_id') or not sample.get('sample_id'):
            raise ValueError('Intrinsic viscosity requires provenance and a matched sample ID')
        if iv.get('solvent')!=solvent['id'] or not math.isclose(positive(iv['temperature_K'],'intrinsic T'),t,abs_tol=1e-6,rel_tol=0):
            raise ValueError('Intrinsic viscosity must match the sample solvent and temperature')
        if iv.get('unit') not in UNIT_INTRINSIC:
            raise ValueError('Unsupported intrinsic viscosity unit')
        intrinsic=positive(iv['value'],'intrinsic viscosity')*UNIT_INTRINSIC[iv['unit']]
        u=c*intrinsic if c is not None else None
    dpn = None
    ends = sample.get('combined_end_group_mass_g_mol')
    if ends is not None:
        ends = finite(ends, 'combined_end_group_mass_g_mol', minimum=0)
        basis = moments['Mn_g_mol']
        if not sample.get('end_group_basis'):
            raise ValueError('Specify the end-group basis that is additional to the repeat representation')
        if basis is not None and basis <= ends:
            raise ValueError('Nominal chain mass must exceed its end-group contribution')
        if basis is not None and struct is not None:
            dpn = (basis - ends) / struct['repeat_mass_g_mol_on_selected_basis']
    return {'sample_id':sample.get('sample_id'),'family':family,'composition_basis':form,
            'number_average_repeat_count':dpn, 'combined_end_group_mass_g_mol':ends,
            'end_group_basis':sample.get('end_group_basis'), 'structure':struct,'molar_masses':moments,'Mn_origin':sample.get('Mn_origin'),
            'solvent':solvent,'temperature_K':t,'mass_concentration_g_mL':c,'c_times_intrinsic':u,
            'intrinsic_viscosity_source':iv.get('provenance') if iv else None}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('sample_json',type=Path)
    args=p.parse_args()
    print(json.dumps(describe(json.loads(args.sample_json.read_text())),indent=2,ensure_ascii=False))
