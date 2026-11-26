"""Regression tests for osmometry units, thermodynamic consistency, and split leakage."""
import json
from pathlib import Path
import sys
from unittest import TestCase, main

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for directory in [HERE, ROOT / 'scripts']:
    if (directory / 'fit_transfer.py').exists():
        sys.path.insert(0, str(directory))
from fit_transfer import (G, VBAR_PS, VREF, predicted_a2, merge_tangent,
                          endpoints_summary, fit, read_cloud)
import numpy as np

def result():
    candidates = [HERE / 'transfer-results.js', ROOT / 'dist/transfer-results.js', ROOT / 'transfer-results.js']
    path = next(p for p in candidates if p.is_file())
    return json.loads(path.read_text().split('=', 1)[1].strip().removesuffix(';'))

class TransferChecks(TestCase):
    def test_a2_is_osmotic_virial(self):
        # Stable series for phi*f'−f with ideal chain term removed.
        for row in result()['models']:
            a,b,c = (row[k] for k in 'ABC')
            self.assertAlmostEqual(predicted_a2(row['theta_K'],a,b,c),0,places=12)
            p = 1e-6
            for t in (285.,310.):
                a0=a+b/t
                excess = p**2*(.5+p/3+p**2/4-a0+c-2*c*p)
                self.assertLess(abs(excess/(p/VBAR_PS)**2/VREF-predicted_a2(t,a,b,c)),3e-9)
            self.assertGreater(predicted_a2(310.,a,b,c),predicted_a2(300.,a,b,c))

    def test_conditional_scores_and_misses(self):
        truth=[dict(phi_lean=.01,phi_rich=.20),dict(phi_lean=.02,phi_rich=.12)]
        pred=[dict(phases=2,phi_lean=.02,phi_rich=.19),dict(phases=1,phi_lean=None,phi_rich=None)]
        score=endpoints_summary(truth,pred)
        self.assertEqual(score['n_missed'],1)
        self.assertAlmostEqual(score['mae_gap'],.06)
        self.assertAlmostEqual(score['rmse_pooled'],.01)
        self.assertIsNone(endpoints_summary(truth[1:],pred[1:])['mae_lean'])

    def test_common_tangent(self):
        pred=merge_tangent(100,300,.65,0,0)
        self.assertEqual(pred['phases'],2)
        self.assertAlmostEqual(pred['phi_lean'],.0097680307278,places=10)
        self.assertAlmostEqual(pred['phi_rich'],.2719683621054,places=10)
        self.assertEqual(merge_tangent(100,300,.48,0,0)['phases'],1)

    def test_holdout_not_fitted(self):
        evidence=result()
        rows=evidence['cloudPoints']
        self.assertEqual(sum(r['split']=='calibration' for r in rows),73)
        self.assertEqual(sum(r['split']=='holdout' for r in rows),24)
        clean_train=[r for r in rows if r['split']=='calibration']
        perturbed=[{**r,'T_K':7.} if r['split']=='holdout' else r for r in rows]
        for cloud in (clean_train,perturbed):
            params=fit(cloud,evidence['osmotic']['points'],1.)
            for k, value in zip('ABC',params):
                self.assertAlmostEqual(value,evidence['models'][1][k],places=6)
        self.assertEqual(len(evidence['osmotic']['points']),5)
        self.assertEqual(len(evidence['osmotic']['holdoutPoints']),3)
        self.assertTrue(not ({r['source'] for r in evidence['osmotic']['points']} &
                             {r['source'] for r in evidence['osmotic']['holdoutPoints']}))

    def test_endpoint_scoring(self):
        e=result()
        points=e['coexistence']['points']
        for m in e['models']:
            pred=[p['predictions'][m['id']] for p in points]
            for p, pr in zip(points,pred):
                self.assertTrue(0<=p['phi_lean']<p['phi_rich']<=1)
                self.assertTrue(np.isfinite(p['T_K']))
                self.assertTrue(pr['phases'] in [1,2])
            expected=endpoints_summary(points,pred)
            self.assertEqual(expected,m['endpointMetrics'])

if __name__=='__main__': main()
