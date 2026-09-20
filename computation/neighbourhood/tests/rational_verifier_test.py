import importlib.util
from pathlib import Path
from fractions import Fraction as Q
import unittest
s=importlib.util.spec_from_file_location('verifier',Path(__file__).parents[1]/'tools/verify_rational.py')
v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
c=importlib.util.spec_from_file_location('checker',Path(__file__).parents[1]/'tools/check_certificate.py')
checker=importlib.util.module_from_spec(c);c.loader.exec_module(checker)

class ExactVerifierTests(unittest.TestCase):
 def test_outward_square_root(self):
  for x in [Q(1,3),Q(1,4),Q(999999,1000000),Q(2,7)]:
   y=v.sqrt_upper(x);self.assertGreaterEqual(y*y,x);self.assertLess((y-Q(1,2**80))**2,x)
 def test_real_separation_all_orientations(self):
  for e in [[0,0],[0,1],[1,0],[1,1]]:
   r=v.verify({'multipliers':[['1/4','0'],['1/4','0']],'orientation':e});self.assertTrue(r['disconnection_certified']);self.assertGreater(Q(r['minimum_squared_margin']),0)
 def test_no_connectedness_from_survival(self):
  for e in [[0,0],[0,1],[1,1]]:
   r=v.verify({'multipliers':[['3/5','0'],['3/5','0']],'orientation':e},6,20000);self.assertEqual(r['status'],'unresolved');self.assertFalse(r['connectedness_certified'])
 def test_quarter_turn_separation(self):
  for e in [[0,0],[0,1],[1,1]]:
   r=v.verify({'multipliers':[['0','3/5'],['0','3/5']],'orientation':e});self.assertTrue(r['disconnection_certified']);self.assertGreater(r['depth_reached'],1)
 def test_connected_quarter_turn_not_rejected(self):
  r=v.verify({'multipliers':[['0','3/4'],['0','3/4']],'orientation':[1,1]},5,10000);self.assertFalse(r['disconnection_certified'])
 def test_budget_is_unresolved(self):
  r=v.verify({'multipliers':[['3/5','0'],['3/5','0']],'orientation':[0,0]},14,1);self.assertFalse(r['disconnection_certified']);self.assertIn('budget',r['reason'])
 def test_float_inputs_rejected(self):
  with self.assertRaises(ValueError):v.verify({'multipliers':[[0.25,'0'],['1/4','0']],'orientation':[0,0]})
 def test_noncontractive_rejected(self):
  with self.assertRaises(ValueError):v.verify({'multipliers':[['1','0'],['1/4','0']],'orientation':[0,0]})
 def test_semilinear_composition(self):
  a=((Q(2),Q(3)),(Q(1,3),Q(2,5)),1,Q(4,5));b=((Q(-1),Q(4)),(Q(2,7),Q(-1,4)),0,Q(3,4));z=(Q(2,9),Q(-3,11))
  def apply(p,z):
   m=v.multiply(p[1],v.conjugate(z,p[2]));return(p[0][0]+m[0],p[0][1]+m[1])
  self.assertEqual(apply(v.compose(a,b),z),apply(a,apply(b,z)))
 def test_independent_exact_core_equality_all_ordered_parities(self):
  # |3/5|^2 + |4i/5|^2 = 1: connected by the closed universal core.
  # This benchmark uses the theorem, not a second copy of the verifier.
  for eps in [[0,0],[0,1],[1,0],[1,1]]:
   with self.subTest(orientation=eps):
    r=v.verify({'multipliers':[['3/5','0'],['0','4/5']],'orientation':eps},6,3000)
    self.assertFalse(r['disconnection_certified'])
 def test_independent_signed_real_contact_all_signs_and_parities(self):
  # Opposite orientations agree on R. Every sign choice with ratios
  # 1/3 and 2/3 gives a connected interval with one-point first contact.
  for sign_a in [-1,1]:
   for sign_b in [-1,1]:
    for eps in [[0,0],[0,1],[1,0],[1,1]]:
     with self.subTest(signs=(sign_a,sign_b),orientation=eps):
      r=v.verify({'multipliers':[[str(Q(sign_a,3)),'0'],[str(Q(2*sign_b,3)),'0']],'orientation':eps},6,3000)
      self.assertFalse(r['disconnection_certified'])
 def test_independent_oo_exact_depth_two_benchmark(self):
  # C2 exact example: sqrt(34)/4 > 9/[2(8-3sqrt(2))].
  # The first disks overlap, whereas total-depth-two disks separate.
  r=v.verify({'multipliers':[['3/8','3/8'],['3/8','3/8']],'orientation':[1,1]},3,1000)
  self.assertTrue(r['disconnection_certified'])
  self.assertEqual(r['depth_reached'],2)
  self.assertGreater(Q(r['minimum_squared_margin']),0)
 def certificate(self):
  return v.verify({'multipliers':[['3/8','3/8'],['3/8','3/8']],'orientation':[1,1]},3,1000)
 def test_independent_saved_certificate_checker_all_parities(self):
  for e in [[0,0],[0,1],[1,0],[1,1]]:
   r=v.verify({'multipliers':[['0','3/5'],['0','3/5']],'orientation':e},4,1000)
   self.assertEqual(checker.check(r)['status'],'verified')
 def test_incomplete_cover_is_rejected(self):
  r=self.certificate();r['separation_cover']['leaves'].pop()
  with self.assertRaisesRegex(ValueError,'incomplete'):checker.check(r)
 def test_duplicate_leaf_is_rejected(self):
  r=self.certificate();r['separation_cover']['leaves'].append(r['separation_cover']['leaves'][0])
  with self.assertRaisesRegex(ValueError,'Duplicate'):checker.check(r)
 def test_forged_gap_is_rejected(self):
  r=self.certificate();r['contact_gap_lower_bound']='100'
  with self.assertRaisesRegex(ValueError,'gap fails'):checker.check(r)
 def test_forged_neighbourhood_is_rejected(self):
  r=self.certificate();r['parameter_neighbourhood']['closed_radius']='1/10'
  with self.assertRaisesRegex(ValueError,'gap is not justified'):checker.check(r)
 def test_forged_input_hash_is_rejected(self):
  r=self.certificate();r['input']['multipliers'][0][0]='1/3'
  with self.assertRaisesRegex(ValueError,'checksum'):checker.check(r)
 def test_inward_modulus_bound_is_rejected(self):
  r=self.certificate();r['norm_upper_bounds'][0]='1/4'
  with self.assertRaisesRegex(ValueError,'modulus bounds'):checker.check(r)
 def test_unresolved_has_no_certificate(self):
  r=v.verify({'multipliers':[['3/5','0'],['3/5','0']],'orientation':[0,0]},4,1000)
  self.assertNotIn('separation_cover',r);self.assertNotIn('parameter_neighbourhood',r)
 def test_real_neighbourhood_with_analytic_oracle(self):
  r=v.verify({'multipliers':[['1/4','0'],['1/4','0']],'orientation':[0,0]})
  eta=Q(r['parameter_neighbourhood']['closed_radius'])
  self.assertGreater(eta,0);self.assertLess(2*(Q(1,4)+eta),1)
  self.assertEqual(checker.check(r)['status'],'verified')
 def test_mixed_depth_complete_cover(self):
  r=v.verify({'multipliers':[['0','-7/10'],['3/10','1/10']],'orientation':[0,0]},5,1000)
  self.assertEqual(set(len(u) for u,_ in r['separation_cover']['leaves']),{2,3,4})
  self.assertEqual(checker.check(r)['checked_leaves'],10)
 def test_cover_ancestor_conflict_is_rejected(self):
  r=self.certificate();r['separation_cover']['leaves'].append(['-','+'])
  with self.assertRaisesRegex(ValueError,'ancestor'):checker.check(r)

if __name__=='__main__':unittest.main()
