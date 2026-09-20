#!/usr/bin/env python3
"""Scientific correctness and regression tests; standard library only."""
import copy
from fractions import Fraction as Q
from itertools import product
import unittest

from exact_exclusion import CQ, BENCHMARKS, benchmark, certify, compose_word, contraction_bound, maps
from verify_certificate import centre, verify


class ExactExclusionTests(unittest.TestCase):
    def test_contraction_bounds(self):
        for a, b in [(CQ(), CQ()), (CQ(Q(1, 2)), CQ(Q(1, 3))),
                     (CQ(Q(1, 3), Q(2, 5)), CQ(Q(-2, 7), Q(1, 4))),
                     (CQ(Q(999_999, 1_000_000)), CQ(Q(1, 9)))]:
            q_bound = contraction_bound(a, b)
            self.assertLess(q_bound, 1)
            self.assertGreaterEqual(q_bound*q_bound, a.norm2())
            self.assertGreaterEqual(q_bound*q_bound, b.norm2())
        with self.assertRaises(ValueError):
            contraction_bound(CQ(Q(1)), CQ())

    def test_composition_against_independent_right_to_left_application(self):
        a, b = CQ(Q(1, 3), Q(2, 5)), CQ(Q(-2, 7), Q(1, 4))
        points = (CQ(), CQ(Q(2, 7), Q(-3, 5)), CQ(Q(-1, 11), Q(5, 13)))
        for p, q in product((0, 1), repeat=2):
            generators = maps(a, b, p, q)
            for length in range(1, 6):
                for letters in product("FG", repeat=length):
                    word = "".join(letters)
                    composed = compose_word(word, generators)
                    independent_centre = centre(word, (a.re, a.im), (b.re, b.im), p, q)
                    self.assertEqual((composed.B.re, composed.B.im), independent_centre)
                    for point in points:
                        direct = point
                        for letter in reversed(word):
                            direct = generators[letter].apply(direct)
                        self.assertEqual(composed.apply(point), direct)
                    self.assertEqual(composed.parity, sum(p if letter == "F" else q for letter in word) % 2)

    def test_all_benchmarks_and_external_interval_fact(self):
        for name in BENCHMARKS:
            result = benchmark(name, max_depth=8, max_pairs=20_000)
            self.assertTrue(verify(result)["verified"])
            if name.startswith("interval_half"):
                self.assertEqual(result["result"], "UNRESOLVED")
                self.assertGreater(len(result["frontier"]), 0)
                # Exact external fact: images of [-2,2] are [-2,0] and [0,2].
                f, g = maps(CQ(Q(1, 2)), CQ(Q(1, 2)), result["family"]["p"], result["family"]["q"]).values()
                self.assertEqual([f.apply(CQ(Q(x))).re for x in (-2, 2)], [Q(-2), Q(0)])
                self.assertEqual([g.apply(CQ(Q(x))).re for x in (-2, 2)], [Q(0), Q(2)])
            else:
                self.assertEqual(result["result"], "DISCONNECTED")

    def test_budget_and_depth_never_claim_connectedness(self):
        for depth, budget in ((1, 20_000), (8, 1), (8, 2)):
            result = benchmark("imaginary_two_thirds_pm", max_depth=depth, max_pairs=budget)
            self.assertEqual(result["result"], "UNRESOLVED")
            self.assertTrue(verify(result)["verified"])

    def test_verifier_rejects_corruption(self):
        genuine = benchmark("imaginary_two_thirds_pm")
        for mutate in (
            lambda x: x["excluded"].pop(),
            lambda x: x["excluded"][0].update({"strict_gap": "0"}),
            lambda x: x["bound"].update({"contraction_bound": "1/10"}),
            lambda x: x["excluded"].append(copy.deepcopy(x["excluded"][0])),
            lambda x: x["excluded"][0].update({"centre_u": ["0", "0"]}),
        ):
            changed = copy.deepcopy(genuine)
            mutate(changed)
            with self.assertRaises(ValueError):
                verify(changed)
        partial = benchmark("interval_half_pm")
        partial["result"] = "CONNECTED"
        with self.assertRaises(ValueError):
            verify(partial)

    def test_zero_multipliers(self):
        result = certify(CQ(), CQ(), 1, 0)
        self.assertEqual(result["result"], "DISCONNECTED")
        self.assertTrue(verify(result)["verified"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
