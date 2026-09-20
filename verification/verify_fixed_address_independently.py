#!/usr/bin/env python3
"""Independent exact checks for the fixed-address port; Python standard library.

Two-variable polynomial arithmetic verifies rational contact identities after
cross multiplication. Gaussian-rational arithmetic checks semilinear composition
and inversion for every labelled orientation type. No manuscript code imported.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json


def padd(a, b):
    ans = dict(a)
    for m, x in b.items():
        ans[m] = ans.get(m, 0) + x
        if not ans[m]:
            del ans[m]
    return ans


def pneg(a):
    return {m: -x for m, x in a.items()}


def pmul(a, b):
    ans = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            m = (i + k, j + l)
            ans[m] = ans.get(m, 0) + x * y
    return {m: x for m, x in ans.items() if x}


class Rational:
    def __init__(self, numerator, denominator=None):
        self.n = numerator
        self.d = {(0, 0): 1} if denominator is None else denominator

    def __add__(self, b):
        return Rational(padd(pmul(self.n, b.d), pmul(b.n, self.d)),
                        pmul(self.d, b.d))

    def __neg__(self):
        return Rational(pneg(self.n), self.d)

    def __sub__(self, b):
        return self + (-b)

    def __mul__(self, b):
        return Rational(pmul(self.n, b.n), pmul(self.d, b.d))

    def __truediv__(self, b):
        assert b.n
        return Rational(pmul(self.n, b.d), pmul(self.d, b.n))


zero = Rational({})
one = Rational({(0, 0): 1})
c = Rational({(1, 0): 1})
d = Rational({(0, 1): 1})


def word_map(word):
    translation, multiplier = zero, one
    for letter in word:
        translation = translation + multiplier
        multiplier = multiplier * {1: c, 2: d}[letter]
    return translation, multiplier


def periodic_value(prefix, period):
    tp, rp = word_map(prefix)
    tv, rv = word_map(period)
    return tp + rp * tv / (one - rv)


cases = [
    ((1,), (2,), (2,), (1,)),
    ((1, 2, 2), (2, 1), (2, 1, 1), (1, 2)),
]
expected = [
    (c - d) * (one - c - d) / ((one - c) * (one - d)),
    (c - d) * (one - c * d - c * d - c * d * (c + d)) / (one - c * d),
]
for (p1, v1, p2, v2), claimed in zip(cases, expected):
    actual = periodic_value(p1, v1) - periodic_value(p2, v2)
    assert not (actual - claimed).n

# Exact differentiation of 2(x^2+y^2)(1+x)-1 in its second variable.
curve = {(2, 0): 2, (0, 2): 2, (3, 0): 2, (1, 2): 2, (0, 0): -1}
derivative_y = {(i, j-1): j*v for (i, j), v in curve.items() if j}
assert derivative_y == {(0, 1): 4, (1, 1): 4}

Z = (F(0), F(0))
ONE = (F(1), F(0))


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def neg(z):
    return -z[0], -z[1]


def conj(z):
    return z[0], -z[1]


def mul(z, w):
    return z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0]


def inv(z):
    n = z[0]**2 + z[1]**2
    return z[0]/n, -z[1]/n


def J(z, parity):
    return conj(z) if parity else z


def apply(g, z):
    translation, multiplier, parity = g
    return add(translation, mul(multiplier, J(z, parity)))


def compose(g, h):
    t, r, e = g
    u, v, f = h
    return add(t, mul(r, J(u, e))), mul(r, J(v, e)), e ^ f


def invert(g, z):
    t, r, e = g
    return J(mul(add(z, neg(t)), inv(r)), e)


a = (F(1, 3), F(1, 4))
b = (F(-2, 5), F(1, 7))
checks = 0
for e, f in product((0, 1), repeat=2):
    generators = [(neg(ONE), a, e), (ONE, b, f)]
    for n in range(1, 6):
        for u in product((0, 1), repeat=n):
            g = (Z, ONE, 0)
            for letter in u:
                g = compose(g, generators[letter])
            for z in [Z, ONE, (F(2, 7), F(-3, 5))]:
                nested = z
                for letter in reversed(u):
                    nested = apply(generators[letter], nested)
                assert apply(g, z) == nested
                assert invert(g, nested) == z
                assert apply(g, invert(g, z)) == z
                checks += 1

record = {
    "method": "Independent polynomial cross multiplication and Gaussian-rational arithmetic; no port code imported.",
    "mirror_identity_residuals": [0, 0],
    "mirror_addresses": ["12^infinity versus 21^infinity", "122(21)^infinity versus 211(12)^infinity"],
    "upper_branch_y_derivative": "4*y*(1+x)",
    "semilinear_exact_checks": checks,
    "word_lengths": [1, 5],
    "labelled_orientation_types": 4,
    "test_parameters": {"a": ["1/3", "1/4"], "b": ["-2/5", "1/7"]},
    "scope": "These algebraic checks supplement the proof review; they are not a numerical membership or global-priority certificate.",
    "result": "PASS",
}
output = Path(__file__).with_name("fixed_address_independent_checks.json")
output.write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
