"""Exact coefficient checks for the two recovered mirror contact identities.

This uses only integer polynomial arithmetic in the independent variables c,d.
Substituting d=conjugate(c) afterward gives the mirror-family formulas. No
floating-point sample or symbolic package is needed.
"""
from dataclasses import dataclass
from pathlib import Path
import json


def clean(p):
    return {k: v for k, v in p.items() if v}


def add(p, q):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + v
    return clean(r)


def mul(p, q):
    r = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            key = (i + k, j + l)
            r[key] = r.get(key, 0) + a * b
    return clean(r)


@dataclass
class Rational:
    numerator: dict
    denominator: dict

    def __add__(self, other):
        return Rational(add(mul(self.numerator, other.denominator),
                            mul(other.numerator, self.denominator)),
                        mul(self.denominator, other.denominator))

    def __neg__(self):
        return Rational({k: -v for k, v in self.numerator.items()}, self.denominator)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return Rational(mul(self.numerator, other.numerator),
                        mul(self.denominator, other.denominator))

    def __truediv__(self, other):
        return Rational(mul(self.numerator, other.denominator),
                        mul(self.denominator, other.numerator))

    def equals(self, other):
        return mul(self.numerator, other.denominator) == mul(other.numerator, self.denominator)


ONE = Rational({(0, 0): 1}, {(0, 0): 1})
ZERO = Rational({}, {(0, 0): 1})
C = Rational({(1, 0): 1}, {(0, 0): 1})
D = Rational({(0, 1): 1}, {(0, 0): 1})


def word_map(word):
    linear, translation = ONE, ZERO
    for symbol in word:
        translation = translation + linear
        linear = linear * ({"1": C, "2": D}[symbol])
    return linear, translation


def periodic_value(prefix, period):
    a, b = word_map(prefix)
    u, v = word_map(period)
    return b + a * v / (ONE - u)


def main():
    first = periodic_value("1", "2") - periodic_value("2", "1")
    first_expected = (C - D) * (ONE - C - D) / ((ONE - C) * (ONE - D))
    second = periodic_value("122", "21") - periodic_value("211", "12")
    second_expected = (C - D) * (ONE - C * D - C * D - C * D * (C + D)) / (ONE - C * D)
    checks = [
        {"addresses": ["12^infinity", "21^infinity"],
         "expected": "(c-d)(1-c-d)/((1-c)(1-d))", "passed": first.equals(first_expected)},
        {"addresses": ["122(21)^infinity", "211(12)^infinity"],
         "expected": "(c-d)(1-2cd-cd(c+d))/(1-cd)", "passed": second.equals(second_expected)},
    ]
    assert all(item["passed"] for item in checks)
    result = {"method": "exact integer bivariate polynomial cross-multiplication", "checks": checks}
    Path(__file__).with_name("contact_identity_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
