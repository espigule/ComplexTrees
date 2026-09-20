"""Exact finite-prefix checks for the restored address recodings.

This is a verification aid for the algebra in the proofs, not an approximation
certificate for infinite parameter-space topology. Fractions avoid roundoff.
"""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z):
    return z[0], -z[1]


def scale(z, scalar):
    return scalar*z[0], scalar*z[1]


def divide(z, w):
    return scale(mul(z, conj(w)), 1/(w[0]*w[0]+w[1]*w[1]))


def code(word, a, family):
    z = (F(0), F(0))
    for e in reversed(word):
        use_conjugate = family in ("OO", "anti") or (family == "DO" and e == 1)
        v = conj(z) if use_conjugate else z
        coefficient = (e * a[0], e * a[1]) if family == "anti" else a
        z = add((F(e), F(0)), mul(coefficient, v))
    return z


def dd_recode(word):
    return [((-1) ** (j // 2)) * e for j, e in enumerate(word)]


def mixed_recode(word):
    state = [1, 1]
    answer = []
    for j, e in enumerate(word):
        state[j % 2] *= e
        answer.append(state[j % 2])
    return answer


def anti_recode(word):
    p = 1
    answer = []
    for e in word:
        p *= e
        answer.append(p)
    return answer


def main():
    checked = 0
    for depth in range(1, 11):
        for word in product((-1, 1), repeat=depth):
            for r in (F(1, 3), F(3, 5), F(4, 5)):
                a = (F(0), r)
                for family, recode in (("DD", dd_recode), ("DO", mixed_recode), ("anti", anti_recode)):
                    mapped = recode(word)
                    assert mapped[0] == word[0]
                    assert code(word, a, family) == code(mapped, a, "OO"), (word, r, family)
                    checked += 1
            a = (F(3, 10), F(2, 5))
            assert code(word, a, "anti") == code(anti_recode(word), a, "OO")
            checked += 1

    def poly_add(*terms):
        answer = {}
        for sign, term in terms:
            for degree, coefficient in term.items():
                answer[degree] = answer.get(degree, 0) + sign * coefficient
        return {degree: coefficient for degree, coefficient in answer.items() if coefficient}

    def poly_mul(left, right):
        answer = {}
        for (i, j), c in left.items():
            for (k, ell), d in right.items():
                degree = i+k, j+ell
                answer[degree] = answer.get(degree, 0) + c*d
        return {degree: coefficient for degree, coefficient in answer.items() if coefficient}

    one = {(0, 0): 1}
    a, b = {(1, 0): 1}, {(0, 1): 1}
    one_a, one_b = poly_add((1, one), (-1, a)), poly_add((1, one), (-1, b))
    identities = []
    for k in range(1, 9):
        ak, bk = {(k, 0): 1}, {(0, k): 1}
        left_numerator = poly_add((-1, poly_mul(poly_add((1, one), (-1, ak)), one_b)), (1, poly_mul(ak, one_a)))
        right_numerator = poly_add((1, poly_mul(poly_add((1, one), (-1, bk)), one_a)), (-1, poly_mul(bk, one_b)))
        target_numerator = poly_mul(poly_add((1, a), (1, b), (-2, one)), poly_add((1, ak), (1, bk), (-1, one)))
        assert not poly_add((1, left_numerator), (-1, right_numerator), (1, target_numerator))
        identities.append(k)

    triangles = []
    zero, unit = (F(0), F(0)), (F(1), F(0))
    for adjacent, opposite, hypotenuse in ((3, 4, 5), (5, 12, 13), (8, 15, 17)):
        cosine, sine = F(adjacent, hypotenuse), F(opposite, hypotenuse)
        phase = cosine, sine
        c = scale(phase, cosine)
        one_c = add(unit, scale(c, -1))
        foot = c[0], F(0)
        raw0 = lambda z: mul(c, conj(z))
        raw1 = lambda z: add(c, mul(one_c, conj(z)))
        assert tuple(map(raw0, (zero, unit, c))) == (zero, c, foot)
        assert tuple(map(raw1, (zero, unit, c))) == (c, unit, foot)
        mean = c[0], c[1]/3
        d = add(raw1(mean), scale(mean, -1))
        multiplier_phase = divide(conj(d), d)
        normalized_a = mul(c, multiplier_phase)
        normalized_b = mul(one_c, multiplier_phase)
        phase_minus3 = mul(mul(conj(phase), conj(phase)), conj(phase))
        assert normalized_a == scale(phase_minus3, -cosine)
        assert normalized_b == mul((F(0), sine), phase_minus3)
        # Directly certify vertex projections for the dust path at a rational t.
        t = F(9, 10)
        vertices0 = [scale(raw0(z), t) for z in (zero, unit, c)]
        vertices1 = [add(c, scale(add(raw1(z), scale(c, -1)), t)) for z in (zero, unit, c)]
        assert max(z[0] for z in vertices0) == t*c[0]
        assert min(z[0] for z in vertices1) == c[0]
        assert t*c[0] < c[0]
        triangles.append(f"cos(theta)={cosine}, sin(theta)={sine}")

    result = {
        "status": "passed",
        "finite_prefix_exact_equalities": checked,
        "maximum_prefix_length": 10,
        "quarter_turn_rational_radii": ["1/3", "3/5", "4/5"],
        "OO_nonreal_rational_coefficient": "3/10+2i/5",
        "cross_contact_symbolic_factorization_k": identities,
        "triangle_exact_vertex_normalization_and_separation_checks": triangles,
        "scope": "Exact finite-prefix and symbolic identity checks complement the written infinite-series proofs; no global numerical frontier claim is inferred."
    }
    out = Path(__file__).with_name("exact_family_checks.json")
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
