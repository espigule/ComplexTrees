#!/usr/bin/env python3
"""Small exact finite-exclusion experiment for two planar similarities.

F(z)=-1+a J_p(z), G(z)=1+b J_q(z), with J_0(z)=z and J_1(z)=conj(z).
All numerical arithmetic uses fractions.Fraction. Finite exhaustion proves
DISCONNECTED. A surviving depth/budget frontier means UNRESOLVED only.

This is a research prototype and does not compute a connectedness locus.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import json


@dataclass(frozen=True)
class CQ:
    re: Q = Q(0)
    im: Q = Q(0)

    def __add__(self, other: CQ) -> CQ:
        return CQ(self.re + other.re, self.im + other.im)

    def __sub__(self, other: CQ) -> CQ:
        return CQ(self.re - other.re, self.im - other.im)

    def __mul__(self, other: CQ) -> CQ:
        return CQ(self.re * other.re - self.im * other.im,
                  self.re * other.im + self.im * other.re)

    def conj(self) -> CQ:
        return CQ(self.re, -self.im)

    def norm2(self) -> Q:
        return self.re * self.re + self.im * self.im

    def as_json(self) -> list[str]:
        return [str(self.re), str(self.im)]


@dataclass(frozen=True)
class Affine:
    """z -> B + A J_parity(z)."""
    A: CQ
    B: CQ
    parity: int

    def apply(self, z: CQ) -> CQ:
        return self.B + self.A * (z.conj() if self.parity else z)

    def after(self, other: Affine) -> Affine:
        """self composed with other (the new letter is appended on the right)."""
        a = other.A.conj() if self.parity else other.A
        b = other.B.conj() if self.parity else other.B
        return Affine(self.A * a, self.B + self.A * b, self.parity ^ other.parity)


IDENTITY = Affine(CQ(Q(1)), CQ(), 0)
BENCHMARKS = {
    "quarter_pp": ("1/4", "0", "1/4", "0", 0, 0,
                   "Separated real similarities; first-level certificate expected."),
    "imaginary_two_thirds_pp": ("0", "2/3", "0", "2/3", 0, 0,
                                "Both maps preserve orientation; nontrivial refinement."),
    "imaginary_two_thirds_pm": ("0", "2/3", "0", "2/3", 0, 1,
                                "Mixed orientation; nontrivial conjugation in cylinders."),
    "imaginary_two_thirds_mm": ("0", "2/3", "0", "2/3", 1, 1,
                                "Both maps reverse orientation."),
    "interval_half_pp": ("1/2", "0", "1/2", "0", 0, 0,
                         "The attractor is exactly [-2,2]; finite survival must remain unresolved."),
    "interval_half_pm": ("1/2", "0", "1/2", "0", 0, 1,
                         "Conjugation fixes the real interval [-2,2]; connected by an external exact argument."),
    "interval_half_mm": ("1/2", "0", "1/2", "0", 1, 1,
                         "Both conjugations fix [-2,2]; finite survival is not a certificate of connectedness."),
}


def contraction_bound(a: CQ, b: CQ, bits: int = 16) -> Q:
    """An exact rational upper bound q_bound<1, using integer square tests only."""
    m = max(a.norm2(), b.norm2())
    if m >= 1:
        raise ValueError("Both multipliers must have squared modulus strictly below 1")
    numerator_root, denominator_root = isqrt(m.numerator), isqrt(m.denominator)
    if numerator_root**2 == m.numerator and denominator_root**2 == m.denominator:
        return Q(numerator_root, denominator_root)
    # floor(sqrt(floor(x))) = floor(sqrt(x)); the following inequality then
    # decides exactly whether a ceiling step is needed.
    while True:
        scale = 1 << bits
        n = isqrt((m.numerator * scale * scale) // m.denominator)
        if n * n * m.denominator < m.numerator * scale * scale:
            n += 1
        bound = Q(n, scale)
        if bound < 1:
            assert bound * bound >= m
            return bound
        bits += 1


def maps(a: CQ, b: CQ, p: int, q: int) -> dict[str, Affine]:
    if p not in (0, 1) or q not in (0, 1):
        raise ValueError("Orientation flags must be 0 (identity) or 1 (conjugation)")
    return {"F": Affine(a, CQ(Q(-1)), p), "G": Affine(b, CQ(Q(1)), q)}


def compose_word(word: str, generators: dict[str, Affine]) -> Affine:
    result = IDENTITY
    for letter in word:
        result = result.after(generators[letter])
    return result


def certify(a: CQ, b: CQ, p: int, q: int, *, max_depth: int = 8,
            max_pairs: int = 20_000, benchmark: str | None = None) -> dict:
    if max_depth < 1 or max_pairs < 1:
        raise ValueError("max_depth and max_pairs must be positive")
    generators = maps(a, b, p, q)
    q_bound = contraction_bound(a, b)
    radius = 1 / (1 - q_bound)
    queue = deque([("F", "G", generators["F"], generators["G"])])
    excluded, frontier, counts = [], [], Counter()
    maximum_depth = 1
    while queue:
        u, v, left, right = queue.popleft()
        depth = len(u)
        maximum_depth = max(maximum_depth, depth)
        if counts["tested_pairs"] >= max_pairs:
            frontier.append({"u": u, "v": v, "reason": "pair_budget"})
            for u, v, left, right in queue:
                maximum_depth = max(maximum_depth, len(u))
                frontier.append({"u": u, "v": v, "reason": "pair_budget"})
            queue.clear()
            break
        counts["tested_pairs"] += 1
        distance2 = (left.B - right.B).norm2()
        threshold2 = (2 * q_bound**depth * radius)**2
        if distance2 > threshold2:
            excluded.append({"u": u, "v": v, "centre_u": left.B.as_json(),
                             "centre_v": right.B.as_json(), "squared_distance": str(distance2),
                             "squared_threshold": str(threshold2),
                             "strict_gap": str(distance2 - threshold2)})
        elif depth >= max_depth:
            frontier.append({"u": u, "v": v, "reason": "depth_limit"})
        else:
            counts["split_pairs"] += 1
            for letter_u in ("F", "G"):
                for letter_v in ("F", "G"):
                    queue.append((u + letter_u, v + letter_v,
                                  left.after(generators[letter_u]), right.after(generators[letter_v])))
    result = "UNRESOLVED" if frontier else "DISCONNECTED"
    return {
        "schema": "rational-planar-ifs-exclusion-1",
        "method": "strict disk separation with equal-length paired cylinder refinement",
        "benchmark": benchmark,
        "family": {"a": a.as_json(), "b": b.as_json(), "p": p, "q": q,
                   "convention": "J_0(z)=z; J_1(z)=conjugate(z); words compose left to right as maps"},
        "bound": {"contraction_bound": str(q_bound), "invariant_disk_radius": str(radius),
                  "invariant_disk_centre": ["0", "0"]},
        "limits": {"max_depth": max_depth, "max_pairs": max_pairs},
        "result": result,
        "statistics": {"tested_pairs": counts["tested_pairs"], "split_pairs": counts["split_pairs"],
                       "excluded_leaf_pairs": len(excluded), "frontier_leaf_pairs": len(frontier),
                       "maximum_depth": maximum_depth},
        "excluded": excluded,
        "frontier": frontier,
        "interpretation": ("Every F/G cylinder pair has been excluded by a strict rational inequality; the attractor is disconnected."
                           if not frontier else "A finite frontier survives the specified limits. No connectedness conclusion is made."),
    }


def benchmark(name: str, **limits) -> dict:
    ar, ai, br, bi, p, q, description = BENCHMARKS[name]
    result = certify(CQ(Q(ar), Q(ai)), CQ(Q(br), Q(bi)), p, q, benchmark=name, **limits)
    result["benchmark_description"] = description
    return result


def save(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="List exact named benchmarks")
    all_parser = sub.add_parser("benchmarks", help="Compute all named benchmarks")
    all_parser.add_argument("--output-dir", type=Path, default=Path("certificates"))
    one_parser = sub.add_parser("benchmark", help="Compute one named benchmark")
    one_parser.add_argument("name", choices=BENCHMARKS)
    one_parser.add_argument("--output", type=Path, required=True)
    custom = sub.add_parser("run", help="Use explicitly supplied rational coefficients")
    for name in ("a-real", "a-imag", "b-real", "b-imag"):
        custom.add_argument("--" + name, type=Q, default=Q(0))
    custom.add_argument("--p", type=int, choices=(0, 1), default=0)
    custom.add_argument("--q", type=int, choices=(0, 1), default=0)
    custom.add_argument("--output", type=Path, required=True)
    for child in (all_parser, one_parser, custom):
        child.add_argument("--max-depth", type=int, default=8)
        child.add_argument("--max-pairs", type=int, default=20_000)
    args = parser.parse_args()
    if args.command == "list":
        for name, values in BENCHMARKS.items():
            print(name + ": " + values[-1])
        return
    limits = {"max_depth": args.max_depth, "max_pairs": args.max_pairs}
    if args.command == "benchmarks":
        table = []
        for name in BENCHMARKS:
            value = benchmark(name, **limits)
            save(args.output_dir / (name + ".json"), value)
            table.append({"benchmark": name, "result": value["result"], **value["statistics"]})
        save(args.output_dir / "benchmark_summary.json", {"limits": limits, "results": table})
        for row in table:
            print(f"{row['benchmark']}: {row['result']}; tested={row['tested_pairs']}, excluded={row['excluded_leaf_pairs']}, frontier={row['frontier_leaf_pairs']}, depth={row['maximum_depth']}")
    else:
        value = (benchmark(args.name, **limits) if args.command == "benchmark" else
                 certify(CQ(args.a_real, args.a_imag), CQ(args.b_real, args.b_imag), args.p, args.q, **limits))
        save(args.output, value)
        print(value["result"] + ": " + value["interpretation"])
        print(json.dumps(value["statistics"], sort_keys=True))


if __name__ == "__main__":
    main()
