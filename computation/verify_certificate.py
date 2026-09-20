#!/usr/bin/env python3
"""Standalone exact verifier for rational-planar-ifs-exclusion-1 JSON.

This file imports no producer code. Cylinder centres are reconstructed by
applying individual letters from right to left to zero. The verifier checks
every strict inequality and the complete four-child prefix-tree covering.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json


def read_complex(value):
    if not isinstance(value, list) or len(value) != 2 or not all(isinstance(x, str) for x in value):
        raise ValueError("Complex values must be pairs of rational strings")
    return Q(value[0]), Q(value[1])


def multiply(a, b):
    return a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0]


def centre(word, a, b, p, q):
    z = Q(0), Q(0)
    for letter in reversed(word):
        multiplier, parity, translation = (a, p, -1) if letter == "F" else (b, q, 1)
        if parity:
            z = z[0], -z[1]
        z = multiply(multiplier, z)
        z = z[0] + translation, z[1]
    return z


def verify(data):
    if data.get("schema") != "rational-planar-ifs-exclusion-1":
        raise ValueError("Unknown certificate schema")
    family, bound = data["family"], data["bound"]
    limits = data["limits"]
    if any(type(limits.get(key)) is not int or limits[key] < 1 for key in ("max_depth", "max_pairs")):
        raise ValueError("Invalid computation limits")
    a, b = read_complex(family["a"]), read_complex(family["b"])
    p, q = family["p"], family["q"]
    if type(p) is not int or type(q) is not int or p not in (0, 1) or q not in (0, 1):
        raise ValueError("Invalid orientation flag")
    q_bound, radius = Q(bound["contraction_bound"]), Q(bound["invariant_disk_radius"])
    if not 0 <= q_bound < 1:
        raise ValueError("Contraction bound must lie in [0,1)")
    if q_bound*q_bound < max(a[0]*a[0]+a[1]*a[1], b[0]*b[0]+b[1]*b[1]):
        raise ValueError("The rational contraction bound is too small")
    if radius != 1 / (1-q_bound) or read_complex(bound["invariant_disk_centre"]) != (Q(0), Q(0)):
        raise ValueError("Invalid invariant disk")
    leaves = {}
    maximum_depth = 0
    for kind in ("excluded", "frontier"):
        for item in data[kind]:
            u, v = item["u"], item["v"]
            if not isinstance(u, str) or not isinstance(v, str) or not u or len(u) != len(v):
                raise ValueError("Cylinder words must have the same positive length")
            if u[0] != "F" or v[0] != "G" or set(u+v) - {"F", "G"}:
                raise ValueError("Invalid F/G root cylinder pair")
            if (u, v) in leaves:
                raise ValueError("Duplicate leaf")
            leaves[(u, v)] = kind
            maximum_depth = max(maximum_depth, len(u))
            if len(u) > limits["max_depth"]:
                raise ValueError("Leaf lies beyond the recorded depth limit")
            if kind == "excluded":
                left, right = centre(u, a, b, p, q), centre(v, a, b, p, q)
                distance2 = (left[0]-right[0])**2 + (left[1]-right[1])**2
                threshold2 = (2*q_bound**len(u)*radius)**2
                if not distance2 > threshold2:
                    raise ValueError("A claimed strict exclusion fails")
                if left != read_complex(item["centre_u"]) or right != read_complex(item["centre_v"]):
                    raise ValueError("Recorded cylinder centre is incorrect")
                if Q(item["squared_distance"]) != distance2 or Q(item["squared_threshold"]) != threshold2:
                    raise ValueError("Recorded squared-distance data are incorrect")
                if Q(item["strict_gap"]) != distance2-threshold2:
                    raise ValueError("Recorded strict gap is incorrect")
            else:
                if item.get("reason") not in ("depth_limit", "pair_budget"):
                    raise ValueError("Unknown frontier reason")
                if item["reason"] == "depth_limit":
                    if len(u) != limits["max_depth"]:
                        raise ValueError("Depth frontier is inconsistent with the recorded limit")
                    left, right = centre(u, a, b, p, q), centre(v, a, b, p, q)
                    if (left[0]-right[0])**2 + (left[1]-right[1])**2 > (2*q_bound**len(u)*radius)**2:
                        raise ValueError("A tested depth frontier should have been excluded")
    if not leaves:
        raise ValueError("Empty cover")
    # Every non-leaf pair must have all four refined children covered. A
    # prefix set avoids exponential traversal through branches missing in the
    # submitted cover. Leftover leaves detect an ancestor/descendant overlap.
    prefixes = set()
    for u, v in leaves:
        for length in range(1, len(u)):
            prefixes.add((u[:length], v[:length]))
    pending, used, split_count = [("F", "G")], set(), 0
    while pending:
        pair = pending.pop()
        if pair in leaves:
            used.add(pair)
        elif pair in prefixes:
            split_count += 1
            u, v = pair
            pending.extend((u+s, v+t) for s in ("F", "G") for t in ("F", "G"))
        else:
            raise ValueError("The leaf list does not cover the complete F/G cylinder-pair tree")
    if len(used) != len(leaves):
        raise ValueError("Redundant descendant leaf under an existing leaf")
    expected = "UNRESOLVED" if data["frontier"] else "DISCONNECTED"
    if data["result"] != expected:
        raise ValueError("Result is inconsistent with the surviving frontier")
    statistics = data["statistics"]
    if statistics["split_pairs"] != split_count or statistics["excluded_leaf_pairs"] != len(data["excluded"]) or statistics["frontier_leaf_pairs"] != len(data["frontier"]):
        raise ValueError("Recorded tree statistics do not match the verified tree")
    if statistics["maximum_depth"] != maximum_depth:
        raise ValueError("Recorded maximum depth is incorrect")
    tested = split_count + len(data["excluded"]) + sum(item["reason"] == "depth_limit" for item in data["frontier"])
    if statistics["tested_pairs"] != tested or tested > limits["max_pairs"]:
        raise ValueError("Recorded number of tested pairs is inconsistent")
    if any(item["reason"] == "pair_budget" for item in data["frontier"]) and tested != limits["max_pairs"]:
        raise ValueError("Budget frontier is inconsistent with the recorded limit")
    return {"verified": True, "result": expected, "excluded_leaves": len(data["excluded"]),
            "frontier_leaves": len(data["frontier"]), "maximum_depth": maximum_depth,
            "meaning": ("The exact finite cover proves disconnectedness." if expected == "DISCONNECTED" else
                        "The partial exclusion cover is valid; no connectedness conclusion follows.")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path, nargs="+")
    args = parser.parse_args()
    for path in args.certificate:
        result = verify(json.loads(path.read_text(encoding="utf-8")))
        print(path.name + ": " + json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
