#!/usr/bin/env python3
"""Independent exact consumer of adaptive-zipper-degree-capture-v1 records.

No producer modules are imported. All geometric decisions use Fraction.
The finite cover is checked by its dyadic interval partition; no search or
adaptive generation is repeated. Winding is counted across a different ray
from the producer (the negative imaginary ray).
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import platform
import time

ZERO = (F(0), F(0))
ONE = (F(1), F(0))
HALF = (F(1, 2), F(0))


def plus(z, w):
    return (z[0]+w[0], z[1]+w[1])


def minus(z, w):
    return (z[0]-w[0], z[1]-w[1])


def times(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])


def conjugate(z):
    return (z[0], -z[1])


def scalar(a, z):
    return (a*z[0], a*z[1])


def square_norm(z):
    return z[0]*z[0]+z[1]*z[1]


def determinant(z, w):
    return z[0]*w[1]-z[1]*w[0]


def evaluate(operator, z):
    shift, coefficient, reversal = operator
    return plus(shift, times(coefficient, conjugate(z) if reversal else z))


def compose(outer, inner):
    t, a, k = outer
    s, b, l = inner
    return (evaluate(outer, s), times(a, conjugate(b) if k else b), k ^ l)


def canonical_maps(family, q):
    if len(family) != 4 or family[:2] not in ('DD', 'DO', 'OO') or family[2:] not in ('00', '01', '10', '11'):
        raise ValueError('Unsupported marked binary zipper family')
    reversals = {'DD': (0, 0), 'DO': (0, 1), 'OO': (1, 1)}[family[:2]]
    bits = tuple(int(s) for s in family[2:])
    first = (scalar(bits[0], q), scalar((-1)**bits[0], q), reversals[0])
    complement = minus(ONE, q)
    second = (plus(q, scalar(bits[1], complement)), scalar((-1)**bits[1], complement), reversals[1])
    return (first, second), bits


def interval_map(bits, word):
    offset, slope = F(0), F(1)
    for char in word:
        digit = int(char)
        offset += slope*F(digit+bits[digit], 2)
        slope *= F((-1)**bits[digit], 2)
    return offset, slope


def interval(bits, word):
    offset, slope = interval_map(bits, word)
    return min(offset, offset+slope), max(offset, offset+slope)


def closest_squared_to_segment(a, b):
    delta = minus(b, a)
    length2 = square_norm(delta)
    if not length2:
        return square_norm(a)
    projection = -(a[0]*delta[0]+a[1]*delta[1])/length2
    projection = max(F(0), min(F(1), projection))
    return square_norm(plus(a, scalar(projection, delta)))


def winding_negative_imaginary_ray(vertices):
    """Apply the half-open crossing rule after rotation by +i.

    The positive real ray in the rotated plane is the negative imaginary ray
    in the original plane. Exact zero-clearance checks precede this count.
    """
    result = 0
    for a, b in zip(vertices, vertices[1:]):
        cross = determinant(a, b)
        if a[0] <= 0 < b[0] and cross > 0:
            result += 1
        elif b[0] <= 0 < a[0] and cross < 0:
            result -= 1
    return result


def verify(path):
    started = time.perf_counter()
    payload = path.read_bytes()
    doc = json.loads(payload)
    assert doc['schema'] == 'adaptive-zipper-degree-capture-v1'
    assert doc['status'] == 'proved_nonembedded_parameter_disk'
    family = doc['family']
    q = (F(doc['q']['real']), F(doc['q']['imag']))
    radius = F(doc['parameter_disk_radius'])
    assert radius >= 0
    generators, bits = canonical_maps(family, q)
    u, v = doc['words']['u'], doc['words']['v']
    assert u and v and not (set(u+v)-set('01'))
    source = [interval(bits, u), interval(bits, v)]
    assert source[0][1] < source[1][0] or source[1][1] < source[0][0]
    assert source == [tuple(map(F, entry)) for entry in doc['source_intervals']]

    # Existence and endpoint compatibility of the canonical continuous curve.
    # h_i sends local endpoints e_i and 1-e_i to global i/2 and (i+1)/2.
    endpoint_targets = ((ZERO, q), (q, ONE))
    for i, (start, end) in enumerate(endpoint_targets):
        assert evaluate(generators[i], ONE if bits[i] else ZERO) == start
        assert evaluate(generators[i], ZERO if bits[i] else ONE) == end

    r0, r1 = map(F, doc['contraction_upper_bounds'])
    assert 0 < r0 < 1 and 0 < r1 < 1
    assert square_norm(q) <= r0*r0
    assert square_norm(minus(ONE, q)) <= r1*r1
    r = max(r0, r1)
    assert r+radius < 1
    R = F(doc['invariant_disk_radius_upper'])
    assert R >= F(1, 2)
    assert r0*R+r1/2 <= R and r1*R+r0/2 <= R
    for i, generator in enumerate(generators):
        # Both map-center distances are bounded by the opposite modulus/2.
        center = evaluate(generator, HALF)
        assert square_norm(minus(center, HALF)) <= ((r1 if i == 0 else r0)/2)**2
    motion = F(doc['uniform_boundary_motion_upper'])
    pointwise_curve_motion = (R+F(1, 2))*radius/(1-r-radius)
    assert motion >= 2*pointwise_curve_motion

    @lru_cache(None)
    def word_operator(word):
        if not word:
            return (ZERO, ONE, 0)
        return compose(word_operator(word[:-1]), generators[int(word[-1])])

    def endpoint(word, side):
        return evaluate(word_operator(word), ONE if side else ZERO)

    roots = (u, v, u, v)
    anchors = (endpoint(v, 0), endpoint(u, 1), endpoint(v, 1), endpoint(u, 0))
    signs = (1, -1, 1, -1)
    by_side = [[] for _ in range(4)]
    for record in doc['boundary_cover']:
        side, word = record['side'], record['word']
        assert type(side) is int and 0 <= side < 4
        assert word.startswith(roots[side]) and not (set(word)-set('01'))
        by_side[side].append(word)
    assert all(by_side)

    polygon = []
    squared_clearances = []
    side_reports = []
    all_cover_words = []
    for side, listed in enumerate(by_side):
        assert len(set(listed)) == len(listed)
        root = roots[side]
        leaves = []
        for word in listed:
            suffix = word[len(root):]
            left, right = interval(bits, suffix)
            leaves.append((left, right, word, suffix))
        leaves.sort()
        assert [leaf[2] for leaf in leaves] == listed
        cursor = F(0)
        for left, right, _, _ in leaves:
            assert left == cursor and right > left
            cursor = right
        assert cursor == 1

        # Verify exhaustiveness also through a complete binary prefix tree.
        suffixes = {entry[3] for entry in leaves}
        maximum_suffix = max(map(len, suffixes))
        tree_nodes = 0
        def complete(prefix):
            nonlocal tree_nodes
            tree_nodes += 1
            if prefix in suffixes:
                return
            assert len(prefix) < maximum_suffix
            complete(prefix+'0')
            complete(prefix+'1')
        complete('')
        assert tree_nodes == 2*len(leaves)-1

        segments = []
        for left, right, word, suffix in leaves:
            operator = word_operator(word)
            center = scalar(signs[side], minus(evaluate(operator, HALF), anchors[side]))
            scale = r0**word.count('0') * r1**word.count('1')
            assert square_norm(operator[1]) <= scale**2
            base_radius = R*scale
            total_radius = base_radius+motion
            clearance = square_norm(center)-total_radius**2
            assert clearance > 0
            squared_clearances.append(clearance)
            _, interval_slope = interval_map(bits, suffix)
            start_bit = 0 if interval_slope > 0 else 1
            a = scalar(signs[side], minus(endpoint(word, start_bit), anchors[side]))
            b = scalar(signs[side], minus(endpoint(word, 1-start_bit), anchors[side]))
            # Endpoints and chord are in the same certified disk as the true subarc.
            assert square_norm(minus(a, center)) <= base_radius**2
            assert square_norm(minus(b, center)) <= base_radius**2
            assert closest_squared_to_segment(a, b) > 0
            segments.append((a, b))
            all_cover_words.append(word)

        if side in (2, 3):
            segments = [(b, a) for a, b in reversed(segments)]
        for a, b in segments:
            if polygon:
                assert polygon[-1] == a
            else:
                polygon.append(a)
            polygon.append(b)
        side_reports.append({'side': side, 'root_word': root, 'leaves': len(leaves), 'tree_nodes': tree_nodes,
                             'exact_unit_interval_partition': True, 'endpoint_chain_verified': True})

    assert polygon[0] == polygon[-1]
    assert len(polygon)-1 == len(doc['boundary_cover']) == doc['accepted_cylinders']
    assert sum(row['tree_nodes'] for row in side_reports) == doc['visited_cylinders']
    assert max(map(len, all_cover_words)) == doc['maximum_word_depth']
    minimum = min(squared_clearances)
    assert minimum == F(doc['minimum_squared_disk_clearance'])
    degree = winding_negative_imaginary_ray(polygon)
    assert degree != 0 and degree == doc['degree']
    polygon_minimum = min(closest_squared_to_segment(a, b) for a, b in zip(polygon, polygon[1:]))
    return {'verified': True, 'family': family, 'q': doc['q'], 'parameter_disk_radius': str(radius),
            'certificate_sha256': hashlib.sha256(payload).hexdigest(),
            'source_intervals': doc['source_intervals'], 'words': doc['words'],
            'all_decisions_exact': True, 'degree': degree, 'winding_ray': 'negative imaginary ray; exact half-open crossing convention',
            'accepted_cylinders': len(doc['boundary_cover']), 'visited_cylinders': doc['visited_cylinders'],
            'maximum_word_depth': doc['maximum_word_depth'], 'sides': side_reports,
            'invariant_disk_radius_upper': str(R), 'uniform_boundary_motion_upper': str(motion),
            'minimum_squared_disk_clearance': str(minimum),
            'minimum_squared_disk_clearance_decimal': float(minimum),
            'minimum_polygon_distance_squared': str(polygon_minimum),
            'minimum_polygon_distance_squared_decimal': float(polygon_minimum),
            'scope': 'Canonical zipper noninjective throughout the stated closed parameter disk; no island assertion',
            'runtime_seconds': time.perf_counter()-started}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('certificates', nargs='*', type=Path)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('ADAPTIVE_DEGREE_INDEPENDENT_REPLAY.json'))
    args = parser.parse_args()
    paths = args.certificates or sorted(Path(__file__).with_name('certificates').glob('*_adaptive_capture.json'))
    assert paths
    results = [verify(path) for path in paths]
    report = {'audit_date': '2026-09-18', 'verified': all(r['verified'] for r in results),
              'implementation': 'Independent stdlib Fraction certificate consumer; no producer imports or numerical decisions',
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'python': platform.python_version(), 'certificates': results}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'verified': report['verified'], 'certificates': [
        {k: r[k] for k in ('family', 'parameter_disk_radius', 'degree', 'accepted_cylinders', 'maximum_word_depth', 'minimum_squared_disk_clearance_decimal', 'runtime_seconds')}
        for r in results], 'output': str(args.output.resolve())}, indent=2))


if __name__ == '__main__':
    main()
