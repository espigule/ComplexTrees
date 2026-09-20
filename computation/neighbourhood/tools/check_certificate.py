#!/usr/bin/env python3
"""Check a finite separation cover without importing or rerunning its search.

Only the Python standard library is required. Complex similarities are evaluated
here as rational 2-by-2 real matrices, independently of the search implementation.
The checker verifies cover completeness, all geometric inequalities, and the
advertised closed parameter neighbourhood. It is not a proof-assistant kernel.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json


def rational(value):
    if not isinstance(value, str) or len(value) > 160000:
        raise ValueError('Expected a bounded exact rational string.')
    result = Q(value)
    if max(result.numerator.bit_length(), result.denominator.bit_length()) > 262144:
        raise ValueError('Rational value exceeds checker resource bound.')
    return result


def require(condition, message):
    if not condition:
        raise ValueError(message)


def matrix_product(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in (0, 1))
                       for j in (0, 1)) for i in (0, 1))


def matrix_vector(a, b):
    return tuple(sum(a[i][k] * b[k] for k in (0, 1)) for i in (0, 1))


def check(report):
    require(report.get('schema') == 'rational-disconnection/2', 'Unsupported certificate schema.')
    spec = report['input']
    eps = spec['orientation']
    require(isinstance(eps, list) and len(eps) == 2 and
            all(type(x) is int and x in (0, 1) for x in eps), 'Invalid orientations.')
    raw = spec['multipliers']
    require(isinstance(raw, list) and len(raw) == 2 and
            all(isinstance(z, list) and len(z) == 2 for z in raw), 'Invalid multipliers.')
    values = [tuple(map(rational, z)) for z in raw]
    norms = [sum(x*x for x in z) for z in values]
    require(all(0 < x < 1 for x in norms), 'Nonzero strict contractions required.')
    canonical = {'multipliers': [[str(x) for x in z] for z in values], 'orientation': eps}
    digest = hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    require(report['input_sha256'] == digest, 'Input checksum does not match rational coefficients.')
    bounds = list(map(rational, report['norm_upper_bounds']))
    require(len(bounds) == 2 and all(0 < r < 1 and r*r >= n for r, n in zip(bounds, norms)),
            'Invalid outward modulus bounds.')
    q = max(bounds)
    radius = rational(report['enclosing_radius'])
    require(radius >= 1/(1-q), 'Containing disk is not invariant.')
    gap = rational(report['contact_gap_lower_bound'])
    require(gap > 0, 'A positive contact gap is required.')
    squared_margin = rational(report['minimum_squared_margin'])
    require(squared_margin > 0, 'A positive squared margin is required.')
    cover = report['separation_cover']
    require(cover.get('alphabet') == '-+' and cover.get('composition') == 'outermost first',
            'Unknown word convention.')
    leaves = cover['leaves']
    require(isinstance(leaves, list) and 0 < len(leaves) <= 1000000, 'Invalid cover size.')

    # Words with first signs - and + are a quaternary tree of paired suffixes.
    # A complete finite prefix antichain covers every infinite address pair.
    trie = {}
    for pair in leaves:
        require(isinstance(pair, list) and len(pair) == 2 and
                all(isinstance(w, str) for w in pair), 'Invalid word pair.')
        u, v = pair
        require(1 <= len(u) == len(v) <= 30 and u[0] == '-' and v[0] == '+' and
                all(c in '-+' for c in u+v), 'Invalid first-level pair or word depth.')
        node = trie
        for a, b in zip(u[1:], v[1:]):
            require('leaf' not in node, 'Cover contains an ancestor and its descendant.')
            digit = 2 * (a == '+') + (b == '+')
            node = node.setdefault(digit, {})
        require(not node, 'Duplicate leaf or ancestor in cover.')
        node['leaf'] = True

    def complete(node):
        if node == {'leaf': True}:
            return
        require(set(node) == {0, 1, 2, 3}, 'Prefix cover is incomplete.')
        for child in node.values():
            complete(child)
    complete(trie)

    identity = ((Q(1), Q(0)), (Q(0), Q(1)))
    matrices = []
    for (a, b), e in zip(values, eps):
        matrices.append(((a, b if e else -b), (b, -a if e else a)))
    cache = {'': ((Q(0), Q(0)), identity, Q(1))}

    def word_data(word):
        if word not in cache:
            t, m, r = word_data(word[:-1])
            i = int(word[-1] == '+')
            shift = matrix_vector(m, (Q(2*i-1), Q(0)))
            cache[word] = (tuple(t[j]+shift[j] for j in (0, 1)),
                           matrix_product(m, matrices[i]), r*bounds[i])
        return cache[word]

    for u, v in leaves:
        a, _, ra = word_data(u)
        b, _, rb = word_data(v)
        d2 = sum((x-y)**2 for x, y in zip(a, b))
        radius_sum = radius*(ra+rb)
        require(d2 >= (radius_sum+gap)**2, 'Advertised gap fails on a cylinder pair.')
        require(d2-radius_sum**2 >= squared_margin, 'Advertised squared margin fails.')

    neighbourhood = report['parameter_neighbourhood']
    require(neighbourhood['orientation'] == eps, 'Neighbourhood orientation differs.')
    require(neighbourhood['metric'] == 'max(|lambda_minus-mu_minus|, |lambda_plus-mu_plus|)',
            'Unexpected parameter metric.')
    eta = rational(neighbourhood['closed_radius'])
    qbar = rational(neighbourhood['contraction_upper_bound'])
    newgap = rational(neighbourhood['contact_gap_lower_bound'])
    require(eta > 0 and q + eta <= qbar < 1, 'Neighbourhood leaves the strict contraction domain.')
    require(all(eta*eta < n for n in norms), 'Neighbourhood can contain a zero multiplier.')
    require(0 < newgap <= gap - 2*eta/(1-qbar)**2, 'Neighbourhood contact gap is not justified.')
    return {'schema': 'rational-separation-check/1', 'status': 'verified',
            'input_sha256': digest, 'checked_leaves': len(leaves),
            'maximum_word_depth': max(len(u) for u, _ in leaves),
            'rational_centre_disconnected': True, 'closed_parameter_neighbourhood_disconnected': True,
            'contact_gap_lower_bound': str(gap), 'parameter_radius': str(eta),
            'arithmetic': 'exact rational; independent real-matrix word evaluation',
            'scope': 'finite cover and stated perturbation inequalities; not formal proof-assistant verification'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate')
    parser.add_argument('--output')
    args = parser.parse_args()
    try:
        report = json.loads(Path(args.certificate).read_text())
        result = check(report)
    except (ValueError, KeyError, TypeError, ZeroDivisionError, RecursionError) as error:
        parser.error(str(error))
    output = json.dumps(result, indent=2)+'\n'
    if args.output:
        Path(args.output).write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
