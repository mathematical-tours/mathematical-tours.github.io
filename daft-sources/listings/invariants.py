"""Exact Reynolds averages for small finite matrix groups (characteristic zero)."""
import sympy as sp


def matrix_group(generators, limit=256):
    n = generators[0].rows
    group = [sp.eye(n)]
    known = {tuple(group[0])}
    for element in group:
        for generator in generators:
            candidate = (element*generator).applyfunc(sp.simplify)
            key = tuple(candidate)
            if key not in known:
                if len(group) >= limit:
                    raise ValueError('The group did not close within the explicit size limit.')
                known.add(key)
                group.append(candidate)
    return group


def exponents_of_degree(variables, degree):
    if variables == 1:
        yield (degree,)
    else:
        for first in range(degree+1):
            for rest in exponents_of_degree(variables-1, degree-first):
                yield (first,)+rest


def reynolds_generators(group, variables):
    """An overcomplete generating set from Noether's degree bound."""
    forms = [A*sp.Matrix(variables) for A in group]  # Inverses give the same average.
    found = {}
    for degree in range(1, len(group)+1):
        for powers in exponents_of_degree(len(variables), degree):
            average = sum(sp.prod(u**j for u, j in zip(row, powers)) for row in forms)/len(group)
            polynomial = sp.Poly(sp.expand(average), *variables, extension=True)
            if not polynomial.is_zero:
                normalized = polynomial.monic().as_expr()
                found.setdefault(normalized, None)
    return list(found)
