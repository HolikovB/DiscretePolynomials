import pytest
from fractions import Fraction

from discpoly.group import DirectSumGroup, CyclicGroup
from discpoly.discpoly import PolynomialFunction

def test_constant_function():
    G3 = CyclicGroup(3)
    vals = {x: Fraction(2, 3) for x in G3.elements()}
    f = PolynomialFunction(G3, vals, use_fraction=True)
    assert f.degree_leq(0)
    assert not f.degree_leq(-1)

def test_linear_function_on_z4():
    G4 = CyclicGroup(4)
    # f(x) = x/4 (mod 1)
    vals = {x: Fraction(x, 4) for x in G4.elements()}
    f_lin = PolynomialFunction(G4, vals, use_fraction=True)
    assert f_lin.degree_leq(1)
    assert not f_lin.degree_leq(0)

def test_quadratic_on_z4():
    G4 = CyclicGroup(4)
    # f(x) = x^2 / 8  (mod 1).  Then degree = 2, not ≤1.
    vals = {}
    for x in G4.elements():
        val = Fraction(x*x, 8)
        vals[x] = val - (val.numerator // val.denominator)  # reduce mod 1
    f_quad = PolynomialFunction(G4, vals, use_fraction=True)
    assert f_quad.degree_leq(2)
    assert not f_quad.degree_leq(1)

def test_multivariate_degree_on_z4x4():
    G = DirectSumGroup([4,4])
    # B_{1,1}(x1,x2) = x1*x2  has total degree = 2
    vals = {}
    for (x1, x2) in G.elements():
        vals[(x1, x2)] = Fraction(x1 * x2, 2)  # integer
    f_b11 = PolynomialFunction(G, vals, use_fraction=True)
    assert f_b11.degree_leq(2)
    assert not f_b11.degree_leq(1)