from fractions import Fraction
from typing import Union

Real = Union[float, Fraction]

def mod1_float(x: float) -> float:
    """
    Reduce a float x into [0,1) by subtracting floor(x). 
    (Equivalent to x % 1 for non‐extreme floats.)
    """
    return x - float(int(x))

def mod1_frac(x: Fraction) -> Fraction:
    """
    Reduce a Fraction x into [0,1) by subtracting its integer part.
    """
    return x - (x.numerator // x.denominator)

def add_mod1(a: Real, b: Real) -> Real:
    """
    Add two elements in R/Z, returning the result mod 1.
    If both are Fraction, do exact rational mod 1; otherwise coerce to float.
    """
    if isinstance(a, Fraction) and isinstance(b, Fraction):
        return mod1_frac(a + b)
    else:
        return mod1_float(float(a) + float(b))

def sub_mod1(a: Real, b: Real) -> Real:
    """
    Compute (a - b) mod 1 in R/Z.
    """
    if isinstance(a, Fraction) and isinstance(b, Fraction):
        return mod1_frac(a - b)
    else:
        return mod1_float(float(a) - float(b))

def zero_mod1(typ=float) -> Real:
    """
    Return the “zero” element in R/Z, either 0.0 (float) or Fraction(0,1).
    """
    if typ is Fraction:
        return Fraction(0, 1)
    return 0.0