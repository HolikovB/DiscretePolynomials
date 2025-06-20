from typing import Any, Dict, List, Union
from fractions import Fraction
import itertools

from discpoly.group import AbelianGroup
from discpoly.utils import add_mod1, sub_mod1, mod1_frac, zero_mod1



class PolynomialFunction:
    """
    Represents a function f: G -> R/Z, where G is an AbelianGroup and
    values are stored in mod‐1 arithmetic.
    Internally: self.values is a dict mapping each g∈G.elements() to a Real in [0,1).
    """
    def __init__(
        self,
        G: AbelianGroup,
        values: Dict[Any, Fraction],
    ):
        self.G = G
        # Ensure every element of G has a value provided
        missing = [g for g in G if g not in values]
        if missing:
            raise ValueError(f"Missing values for elements: {missing!r}")
        
        self.values: Dict[Any, Fraction] = {}
        for g, v in values.items():
            fv = Fraction(v) if not isinstance(v, Fraction) else v
            # reduce fv to [0,1)
            fv_mod = fv - (fv.numerator // fv.denominator)
            self.values[g] = fv_mod
        self._zero = zero_mod1(Fraction)


    def __call__(self, g: Any) -> Fraction:
        return self.values[g]

    def __add__(self, other: "PolynomialFunction") -> "PolynomialFunction":
        """
        Pointwise addition modulo 1: (f + g)(x) = f(x) + g(x) (mod 1).
        Both must live on the same group G, and both must agree on use_fraction.
        """
        if self.G is not other.G:
            raise ValueError("Cannot add PolynomialFunctions on different groups.")
        # Check that they use the same numeric type


        new_vals: Dict[Any, Fraction] = {}
        for x in self.G:
            a = self.values[x]
            b = other.values[x]
            new_vals[x] = add_mod1(a, b)
        return PolynomialFunction(self.G, new_vals)
    
    def __mul__(self, other: Union["PolynomialFunction","Fraction"]) -> "PolynomialFunction":
        """
        Pointwise multiplication modulo 1: (f * g)(x) = [f(x) * g(x)] (mod 1).
        Both must live on the same group G, and both must agree on use_fraction.
        """
        if isinstance(other, PolynomialFunction):
            if self.G is not other.G:
                raise ValueError("Cannot multiply PolynomialFunctions on different groups.")

            new_vals: Dict[Any, Fraction] = {}

            # exact rational multiplication
            for x in self.G:
                a: Fraction = self.values[x]
                b: Fraction = other.values[x]
                prod = a * b
                new_vals[x] = mod1_frac(prod)

            return PolynomialFunction(self.G, new_vals)
        elif isinstance(other, Fraction):
            use_frac_self = isinstance(self._zero, Fraction)
            use_frac_scalar = isinstance(other, Fraction)
            if use_frac_self != use_frac_scalar:
                raise ValueError(
                    "Cannot multiply: scalar type must match function’s numeric type "
                    "(both float or both Fraction)."
                )
            use_fraction = use_frac_self

            new_vals: Dict[Any, Fraction] = {}
            if use_fraction:
                # exact Fraction * Fraction mod 1
                scalar: Fraction = other
                for x in self.G:
                    a: Fraction = self.values[x]
                    prod = scalar * a
                    new_vals[x] = mod1_frac(prod)
            return PolynomialFunction(self.G, new_vals)

        else:
            raise ValueError("Operand for multiplication must be either a PolynomialFunction or a numeric scalar.")
    
    def difference(self, h: Any) -> "PolynomialFunction":
        """
        Compute Δ_h f:  x ↦ f(x + h) − f(x)  (mod 1).
        Returns a new PolynomialFunction on the same group.
        """
        new_vals: Dict[Any, Fraction] = {}
        for x in self.G:
            xh = self.G.add(x, h)
            a = self.values[xh]
            b = self.values[x]
            new_vals[x] = sub_mod1(a, b)
        return PolynomialFunction(self.G, new_vals)


    def iterated_difference(self, hs: List[Any]) -> "PolynomialFunction":
        """
        Compute Δ_{h1, …, hk} f by applying .difference(h1), then .difference(h2), …
        """
        f_curr = self
        for h in hs:
            f_curr = f_curr.difference(h)
        return f_curr

    def is_zero_function(self) -> bool:
        """
        Check if f(x) ≡ 0 in R/Z for all x∈G.
        """
        return all(val == 0 for val in self.values.values())


    def degree_leq(self, d: int) -> bool:
        """
        If it passes this test, check that for every (d+1)-tuple (h1,…,h_{d+1})∈G^{d+1},
        Δ_{h1,…,h_{d+1}} f is the zero function.  Return True iff so.
        Warning: runtime is |G|^{(d+1)}, so only small d/|G| feasible.
        """
        if d < 0:
            return False
        elems = list(self.G.elements())
        for hs in itertools.product(elems, repeat=d+1):
            if not self.iterated_difference(list(hs)).is_zero_function():
                return False
        return True

    def explain(self) -> None:
        """
        Print a table of f(g) for each g ∈ G.
        """
        print("PolynomialFunction values (mod 1):")
        for g in self.G:
            print(f"  f({g}) = {self.values[g]}")
