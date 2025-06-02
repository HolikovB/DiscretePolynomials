from typing import Any, Dict, List, Union
from fractions import Fraction
import itertools

from discpoly.group import AbelianGroup
from discpoly.utils import add_mod1, sub_mod1, mod1_frac, mod1_float, zero_mod1

Real = Union[float, Fraction]

class PolynomialFunction:
    """
    Represents a function f: G -> R/Z, where G is an AbelianGroup and
    values are stored in mod‐1 arithmetic (float or Fraction).
    Internally: self.values is a dict mapping each g∈G.elements() to a Real in [0,1).
    """
    def __init__(
        self,
        G: AbelianGroup,
        values: Dict[Any, Real],
        use_fraction: bool = False
    ):
        self.G = G
        # Ensure every element of G has a value provided
        missing = [g for g in G if g not in values]
        if missing:
            raise ValueError(f"Missing values for elements: {missing!r}")

        # Coerce values into float or Fraction, reduce mod1
        if use_fraction:
            self.values: Dict[Any, Fraction] = {}
            for g, v in values.items():
                fv = Fraction(v) if not isinstance(v, Fraction) else v
                # reduce fv to [0,1)
                fv_mod = fv - (fv.numerator // fv.denominator)
                self.values[g] = fv_mod
            self._zero = zero_mod1(Fraction)
        else:
            self.values: Dict[Any, float] = {}
            for g, v in values.items():
                fv = float(v)
                self.values[g] = fv - float(int(fv))  # mod1
            self._zero = zero_mod1(float)

    def __call__(self, g: Any) -> Real:
        return self.values[g]

    def __add__(self, other: "PolynomialFunction") -> "PolynomialFunction":
        """
        Pointwise addition modulo 1: (f + g)(x) = f(x) + g(x) (mod 1).
        Both must live on the same group G, and both must agree on use_fraction.
        """
        if self.G is not other.G:
            raise ValueError("Cannot add PolynomialFunctions on different groups.")
        # Check that they use the same numeric type
        use_frac_self = isinstance(self._zero, Fraction)
        use_frac_other = isinstance(other._zero, Fraction)
        if use_frac_self != use_frac_other:
            raise ValueError("Cannot add: one PolynomialFunction uses Fraction, the other float.")
        use_fraction = use_frac_self

        new_vals: Dict[Any, Real] = {}
        for x in self.G:
            a = self.values[x]
            b = other.values[x]
            # safe to call add_mod1 for both Fraction&Fraction or float&float
            new_vals[x] = add_mod1(a, b)
        return PolynomialFunction(self.G, new_vals, use_fraction=use_fraction)
    
    def __mul__(self, other: "PolynomialFunction") -> "PolynomialFunction":
        """
        Pointwise multiplication modulo 1: (f * g)(x) = [f(x) * g(x)] (mod 1).
        Both must live on the same group G, and both must agree on use_fraction.
        """
        if self.G is not other.G:
            raise ValueError("Cannot multiply PolynomialFunctions on different groups.")
        use_frac_self = isinstance(self._zero, Fraction)
        use_frac_other = isinstance(other._zero, Fraction)
        if use_frac_self != use_frac_other:
            raise ValueError("Cannot multiply: one PolynomialFunction uses Fraction, the other float.")
        use_fraction = use_frac_self

        new_vals: Dict[Any, Real] = {}
        if use_fraction:
            # exact rational multiplication
            for x in self.G:
                a: Fraction = self.values[x]
                b: Fraction = other.values[x]
                prod = a * b
                new_vals[x] = mod1_frac(prod)
        else:
            # float multiplication
            for x in self.G:
                a: float = self.values[x]  # already in [0,1)
                b: float = other.values[x]
                prod = a * b
                new_vals[x] = mod1_float(prod)
        return PolynomialFunction(self.G, new_vals, use_fraction=use_fraction)
    
    def difference(self, h: Any) -> "PolynomialFunction":
        """
        Compute Δ_h f:  x ↦ f(x + h) − f(x)  (mod 1).
        Returns a new PolynomialFunction on the same group.
        """
        new_vals: Dict[Any, Real] = {}
        use_frac = isinstance(self._zero, Fraction)
        for x in self.G:
            xh = self.G.add(x, h)
            a = self.values[xh]
            b = self.values[x]
            new_vals[x] = sub_mod1(a, b)
        return PolynomialFunction(self.G, new_vals, use_fraction=use_frac)

    def iterated_difference(self, hs: List[Any]) -> "PolynomialFunction":
        """
        Compute Δ_{h1, …, hk} f by applying .difference(h1), then .difference(h2), …
        """
        f_curr = self
        for h in hs:
            f_curr = f_curr.difference(h)
        return f_curr

    def is_zero_function(self, tol: float = 1e-9) -> bool:
        """
        Check if f(x) ≡ 0 in R/Z for all x∈G.
        If float: allow |val|<tol or |val−1|<tol. If Fraction: exact eq = 0.
        """
        if isinstance(self._zero, Fraction):
            return all(val == 0 for val in self.values.values())
        else:
            for val in self.values.values():
                if not (abs(val) < tol or abs(val - 1.0) < tol):
                    return False
            return True

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
