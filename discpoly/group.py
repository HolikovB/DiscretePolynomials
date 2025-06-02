from typing import Any, Callable, Iterable, List, Tuple

class AbelianGroup:
    """
    Minimal interface for a (finite or finitely generated) abelian group G.
    Subclasses must implement:
      - elements() -> Iterable[Any]
      - zero() -> Any
      - add(x,y) -> Any
      - neg(x) -> Any
    """
    def elements(self) -> Iterable[Any]:
        raise NotImplementedError

    def zero(self) -> Any:
        raise NotImplementedError

    def add(self, x: Any, y: Any) -> Any:
        raise NotImplementedError

    def neg(self, x: Any) -> Any:
        raise NotImplementedError

    def __iter__(self):
        return iter(self.elements())


class CyclicGroup(AbelianGroup):
    """
    The cyclic group Z/mZ, elements are {0,1,…,m-1} with addition mod m.
    """
    def __init__(self, m: int):
        assert m >= 1
        self.m = m
        self._elements = list(range(m))

    def elements(self):
        return list(self._elements)

    def zero(self):
        return 0

    def add(self, x: int, y: int) -> int:
        return (x + y) % self.m

    def neg(self, x: int) -> int:
        return (-x) % self.m


class DirectSumGroup(AbelianGroup):
    """
    The direct sum Z/m1Z × Z/m2Z × … × Z/mnZ.  Represent each element as an n‐tuple,
    each coordinate mod the corresponding modulus.
    """
    def __init__(self, moduli: List[int]):
        assert all(mi >= 1 for mi in moduli)
        self.moduli = tuple(moduli)
        self._dims = len(moduli)

        # Precompute all elements as tuples
        def _generate_all(mods):
            if not mods:
                yield ()
            else:
                m0, *rest = mods
                for tail in _generate_all(rest):
                    for i in range(m0):
                        yield (i,) + tail

        self._elements = list(_generate_all(self.moduli))

    def elements(self):
        return list(self._elements)

    def zero(self) -> Tuple[int, ...]:
        return tuple(0 for _ in self.moduli)

    def add(self, x: Tuple[int, ...], y: Tuple[int, ...]) -> Tuple[int, ...]:
        return tuple((x[i] + y[i]) % self.moduli[i] for i in range(self._dims))

    def neg(self, x: Tuple[int, ...]) -> Tuple[int, ...]:
        return tuple((-x[i]) % self.moduli[i] for i in range(self._dims))