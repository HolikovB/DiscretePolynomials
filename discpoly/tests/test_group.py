import pytest
from discpoly.group import CyclicGroup, DirectSumGroup

def test_cyclic_group():
    G5 = CyclicGroup(5)
    assert G5.zero() == 0
    assert G5.add(2, 3) == 0  # mod 5
    assert G5.neg(3) == 2

def test_direct_sum_group_4x4():
    G44 = DirectSumGroup([4,4])
    zero = G44.zero()
    assert zero == (0, 0)

    a = (1, 2)
    b = (3, 3)
    s = G44.add(a, b)
    # (1+3)%4=0, (2+3)%4=1
    assert s == (0, 1)

    # check inverses
    assert G44.neg((1,2)) == ((-1) % 4, (-2) % 4) == (3, 2)
    assert G44.neg((0,0)) == (0, 0)