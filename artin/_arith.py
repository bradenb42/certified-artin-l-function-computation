"""Small arithmetic helpers shared across the computational package."""

from math import gcd


def multiplicative_order(a: int, modulus: int) -> int:
    """Return the multiplicative order of ``a`` modulo ``modulus``."""
    if modulus == 1:
        return 1
    if gcd(a, modulus) != 1:
        raise ValueError("multiplicative order requires coprime arguments")

    order, value = 1, a % modulus
    while value != 1:
        value = value * a % modulus
        order += 1
    return order
