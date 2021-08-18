import collections
import random


EllipticCurve = collections.namedtuple('EllipticCurve', 'name p a b g n h')

ecc_curve = EllipticCurve(
    'secp256k1',
    # Field characteristic.
    p=0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f,
    # Curve coefficients.
    a=0,
    b=7,
    # Base point.
    g=(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
       0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8),
    # Subgroup order.
    n=0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141,
    # Subgroup cofactor.
    h=1,
)

# Below is Modular arithmetic ##########################################################


def inverse_modulus(k, p):
    """Function will returns the inverse of k modulo p.
    This function returns the only integer x such that (x * k) % p == 1.
    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError('division by zero')

    if k < 0:
        # k ** -1 = p - (-k) ** -1  (mod p)
        return p - inverse_modulus(-k, p)

    # Below is Extended Euclidean Algorithm
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    gcd, x, y = old_r, old_s, old_t

    assert gcd == 1
    assert (k * x) % p == 1

    return x % p


# Below is the Functions that will work on ecc_curve points #########################################


def point_on_ecc_curve(point):
    """ Function will Returns True if the provided point is on the elliptic ecc_curve."""
    if point is None:
        # None represents the point at infinity.
        return True

    x, y = point

    return (y * y - x * x * x - ecc_curve.a * x - ecc_curve.b) % ecc_curve.p == 0


def point_addition(point1, point2):
    """Function will return the result of point1 + point2 by the principle of group law."""
    assert point_on_ecc_curve(point1)
    assert point_on_ecc_curve(point2)

    if point1 is None:
        # 0 + point2 = point2
        return point2
    if point2 is None:
        # point1 + 0 = point1
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        # point1 + (-point1) = 0
        return None

    if x1 == x2:
        # This is the case point1 == point2.
        m = (3 * x1 * x1 + ecc_curve.a) * inverse_modulus(2 * y1, ecc_curve.p)
    else:
        # This is the case point1 != point2.
        m = (y1 - y2) * inverse_modulus(x1 - x2, ecc_curve.p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % ecc_curve.p, -y3 % ecc_curve.p)

    assert point_on_ecc_curve(result)

    return result


def point_negative(point):
    pass


def point_multiplication(k, point):
    """Function will Return k * point calculated using the double and point_addition algorithm."""
    assert point_on_ecc_curve(point)

    if k % ecc_curve.n == 0 or point is None:
        return None

    if k < 0:
        # k * point = -k * (-point)
        return point_multiplication(-k, point_negative(point))

    result = None
    addend = point

    while k:
        if k & 1:
            # Add.
            result = point_addition(result, addend)

        # Double.
        addend = point_addition(addend, addend)

        k >>= 1

    assert point_on_ecc_curve(result)

    return result


# Keypair generation Public key and private key generation ################################################


def make_keypair_pr_pb():
    """Generating a random private key and public key pair for use"""
    private_key = random.randrange(1, ecc_curve.n)
    public_key = point_multiplication(private_key, ecc_curve.g)

    return private_key, public_key


print("Generator/Basepoint:\t", ecc_curve.g)

alice_secret_key, alice_public_key = make_keypair_pr_pb()
bob_secret_key, bob_public_key = make_keypair_pr_pb()

print("Alice\'s secret key:\t", alice_secret_key)
print("Alice\'s public key:\t", alice_public_key)
print("Bob\'s secret key:\t", bob_secret_key)
print("Bob\'s public key:\t", bob_public_key)

print("==========================")

shared_secret1 = point_multiplication(bob_secret_key, alice_public_key)
shared_secret2 = point_multiplication(alice_secret_key, bob_public_key)

print("==========================")
print("Alice\'s shared key:\t", shared_secret1)
print("Bob\'s shared key:\t", shared_secret2)

print("==========================")
print("The shared secret value/key is the x-value: \t", (shared_secret1[0]))