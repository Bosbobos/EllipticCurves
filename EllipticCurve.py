import math
import random


class EllipticCurve:
    def __init__(self, p, a, b):
        self.p = p
        self.a = a
        self.b = b
        if (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p == 0:
            raise ValueError("Кривая не удовлетворяет условию 4a³ + 27b² ≠ 0 mod p.")

    def is_on_curve(self, point):
        if point == 'inf':
            return True
        x, y = point
        return (pow(y, 2, self.p) - (pow(x, 3, self.p) + self.a * x + self.b) % self.p) % self.p == 0


class ECPoint:
    def __init__(self, curve, x, y):
        self.curve = curve
        self.x = x % curve.p
        self.y = y % curve.p
        if not curve.is_on_curve((self.x, self.y)):
            raise ValueError("Точка не принадлежит кривой.")

    def __eq__(self, other):
        if other == 'inf':
            return False
        return self.x == other.x and self.y == other.y and self.curve == other.curve

    def __add__(self, other):
        if other == 'inf':
            return self
        if self == 'inf':
            return other
        if self.x == other.x and (self.y + other.y) % self.curve.p == 0:
            return 'inf'
        if self != other:
            m = ((other.y - self.y) * mod_inverse(other.x - self.x, self.curve.p)) % self.curve.p
        else:
            m = ((3 * pow(self.x, 2, self.curve.p) + self.curve.a) * mod_inverse(2 * self.y,
                                                                                 self.curve.p)) % self.curve.p
        x3 = (pow(m, 2, self.curve.p) - self.x - other.x) % self.curve.p
        y3 = (m * (self.x - x3) - self.y) % self.curve.p
        return ECPoint(self.curve, x3, y3)

    def __radd__(self, other):
        return self + other

    def __mul__(self, scalar):
        result = 'inf'
        current = self
        while scalar > 0:
            if scalar % 2 == 1:
                result = result + current
            current = current + current
            scalar = scalar // 2
        return result

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __repr__(self):
        return f"({self.x}, {self.y})"


def mod_inverse(a, p):
    g, x, y = extended_gcd(a, p)
    if g != 1:
        raise ValueError("Обратный элемент не существует.")
    return x % p


def extended_gcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = extended_gcd(b % a, a)
        return (g, x - (b // a) * y, y)


def legendre_symbol(a, p):
    ls = pow(a, (p - 1) // 2, p)
    if ls == p - 1:
        return -1
    return ls


def tonelli_shanks(n, p):
    if legendre_symbol(n, p) != 1:
        return None
    if n == 0:
        return 0
    if p == 2:
        return p
    if p % 4 == 3:
        x = pow(n, (p + 1) // 4, p)
        return x
    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1
    z = 2
    while legendre_symbol(z, p) != -1:
        z += 1
    c = pow(z, Q, p)
    x = pow(n, (Q + 1) // 2, p)
    t = pow(n, Q, p)
    m = S
    while t != 1:
        i, temp = 0, t
        while temp != 1 and i < m:
            temp = pow(temp, 2, p)
            i += 1
        if i == m:
            return None
        b = pow(c, 1 << (m - i - 1), p)
        x = (x * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        m = i
    return x


def find_points(curve):
    points = ['inf']
    for x in range(curve.p):
        rhs = (pow(x, 3, curve.p) + curve.a * x + curve.b) % curve.p
        if legendre_symbol(rhs, curve.p) == 1:
            y = tonelli_shanks(rhs, curve.p)
            points.append(ECPoint(curve, x, y))
            points.append(ECPoint(curve, x, (-y) % curve.p))
        elif rhs == 0:
            points.append(ECPoint(curve, x, 0))
    return points


def naive_order(curve):
    return len(find_points(curve))


def bsgs(curve, P, Q):
    m = int(math.ceil(math.sqrt(curve.p)))
    table = {}
    for j in range(m):
        point = j * P
        table[(point.x, point.y)] = j
    m_point = m * P
    current = Q
    for k in range(m):
        if (current.x, current.y) in table:
            j = table[(current.x, current.y)]
            return (m * k - j) % curve.p
        current = current + m_point
    return None


def curve_order(curve, max_trials=10):
    points = find_points(curve)
    #if len(points) < 2 ** 10 or len(points) > 2 ** 512:
    #    raise ValueError("Порядок кривой не входит в требуемый диапазон.")
    return len(points)


def find_prime_subgroups(curve):
    order = curve_order(curve)
    factors = {}
    n = order
    i = 2
    while i * i <= n:
        while n % i == 0:
            factors[i] = factors.get(i, 0) + 1
            n = n // i
        i += 1
    if n > 1:
        factors[n] = 1
    prime_factors = list(factors.keys())
    subgroups = []
    for p in prime_factors:
        for point in find_points(curve):
            if point != 'inf':
                candidate = (order // p) * point
                if candidate == 'inf':
                    subgroups.append(p)
                    break
    return subgroups


# Пример использования
p = 11
a = 1
b = 6
curve = EllipticCurve(p, a, b)
points = find_points(curve)
print("Точки кривой:", points)
print("Порядок кривой (наивный метод):", naive_order(curve))
print("Порядок кривой (BSGS):", curve_order(curve))
P = ECPoint(curve, 2, 7)
print("2P:", P + P)
print("3P:", 3 * P)
print("Подгруппы простого порядка:", find_prime_subgroups(curve))
