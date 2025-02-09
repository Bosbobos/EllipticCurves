from ECPointInf import ECPointInf

class ECPoint:
    def __init__(self, curve, x, y):
        self.curve = curve
        self.x = x % curve.p
        self.y = y % curve.p
        if not curve.is_on_curve((self.x, self.y)):
            raise ValueError("Точка не принадлежит кривой.")

    def __add__(self, other):
        if isinstance(other, ECPointInf):
            return self
        # Случаи сложения точек
        if self.x == other.x:
            if (self.y + other.y) % self.curve.p == 0:
                return ECPointInf(self.curve)
            else:
                # Удвоение точки
                return self.double()
        # Формулы сложения
        m = ((other.y - self.y) * mod_inverse(other.x - self.x, self.curve.p)) % self.curve.p
        x3 = (m**2 - self.x - other.x) % self.curve.p
        y3 = (m * (self.x - x3) - self.y) % self.curve.p
        return ECPoint(self.curve, x3, y3)

    def double(self):
        if (self.y % self.curve.p == 0):
            return ECPointInf(self.curve)

        m = ((3 * self.x**2 + self.curve.a) * mod_inverse(2 * self.y, self.curve.p)) % self.curve.p
        x3 = (m**2 - 2 * self.x) % self.curve.p
        y3 = (m * (self.x - x3) - self.y) % self.curve.p
        return ECPoint(self.curve, x3, y3)

    def __radd__(self, other):
        return self + other

    def __mul__(self, scalar):
        if scalar == 0:
            return ECPointInf(self.curve)
        result = ECPointInf(self.curve)
        current = self
        while scalar > 0:
            if scalar % 2 == 1:
                result = result + current
            current = current.double()
            if isinstance(current, ECPointInf):
                break
            scalar = scalar // 2
        return result

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __repr__(self):
        return f"({self.x}, {self.y})"

    def __neg__(self):
        return point_neg(self)


def mod_inverse(a, p):
    a %= p
    g, x, y = extended_gcd(a, p)
    if g != 1:
        raise ValueError("Обратный элемент не существует.")
    return x % p

def extended_gcd(a, b):
    """
    Вычисляет расширенный алгоритм Евклида для нахождения наибольшего общего делителя
    двух чисел a и b, а также коэффициентов x и y, таких что a * x + b * y = gcd(a, b).

    Аргументы:
    a -- Первое целое число
    b -- Второе целое число

    Возвращает:
    Кортеж (g, x, y), где:
    g -- Наибольший общий делитель a и b,
    x -- Коэффициент для a в линейном представлении НОД,
    y -- Коэффициент для b в линейном представлении НОД.
    """
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = extended_gcd(b % a, a)
        return (g, x - (b // a) * y, y)

def point_neg(P):
    """
    Возвращает отрицание точки P.
    Для точки в бесконечности возвращается она же.
    """
    if isinstance(P, ECPointInf):
        return P
    return ECPoint(P.curve, P.x, (-P.y) % P.curve.p)