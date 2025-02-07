from ECPointInf import ECPointInf


class EllipticCurve:
    def __init__(self, p, a, b):
        self.p = p
        self.a = a % p
        self.b = b % p
        # Проверка условия 4a³ + 27b² ≠ 0 mod p
        if (4 * pow(self.a, 3, p) + 27 * pow(self.b, 2, p)) % p == 0:
            raise ValueError("Кривая не удовлетворяет условию 4a³ + 27b² ≠ 0 mod p.")

    def is_on_curve(self, point):
        if isinstance(point, ECPointInf):
            return True
        x, y = point
        return (pow(y, 2, self.p) - (pow(x, 3, self.p) + self.a * x + self.b) % self.p) % self.p == 0
