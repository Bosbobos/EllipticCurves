from EllipticCurve import EllipticCurve
from Tools import *


if __name__ == '__main__':
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
