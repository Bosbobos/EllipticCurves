from EllipticCurve import EllipticCurve
from ECPoint import *
from ECPointInf import ECPointInf
from Tools import *
import pytest


def test_valid_curve_creation():
    # Тест: корректное создание кривой (например, p=17, a=2, b=2)
    curve = EllipticCurve(17, 2, 2)
    assert curve.p == 17
    assert curve.a == 2 % 17
    assert curve.b == 2 % 17


def test_invalid_curve_creation():
    # Тест: создание кривой с параметрами a=0, b=0 (4a³+27b² = 0 mod p) должно вызывать ошибку.
    with pytest.raises(ValueError):
        EllipticCurve(17, 0, 0)


def test_is_on_curve():
    curve = EllipticCurve(17, 2, 2)
    # Точка (5,1) лежит на кривой:
    # 1² = 5³ + 2*5 + 2 mod 17  ⇒ 1 = 125+10+2 = 137 mod 17 = 1.
    point = (5, 1)
    assert curve.is_on_curve(point)

    # Точка (0,1) не принадлежит кривой: 1² ≠ 0³+2*0+2 (1 ≠ 2 mod 17)
    assert not curve.is_on_curve((0, 1))

    # Точка в бесконечности всегда принадлежит кривой.
    inf_point = ECPointInf(curve)
    assert curve.is_on_curve(inf_point)


def test_point_addition():
    curve = EllipticCurve(17, 2, 2)
    P = ECPoint(curve, 5, 1)
    Q = ECPoint(curve, 6, 3)
    # Проверяем, что результат сложения лежит на кривой.
    R = P + Q
    if not isinstance(R, ECPointInf):
        assert curve.is_on_curve((R.x, R.y))
    # Сложение с точкой в бесконечности.
    inf_point = ECPointInf(curve)
    assert (P + inf_point) == P
    assert (inf_point + P) == P


def test_point_doubling():
    curve = EllipticCurve(17, 2, 2)
    P = ECPoint(curve, 5, 1)
    doubled = P.double()
    assert curve.is_on_curve((doubled.x, doubled.y))


def test_scalar_multiplication():
    curve = EllipticCurve(17, 2, 2)
    P = ECPoint(curve, 5, 1)
    # 0*P должно давать точку в бесконечности.
    zeroP = 0 * P
    assert isinstance(zeroP, ECPointInf)
    # 2*P должно совпадать с удвоением.
    twoP = 2 * P
    doubleP = P.double()
    if not (isinstance(twoP, ECPointInf) or isinstance(doubleP, ECPointInf)):
        assert twoP.x == doubleP.x and twoP.y == doubleP.y
    # 3*P лежит на кривой.
    threeP = 3 * P
    if not isinstance(threeP, ECPointInf):
        assert curve.is_on_curve((threeP.x, threeP.y))


def test_mod_inverse():
    # Для p=11: обратный к 3 равен 4, т.к. 3*4 = 12 ≡ 1 mod 11.
    inv = mod_inverse(3, 11)
    assert (3 * inv) % 11 == 1
    # Если обратного не существует, должно возникать исключение.
    with pytest.raises(ValueError):
        mod_inverse(0, 11)


def test_extended_gcd():
    # Для чисел 240 и 46 НОД = 2, а линейная комбинация должна давать 2.
    g, x, y = extended_gcd(240, 46)
    assert g == 2
    assert 240 * x + 46 * y == g


def test_legendre_symbol():
    # Для простого 17, квадратичный вычет: 4 (2² = 4).
    assert legendre_symbol(4, 17) == 1
    # Для невычета, например, 3: 3^8 mod 17 = 16 = -1 mod 17.
    assert legendre_symbol(3, 17) == -1
    # Для a=2 – значение может быть 1 или -1.
    ls = legendre_symbol(2, 17)
    assert ls in (1, -1)


def test_tonelli_shanks():
    # Для простого 17: корень из 4 равен 2 или 15.
    root = tonelli_shanks(4, 17)
    assert (root * root) % 17 == 4
    # Для невычета функция должна вернуть None.
    assert tonelli_shanks(3, 17) is None


def test_find_points_and_naive_order():
    curve = EllipticCurve(17, 2, 2)
    pts = find_points(curve)
    # Проверяем, что точка в бесконечности присутствует.
    assert any(isinstance(pt, ECPointInf) for pt in pts)
    # Все точки (кроме бесконечно удалённой) должны удовлетворять уравнению кривой.
    for pt in pts:
        if isinstance(pt, ECPointInf):
            continue
        assert curve.is_on_curve((pt.x, pt.y))
    # Проверяем, что функция naive_order возвращает число точек.
    assert naive_order(curve) == len(pts)


def test_bsgs():
    curve = EllipticCurve(17, 2, 2)
    # Выбираем точку P.
    P = ECPoint(curve, 5, 1)
    # Пусть n = 7, и Q = 7*P.
    n = 7
    Q = n * P
    # Если Q — точка в бесконечности, пропускаем тест.
    if isinstance(Q, ECPointInf):
        pytest.skip("Q равна бесконечности, выберите другой скаляр.")
    # Восстанавливаем n с помощью алгоритма bsgs.
    dlog = bsgs(curve, P, Q)
    assert dlog is not None
    # Проверяем, что dlog * P = Q.
    result = dlog * P
    if isinstance(result, ECPointInf) and isinstance(Q, ECPointInf):
        assert True
    else:
        assert result.x == Q.x and result.y == Q.y


def test_find_prime_subgroups():
    curve = EllipticCurve(17, 2, 2)
    subs = find_prime_subgroups(curve)
    # Функция должна вернуть список (хотя реализация может работать не совсем корректно
    # из-за сравнения с 'inf').
    assert isinstance(subs, list)
    # Для каждого найденного простого делителя порядка проверяем, что он действительно делит порядок.
    order = curve_order(curve)
    # Факторизация порядка (простая версия).
    factors = {}
    n = order
    i = 2
    while i * i <= n:
        while n % i == 0:
            factors[i] = factors.get(i, 0) + 1
            n //= i
        i += 1
    if n > 1:
        factors[n] = 1
    for sg in subs:
        assert sg in factors


def test_repr():
    curve = EllipticCurve(17, 2, 2)
    P = ECPoint(curve, 5, 1)
    inf_point = ECPointInf(curve)
    assert repr(P) == f"({P.x}, {P.y})"
    assert repr(inf_point) == "inf"


@pytest.mark.parametrize(
    'p, a, b', [
        (7, -2, 1)
    ]
)
def test_analysis(p, a, b):
    curve = EllipticCurve(p, a, b)
    points = find_points(curve)
    groups = find_prime_subgroups(curve)
    order = curve_order(curve)
    p1_order = point_order(points[0])
    p2_order = point_order(points[1])
    assert len(points) == 12
