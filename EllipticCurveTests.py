import pytest
from EllipticCurve import EllipticCurve, ECPoint, mod_inverse, find_points, curve_order, find_prime_subgroups


# Фикстуры для тестовых данных
@pytest.fixture
def small_curve():
    return EllipticCurve(p=11, a=1, b=6)


@pytest.fixture
def small_curve_point(small_curve):
    return ECPoint(small_curve, x=2, y=7)


@pytest.fixture
def large_curve():
    return EllipticCurve(p=23, a=1, b=1)


@pytest.fixture
def large_curve_point(large_curve):
    return ECPoint(large_curve, x=3, y=10)


# Тесты валидации кривой
def test_valid_curve_creation():
    with pytest.raises(ValueError):
        EllipticCurve(p=11, a=4, b=4)  # 4*4³ + 27*4² = 256 + 432 = 688 ≡ 688%11=6 ≠ 0
    EllipticCurve(p=11, a=1, b=6)  # Должно пройти без ошибок


# Тесты операций с точками
def test_point_addition(small_curve, small_curve_point):
    P = small_curve_point
    inf = ECPoint(small_curve, None, None)  # Бесконечно удаленная точка

    # P + inf = P
    assert P + inf == P

    # P + P
    P2 = P + P
    assert P2.x == 5 and P2.y == 2

    # Обратная точка
    Q = ECPoint(small_curve, 2, 4)
    assert P + Q == inf


def test_scalar_multiplication(small_curve_point):
    P = small_curve_point
    assert 2 * P == P + P
    assert 3 * P == P + P + P
    assert 10 * P + P == ECPoint(small_curve_point.curve, None, None)  # Порядок точки


# Тесты арифметики поля
def test_mod_inverse():
    assert mod_inverse(3, 11) == 4  # 3*4=12 ≡ 1 mod 11
    assert mod_inverse(7, 31) == 9  # 7*9=63 ≡ 1 mod 31


# Тесты подсчета точек
def test_naive_order(small_curve):
    assert len(find_points(small_curve)) == 12  # Проверьте реальное значение для вашей кривой


def test_curve_order(large_curve):
    assert curve_order(large_curve) == 28  # Для кривой y² = x³ + x + 1 над F_23


# Тесты подгрупп
def test_prime_subgroups(small_curve):
    subgroups = find_prime_subgroups(small_curve)
    assert set(subgroups) == {2, 3}  # Для порядка 12 = 2²*3


# Тесты исключений
def test_invalid_point(small_curve):
    with pytest.raises(ValueError):
        ECPoint(small_curve, 1, 1)  # Не на кривой: 1² ≢ 1³ + 1*1 + 6 mod 11


def test_zero_multiplication(small_curve_point):
    assert 0 * small_curve_point == ECPoint(small_curve_point.curve, None, None)


# Параметризованный тест для операций
@pytest.mark.parametrize("k,expected", [
    (0, "inf"),
    (1, (2, 7)),
    (2, (5, 2)),
    (3, (8, 8)),
])
def test_multiplication_results(k, expected, small_curve_point):
    result = k * small_curve_point
    if expected == "inf":
        assert result.x is None and result.y is None
    else:
        assert (result.x, result.y) == expected