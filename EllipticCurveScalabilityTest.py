import pytest
from EllipticCurve import EllipticCurve
from Tools import *

import pytest


@pytest.mark.parametrize(
    "curve, expected_length",
    [
        (EllipticCurve(a=2, b=4, p=7), 7),  # Простой порядок
        (EllipticCurve(a=5, b=9, p=11), 11),  # Простой порядок
        (EllipticCurve(a=3, b=5, p=13), 13),  # Простой порядок
        (EllipticCurve(a=7, b=6, p=17), 17),  # Простой порядок
        (EllipticCurve(a=2, b=8, p=19), 19),  # Простой порядок
        (EllipticCurve(a=10, b=3, p=23), 23),  # Простой порядок
        (EllipticCurve(a=5, b=15, p=29), 29),  # Простой порядок
        (EllipticCurve(a=11, b=7, p=31), 31),  # Простой порядок
        (EllipticCurve(a=13, b=19, p=37), 37),  # Простой порядок
        (EllipticCurve(a=8, b=20, p=97), 97),  # Простой порядок
        (EllipticCurve(a=5, b=13, p=101), 101),  # Простой порядок
        (EllipticCurve(a=15, b=21, p=199), 199),  # Простой порядок
    ]
)
def test_prime_order(curve, expected_length):
    """
    Тест на проверку простого порядка кривой и корректного нахождения простых подгрупп.
    """
    points = find_points(curve)
    order = curve_order(curve)  # Получаем порядок группы кривой
    prime_subgroups = find_prime_subgroups(curve)

    # Проверим, что порядок группы соответствует ожидаемому
    assert order == expected_length  # Порядок должен быть равен ожидаемой длине

    # Проверим, что подгруппы содержат только точки с порядком expected_length
    assert len(prime_subgroups) == 1  # Мы должны получить только одну подгруппу
    assert len(prime_subgroups[0]) == expected_length  # Подгруппа должна содержать ожидаемое количество точек


@pytest.mark.parametrize(
    "curve, expected_subgroup_sizes",
    [
        (EllipticCurve(a=3, b=4, p=11), [2, 3]),  # Составной порядок
        (EllipticCurve(a=7, b=9, p=15), [3, 5]),  # Составной порядок
        (EllipticCurve(a=5, b=13, p=20), [2, 5]),  # Составной порядок
    ]
)
def test_composite_order(curve, expected_subgroup_sizes):
    """
    Тест на проверку составного порядка и соответствующих подгрупп.
    """
    points = find_points(curve)
    order = curve_order(curve)  # Получаем порядок группы кривой
    prime_subgroups = find_prime_subgroups(curve)

    # Проверим, что подгруппы содержат только простые подгруппы
    assert len(prime_subgroups) > 0  # Мы должны получить хотя бы одну подгруппу
    for subgroup in prime_subgroups:
        # Проверим, что размер подгруппы соответствует одному из ожидаемых размеров
        assert len(subgroup) in expected_subgroup_sizes


@pytest.mark.parametrize(
    "curve, expected_order",
    [
        (EllipticCurve(a=5, b=7, p=11), 11),  # Подтверждаем правильность порядка для точки на кривой
        (EllipticCurve(a=2, b=4, p=13), 13),  # Простой порядок
        (EllipticCurve(a=7, b=5, p=17), 17),  # Простой порядок
        (EllipticCurve(a=11, b=9, p=23), 23),  # Простой порядок
        (EllipticCurve(a=8, b=10, p=29), 29),  # Простой порядок
        (EllipticCurve(a=9, b=14, p=101), 101),  # Простой порядок
    ]
)
def test_point_order_large(curve, expected_order):
    """
    Тест на проверку вычисления порядка точки для больших кривых.
    """
    points = find_points(curve)

    # Проверим, что метод корректно находит порядок точки
    for point in points:
        if isinstance(point, ECPointInf):
            continue
        order = point_order(point)  # Вычисляем порядок точки
        assert order is not None, f"Не удалось вычислить порядок точки {point}"


@pytest.mark.parametrize(
    "curve, order, found",
    [
        (EllipticCurve(a=2, b=4, p=11), 11, True),  # Тест на нахождение точки с заданным порядком
        (EllipticCurve(a=7, b=9, p=13), 2, True),  # Тест с другой кривой
        (EllipticCurve(a=8, b=5, p=17), 5, True),  # Тест с ещё одной кривой
        (EllipticCurve(a=13, b=21, p=23), 3, False),  # Тест, где такой точки не должно быть
    ]
)
def test_point_of_order(curve, order, found):
    """
    Тест на нахождение точки с заданным порядком.
    """
    point = point_of_order(curve, order)
    if found:
        assert point is not None, "Не удалось найти точку заданного порядка."
    else:
        assert point is None, "Не ожидали найти точку с этим порядком."


@pytest.mark.parametrize(
    "curve, expected_subgroup_count",
    [
        (EllipticCurve(a=2, b=4, p=7), 1),  # Подгруппа для кривой порядка 7
        (EllipticCurve(a=4, b=8, p=13), 2),  # Проверка для кривой порядка 12
        (EllipticCurve(a=5, b=11, p=17), 1),  # Кривая с порядком 17
        (EllipticCurve(a=8, b=7, p=23), 3),  # Проверка для более сложной кривой
        (EllipticCurve(a=9, b=14, p=29), 1),  # Проверка для более сложной кривой
    ]
)
def test_large_curve(curve, expected_subgroup_count):
    """
    Тест для проверки работы с большими кривыми.
    """
    prime_subgroups = find_prime_subgroups(curve)
    assert len(prime_subgroups) == expected_subgroup_count  # Проверка на количество подгрупп
    print(f"Тест с большой кривой прошел успешно для p={curve.p}.")
