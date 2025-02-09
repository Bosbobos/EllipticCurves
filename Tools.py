import math
from ECPoint import ECPoint
from ECPointInf import ECPointInf


def point_neg(P):
    """
    Возвращает отрицание точки P.
    Для точки в бесконечности возвращается она же.
    """
    if isinstance(P, ECPointInf):
        return P
    return ECPoint(P.curve, P.x, (-P.y) % P.curve.p)


def legendre_symbol(a, p):
    """
       Вычисляет символ Лежандра (a / p), который равен:
       -1, если a является квадратичным вычетом по модулю p,
       1, если a не является квадратичным вычетом по модулю p,
       и 0, если a делится на p.

       Аргументы:
       a -- Целое число, для которого вычисляется символ Лежандра
       p -- Простое число, модуль

       Возвращает:
       1, -1 или 0 в зависимости от значения символа Лежандра.
       """
    ls = pow(a, (p - 1) // 2, p)
    if ls == p - 1:
        return -1
    return ls


def tonelli_shanks(n, p):
    """
      Алгоритм Тонелли-Шэнкса для нахождения квадратичного корня по модулю простого числа p
      при условии, что корень существует (символ Лежандра равен 1).

      Аргументы:
      n -- Число, для которого нужно найти квадратичный корень
      p -- Простое число, модуль

      Возвращает:
      Квадратичный корень n по модулю p, если он существует, иначе None.
      """
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
    points = [ECPointInf(curve)]
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


def bsgs(curve, P, Q, n=None):
    """
    Решает задачу дискретного логарифма в циклической подгруппе,
    порождённой точкой P, то есть ищет такое целое d, что
        d * P = Q.
    Если порядок n подгруппы не задан (n=None), то он вычисляется наивно.

    Алгоритм:
      1. Выбирается m = ceil(sqrt(n)).
      2. Вычисляются baby-шаги: для j от 0 до m-1 сохраняется значение j*P.
      3. Выполняются giant-шаги: ищется такое i, что
             Q - i*(m*P)
         встречается среди baby-шагов.
         Тогда d = i*m + j.
    """
    # Если порядок n не задан, вычисляем его наивно (подойдёт для небольших кривых).
    if n is None:
        current = P
        n = 1
        while not isinstance(current, ECPointInf):
            current = current + P
            n += 1
            # Ограничение для предотвращения зацикливания на больших группах.
            if n > curve.p + 5:
                break

    m = int(math.ceil(math.sqrt(n)))

    # Baby-шаги: для j от 0 до m-1 вычисляем j * P.
    baby_steps = {}
    for j in range(m):
        point = j * P
        key = 'inf' if isinstance(point, ECPointInf) else (point.x, point.y)
        baby_steps[key] = j

    mP = m * P
    neg_mP = point_neg(mP)

    # Giant-шаги: ищем i от 0 до m-1, для которого
    # Q - i*(mP) встречается в baby_steps.
    gamma = Q
    for i in range(m):
        key = 'inf' if isinstance(gamma, ECPointInf) else (gamma.x, gamma.y)
        if key in baby_steps:
            j = baby_steps[key]
            d = i * m + j
            return d % n  # возвращаем наименьшее неотрицательное решение
        gamma = gamma + neg_mP  # эквивалентно: gamma = gamma - mP

    return None  # если решение не найдено


def curve_order(curve, max_trials=10):
    points = find_points(curve)
    #if len(points) < 2 ** 10 or len(points) > 2 ** 512:
    #    raise ValueError("Порядок кривой не входит в требуемый диапазон.")
    return len(points)


def find_prime_subgroups_orders(curve):
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
            if not isinstance(point, ECPointInf):
                candidate = (order // p) * point
                if isinstance(candidate, ECPointInf):
                    subgroups.append(p)
                    break
    return subgroups


def find_prime_subgroups(curve):
    """
    Нахождение простых подгрупп кривой.
    Возвращает список подгрупп, каждая из которых состоит из точек данной кривой,
    которые соответствуют простым делителям порядка кривой.
    """
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
        if p == order:  # Исключаем порядок, равный самому порядку кривой, если он не является простым
            continue
        subgroup = [ECPointInf(curve)]  # Подгруппа обязательно включает точку в бесконечности
        for point in find_points(curve):
            if isinstance(point, ECPointInf):
                continue
            if point_order(point) == p:
                subgroup.append(point)
        if len(subgroup) == p:  # В подгруппе должно быть p точек, включая точку в бесконечности
            subgroups.append(subgroup)
    return subgroups


def point_order(P):
    """
    Нахождение порядка точки P на эллиптической кривой.
    Порядок точки - минимальное целое число k, такое что k * P = O (точка в бесконечности).
    Если точка P имеет бесконечный порядок, возвращается None.
    """
    if isinstance(P, ECPointInf):
        return 1
    original = P
    k = 1
    while True:
        P = P + original
        k += 1
        if isinstance(P, ECPointInf):
            return k
#        if k > 10000:  # Ограничение для предотвращения бесконечных циклов
#           break
    return None


def point_of_order(curve, order):
    """
    Нахождение точки заданного порядка на кривой.
    Возвращает точку P на кривой, такую что P * order = O.
    Если такой точки не существует, возвращает None.
    """
    points = find_points(curve)
    for point in points:
        if isinstance(point, ECPointInf):
            continue
        if point_order(point) == order:
            return point
    return None
