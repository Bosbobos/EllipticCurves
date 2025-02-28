from ECPoint import ECPoint
from ECPointInf import ECPointInf
from EllipticCurve import EllipticCurve
from Tools import find_points, is_prime, find_prime_subgroups, curve_order, point_order, point_of_order


def introduce_curve():
    # Ввод данных кривой в одной строке
    p, a, b = map(int, input("Введите параметры p, a, b через пробел: ").split())

    if not is_prime(p):
        print('p должен быть простым')
        return

    # Создание объекта кривой
    curve = EllipticCurve(p, a, b)

    # Получаем и выводим точки кривой
    points = find_points(curve)
    print(f"Точки на кривой:")
    for point in points:
        print(f"Точка: {point}")

    # Вычисление порядка кривой
    order = curve_order(curve)
    print(f"Порядок кривой: {order}")

    return curve


def compute_multiple(curve):
    # Ввод точки на кривой в одной строке
    x_P, y_P = map(int, input("Введите координаты точки P (x, y) через пробел: ").split())

    # Ввод кратности для вычисления точки
    scalar = int(input("Введите кратность для вычисления P * k: "))

    P = ECPoint(curve, x_P, y_P)
    result_point = scalar * P
    print(f"Точка P * {scalar} = {result_point}")


def find_subgroups(curve):
    # Находим и выводим подгруппы простого порядка
    subgroups = find_prime_subgroups(curve)
    print("Подгруппы простого порядка:")
    for subgroup in subgroups:
        print(f"Подгруппа: {[str(point) for point in subgroup]}")


def compute_order(curve):
    # Ввод точки на кривой в одной строке
    x_P, y_P = map(int, input("Введите координаты точки P (x, y) через пробел: ").split())

    # Создание объекта точки
    P = ECPoint(curve, x_P, y_P)

    # Вычисление порядка точки
    order = point_order(P)
    print(f"Порядок точки P = {order}")


def find_point_by_order(curve):
    # Ввод порядка точки
    order = int(input("Введите порядок точки: "))

    # Находим точку с заданным порядком
    result_point = point_of_order(curve, order)
    print(f"Точка с порядком {order}: {result_point}")


def add_points(curve):
    # Ввод двух точек на кривой в одной строке
    x_P1, y_P1 = map(int, input("Введите координаты первой точки P1 (x, y) через пробел: ").split())
    x_P2, y_P2 = map(int, input("Введите координаты второй точки P2 (x, y) через пробел: ").split())

    # Создание объектов точек
    P1 = ECPoint(curve, x_P1, y_P1)
    P2 = ECPoint(curve, x_P2, y_P2)

    # Складываем точки
    result_point = P1 + P2
    print(f"Сумма точек P1 + P2 = {result_point}")


def main():
    print("Добро пожаловать в программу для работы с эллиптическими кривыми!")

    curve = introduce_curve()

    while True:
        print("\nВыберите операцию:")
        print("1. Построить кривую и вывести все точки, вычислить её порядок.")
        print("2. Вычислить точку кратности для заданной точки P.")
        print("3. Найти подгруппы простого порядка.")
        print("4. Вывести порядок заданной точки.")
        print("5. Найти точку заданного порядка.")
        print("6. Сложить две точки.")
        print("7. Выход.")

        choice = input("Ваш выбор: ")

        if choice == "1":
            curve = introduce_curve()
        elif choice == "2":
            if curve:
                compute_multiple(curve)
            else:
                print("Сначала создайте кривую.")
        elif choice == "3":
            if curve:
                find_subgroups(curve)
            else:
                print("Сначала создайте кривую.")
        elif choice == "4":
            if curve:
                compute_order(curve)
            else:
                print("Сначала создайте кривую.")
        elif choice == "5":
            if curve:
                find_point_by_order(curve)
            else:
                print("Сначала создайте кривую.")
        elif choice == "6":
            if curve:
                add_points(curve)
            else:
                print("Сначала создайте кривую.")
        elif choice == "7":
            print("Выход из программы...")
            break
        else:
            print("Неверный выбор. Попробуйте снова.")


if __name__ == "__main__":
    main()
