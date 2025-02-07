class ECPointInf:
    """Класс для бесконечно удаленной точки"""
    def __init__(self, curve):
        self.curve = curve

    def __add__(self, other):
        return other

    def __radd__(self, other):
        return other

    def __mul__(self, other):
        return ECPointInf(self.curve)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        return isinstance(other, ECPointInf)

    def __repr__(self):
        return "inf"
    