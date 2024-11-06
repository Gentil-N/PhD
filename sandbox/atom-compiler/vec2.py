from typing import Self

class Vec2:

    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y

    def add(self, other: Self) -> Self:
        self.x += other.x
        self.y += other.y
        return self

    def sub(self, other: Self) -> Self:
        self.x -= other.x
        self.y -= other.y
        return self

    def mul(self, factor) -> Self:
        self.x *= factor
        self.y *= factor
        return self

    def set(self, other: Self) -> Self:
        self.x = other.x
        self.y = other.y
        return self

    def set_xy(self, x, y) -> Self:
        self.x = x
        self.y = y
        return self

    def fsqlength(self) -> float:
        return float((self.x * self.x) + (self.y * self.y))

    def flength(self) -> float:
        return np.sqrt(self.fsqlength())


def vec2_from_points(vec_a: Vec2, vec_b: Vec2):
    return Vec2(vec_b.x - vec_a.x, vec_b.y - vec_a.y)
