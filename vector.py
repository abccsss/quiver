def vector_add(d1: list[int], d2: list[int]) -> list[int]:
    n = min(len(d1), len(d2))
    result = [0] * n

    for i in range(n):
        result[i] = d1[i] + d2[i]

    return result


def vector_minus(d1: list[int], d2: list[int]) -> list[int]:
    n = min(len(d1), len(d2))
    result = [0] * n

    for i in range(n):
        result[i] = d1[i] - d2[i]

    return result


def vector_leq(d1: list[int], d2: list[int]) -> bool:
    n = min(len(d1), len(d2))
    for i in range(n):
        if d1[i] > d2[i]:
            return False

    return True


def dot_product(d1: list[int], d2: list[int]) -> int:
    result = 0

    for i in range(min(len(d1), len(d2))):
        result += d1[i] * d2[i]

    return result


def all_proper_sub(d: list[int]) -> list[list[int]]:
    result = []

    n = len(d)
    level = n - 1
    current = [0] * n

    while True:
        while current[level] == d[level] and level >= 0:
            level -= 1
        if level == -1:
            break

        current[level] += 1
        for i in range(level + 1, n):
            current[i] = 0
        level = n - 1

        result.append(current.copy())

    if len(result) > 0:
        result.pop()  # Remove d itself

    return result
