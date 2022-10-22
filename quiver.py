from sympy import Expr, Integer, Rational, Symbol
from vector import *


def motive_gl(n: int) -> Expr:
    if (n < 0):
        return Integer(0)

    result = Integer(1)
    L = Symbol('L')
    for i in range(n):
        result = result * (L ** n - L ** i)

    return result


def motive_o(n: int) -> Expr:
    if (n < 0):
        return Integer(0)

    result = Integer(1)
    L = Symbol('L')
    half_n = int(n / 2)

    for i in range(half_n):
        result = result * (L ** (2 * half_n) - L ** (2 * i))

    if n % 2 == 0:
        result = result / (L ** half_n)
    else:
        result = result * (L ** half_n)

    return result


def motive_sp(n: int) -> Expr:
    if (n < 0 or n % 2 == 1):
        return Integer(0)

    return motive_o(n + 1)


def epsilon_sd_coeff(n: int) -> Rational:
    p = 1
    q = 1

    for i in range(n):
        p *= -(2 * i + 1)
        q *= 2 * (i + 1)

    return Rational(p, q)


class Quiver:
    num_vertices: int
    edges: list[list[int]]

    def __init__(self, num_vertices: int, edges: list[list[int]]) -> None:
        self.num_vertices = num_vertices
        self.edges = edges

    def chi(self, d1: list[int], d2: list[int]) -> int:
        result = dot_product(d1, d2)

        for e in self.edges:
            result -= d1[e[0]] * d2[e[1]]

        return result

    def chi_multi(self, d: list[list[int]]) -> int:
        if (len(d) == 0):
            return 0

        result = 0
        n = self.num_vertices
        accum = [0] * n

        for i in range(len(d)):
            if i != 0:
                result += self.chi(d[i], accum)

            for j in range(n):
                accum[j] += d[i][j]

        return result

    def delta_zero(self, dim_vector: list[int]) -> Expr:
        freedom = 0
        d = dim_vector
        for e in self.edges:
            freedom += d[e[0]] * d[e[1]]

        L = Symbol('L')
        result = L ** freedom

        for i in range(self.num_vertices):
            result = result / motive_gl(d[i])

        return result

    def enumerate_wcf_terms(self, stability: list[int], dim_vector: list[int], equal_slope: bool = False) -> list[list[list[int]]]:
        results: list[list[list[int]]] = [[dim_vector]]
        n = self.num_vertices
        slope = dot_product(stability, dim_vector) / sum(dim_vector)

        accum: list[list[int]] = [[0] * n]
        all_sub: list[list[list[int]]] = [all_proper_sub(dim_vector)]
        level = 0
        i: list[int] = [0]

        while True:
            if level == -1:
                break

            if i[level] == len(all_sub[level]):
                level -= 1
                i.pop()
                accum.pop()
                all_sub.pop()

                if level >= 0:
                    i[level] += 1
                continue

            # Compute accum[level] + new term, see if it is valid
            new_accum = vector_add(accum[level], all_sub[level][i[level]])
            new_slope = dot_product(stability, new_accum) / sum(new_accum)
            if (not equal_slope and new_slope <= slope) or (equal_slope and new_slope != slope):
                # Term not in WCF
                i[level] += 1
                continue

            accum.append(new_accum)
            remaining = vector_minus(dim_vector, new_accum)

            # Append term [..., new term, remaining] to result
            new_result: list[list[int]] = []
            for j in range(level + 1):
                new_result.append(all_sub[j][i[j]])
            new_result.append(remaining)
            results.append(new_result)

            # Compute all_sub[level + 1], and advance level
            all_sub.append(all_proper_sub(remaining))
            i.append(0)
            level += 1

        return results

    def delta(self, stability: list[int], dim_vector: list[int]) -> Expr:
        result = Integer(0)
        wcf_terms = self.enumerate_wcf_terms(stability, dim_vector)
        L = Symbol('L')

        for term in wcf_terms:
            m = len(term)
            coeff = Integer(-1) ** (m - 1)
            adjust = L ** (-self.chi_multi(term))

            expr = coeff * adjust
            for factor in term:
                expr *= self.delta_zero(factor)

            result += expr

        return result

    def epsilon(self, stability: list[int], dim_vector: list[int]) -> Expr:
        result = Integer(0)
        wcf_terms = self.enumerate_wcf_terms(
            stability, dim_vector, equal_slope=True)
        L = Symbol('L')

        for term in wcf_terms:
            m = len(term)
            coeff = Rational((-1) ** m, m)
            adjust = L ** (-self.chi_multi(term))

            expr = coeff * adjust
            for factor in term:
                expr *= self.delta(stability, factor)

            result += expr

        return result


class SymmetricQuiver(Quiver):
    vertex_involution: list[int]
    edge_involution: list[int]
    vertex_sign: list[int]
    edge_sign: list[int]

    def __init__(
            self, num_vertices: int, edges: list[list[int]],
            vertex_involution: list[int],
            edge_involution: list[int],
            vertex_sign: list[int] | None = None,
            edge_sign: list[int] | None = None) -> None:
        self.num_vertices = num_vertices
        self.edges = edges
        self.vertex_involution = vertex_involution
        self.edge_involution = edge_involution
        self.vertex_sign = vertex_sign or [1] * num_vertices
        self.edge_sign = edge_sign or [1] * len(edges)

    def symmetrize(self, dim_vector: list[int]) -> list[int]:
        n = self.num_vertices
        result = [0] * n

        for i in range(n):
            result[i] = dim_vector[self.vertex_involution[i]] + dim_vector[i]

        return result

    def all_sd_sub(self, dim_vector: list[int]) -> list[list[int]]:
        return list(filter(
            lambda d: vector_leq(self.symmetrize(d), dim_vector),
            all_proper_sub(dim_vector)))

    def chi_sd(self, d1: list[int], d2: list[int]) -> int:
        result = self.chi(d1, d2)
        n = self.num_vertices

        for i in range(n):
            if self.vertex_involution[i] == i:
                result += d2[i] * (d2[i] - self.vertex_sign[i]) / 2
            elif self.vertex_involution[i] > i:
                result += d2[i] * d2[self.vertex_involution[i]]

        j = 0
        for e in self.edges:
            if (self.edge_involution[j] == j):
                result -= d2[e[1]] * (d2[e[1]] + self.edge_sign[j]
                                      * self.vertex_sign[e[1]]) / 2
            elif self.edge_involution[j] > j:
                result -= d2[self.vertex_involution[e[0]]] * d2[e[1]]
            j += 1

        return int(result)

    def chi_sd_multi(self, d: list[list[int]]) -> int:
        if (len(d) == 0):
            return 0

        result = 0
        n = self.num_vertices
        accum = [0] * n

        for i in range(len(d)):
            if i != 0:
                if i != len(d) - 1:
                    result += self.chi(d[i], accum)
                else:
                    result += self.chi_sd(d[i], accum)

            for j in range(n):
                accum[j] += d[i][j]

        return result

    def enumerate_sd_wcf_terms(self, stability: list[int], dim_vector: list[int], equal_slope: bool = False) -> list[list[list[int]]]:
        results: list[list[list[int]]] = [[dim_vector]]
        n = self.num_vertices

        accum: list[list[int]] = [[0] * n]
        all_sub: list[list[list[int]]] = [self.all_sd_sub(dim_vector)]
        level = 0
        i: list[int] = [0]

        while True:
            if level == -1:
                break

            if i[level] == len(all_sub[level]):
                level -= 1
                i.pop()
                accum.pop()
                all_sub.pop()

                if level >= 0:
                    i[level] += 1
                continue

            # Compute accum[level] + new term, see if it is valid
            new_accum = vector_add(accum[level], all_sub[level][i[level]])
            new_slope = dot_product(stability, new_accum) / sum(new_accum)
            if (not equal_slope and new_slope <= 0) or (equal_slope and new_slope != 0):
                # Term not in WCF
                i[level] += 1
                continue

            accum.append(new_accum)
            remaining = vector_minus(dim_vector, self.symmetrize(new_accum))

            # Append term [..., new term, remaining] to result
            new_result: list[list[int]] = []
            for j in range(level + 1):
                new_result.append(all_sub[j][i[j]])
            new_result.append(remaining)
            results.append(new_result)

            # Compute all_sub[level + 1], and advance level
            all_sub.append(self.all_sd_sub(remaining))
            i.append(0)
            level += 1

        return results

    def delta_sd_zero(self, dim_vector: list[int]) -> Expr:
        freedom: int = 0
        d = dim_vector

        j = 0
        for e in self.edges:
            if self.edge_involution[j] == j:
                freedom += int(d[e[1]] * (d[e[1]] + self.edge_sign[j]) / 2)
            elif self.edge_involution[j] > j:
                freedom += d[e[0]] * d[e[1]]
            j += 1

        L = Symbol('L')
        result = L ** freedom

        for i in range(self.num_vertices):
            if self.vertex_involution[i] == i:
                if self.vertex_sign[i] == 1:
                    result = result / motive_o(d[i])
                else:
                    if d[i] % 2 != 0:
                        return Integer(0)
                    result = result / motive_sp(d[i])
            elif self.vertex_involution[i] > i:
                result = result / motive_gl(d[i])

        return result

    def delta_sd(self, stability: list[int], dim_vector: list[int]) -> Expr:
        result = Integer(0)
        wcf_terms = self.enumerate_sd_wcf_terms(stability, dim_vector)
        L = Symbol('L')

        for term in wcf_terms:
            m = len(term)
            coeff = Integer(-1) ** (m - 1)
            adjust = L ** (-self.chi_sd_multi(term))

            expr = coeff * adjust
            for i in range(len(term)):
                if i != len(term) - 1:
                    expr *= self.delta_zero(term[i])
                else:
                    expr *= self.delta_sd_zero(term[i])

            result += expr

        return result

    def epsilon_sd(self, stability: list[int], dim_vector: list[int]) -> Expr:
        result = Integer(0)
        wcf_terms = self.enumerate_sd_wcf_terms(
            stability, dim_vector, equal_slope=True)
        L = Symbol('L')

        for term in wcf_terms:
            m = len(term)
            coeff = epsilon_sd_coeff(m - 1)
            adjust = L ** (-self.chi_sd_multi(term))

            expr = coeff * adjust
            for i in range(len(term)):
                if i != len(term) - 1:
                    expr *= self.delta(stability, term[i])
                else:
                    expr *= self.delta_sd(stability, term[i])

            result += expr

        return result
