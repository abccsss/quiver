from vector import *
from quiver import *
from sympy import pretty_print, simplify


def main():
    quiver_a1 = SymmetricQuiver(
        num_vertices=2,
        edges=[[0, 1], [0, 1]],
        vertex_involution=[1, 0],
        edge_involution=[0, 1])

    print('A1 tilde quiver, orthogonal invariants:')
    print()
    for i in range(6):
        pretty_print(simplify(quiver_a1.epsilon_sd([1, -1], [i, i])))
        print()

    quiver_a1_symp = SymmetricQuiver(
        num_vertices=2,
        edges=[[0, 1], [0, 1]],
        vertex_involution=[1, 0],
        edge_involution=[0, 1],
        vertex_sign=[-1, -1])

    print('A1 tilde quiver, symplectic invariants:')
    print()
    for i in range(6):
        pretty_print(simplify(quiver_a1_symp.epsilon_sd([1, -1], [i, i])))
        print()

    quiver_one_loop = SymmetricQuiver(
        num_vertices=1,
        edges=[[0, 0]],
        vertex_involution=[0],
        edge_involution=[0])

    print('One-loop quiver, orthogonal invariants:')
    print()
    for i in range(10):
        pretty_print(simplify(quiver_one_loop.epsilon_sd([0], [i])))
        print()

    quiver_one_loop_symp = SymmetricQuiver(
        num_vertices=1,
        edges=[[0, 0]],
        vertex_involution=[0],
        edge_involution=[0],
        vertex_sign=[-1])

    print('One-loop quiver, symplectic invariants:')
    print()
    for i in range(10):
        pretty_print(simplify(quiver_one_loop_symp.epsilon_sd([0], [i])))
        print()

    quiver_two_loops = SymmetricQuiver(
        num_vertices=1,
        edges=[[0, 0], [0, 0]],
        vertex_involution=[0],
        edge_involution=[0, 1])

    print('Two-loop quiver, orthogonal invariants:')
    print()
    for i in range(10):
        pretty_print(simplify(quiver_two_loops.epsilon_sd([0], [i])))
        print()


main()
