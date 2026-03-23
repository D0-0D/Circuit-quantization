import numpy
import scqubits as scq


zp_yaml = """# zero-pi
branches:
- [JJ, 0, 1, EJ]
- [L, 1,2, L1]
- [L, 0, 2, L2]
- [C, 0,1, C1]
- [C, 1,4, C2]
"""

zero_pi = scq.Circuit(zp_yaml, from_file=False)

print(zero_pi.sym_hamiltonian())

"""my_components = [
        Component(0, 1, 'JJ', 'J'),
        Component(1, 2, 'L', 'L1'),
        Component(2, 0, 'L', 'L2'),
        Component(1, 0, 'C', 2e-6),
    ]"""