"""
Ejemplo 1: Construcción de cuerpos finitos F_p y F_q.

- F_p como Z/pZ.
- F_q = F_p[a]/<g(a)> con g irreducible.
- Operaciones: suma, producto, inverso, potencia.
- Tablas de operaciones y cuadrado latino.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from pyalgcomputacional.cuerpo_fp import cuerpo_fp
from pyalgcomputacional.cuerpo_fq import cuerpo_fq
from pyalgcomputacional.anillo_fp_x import anillo_fp_x


def separador(titulo):
    print("\n" + "=" * 60)
    print(titulo)
    print("=" * 60)


# ---------------------------------------------------------------
# 1. Cuerpo primo F_7
# ---------------------------------------------------------------
separador("1) Cuerpo primo F_7")
F7 = cuerpo_fp(7)
print(f"Cuerpo: {F7}")

a, b = 5, 6
print(f"{a} + {b} = {F7.suma(a, b)}")
print(f"{a} * {b} = {F7.mult(a, b)}")
print(f"-{a}      = {F7.inv_adit(a)}")
print(f"{a}^-1    = {F7.inv_mult(a)}")
print(f"{a}^100   = {F7.pot(a, 100)}")

print("\nTabla de suma en F_7:")
for fila in F7.tabla_suma():
    print(fila)

print("\nTabla de inversos multiplicativos en F_7:")
print(F7.tabla_inv_mult())


# ---------------------------------------------------------------
# 2. Extensión F_49 = F_7[a]/<a^2 + 1>
# ---------------------------------------------------------------
separador("2) Extensión F_49 = F_7[a]/<a^2 + 1>")

F7_anillo = anillo_fp_x(F7, var='a')
g = F7_anillo.elem_de_str("a^2 + 1")   # x^2+1 es irreducible en F_7
F49 = cuerpo_fq(F7, g, var='a')
print(f"Cuerpo: {F49}")
print(f"q = {F49.q}, p = {F49.p}, n = {F49.n}")

# Elementos mediante strings
alpha = F49.elem_de_str("a")       # la clase de a
beta  = F49.elem_de_str("3a + 2")

print(f"alpha       = {F49.conv_a_str(alpha)}")
print(f"beta        = {F49.conv_a_str(beta)}")
print(f"alpha+beta  = {F49.conv_a_str(F49.suma(alpha, beta))}")
print(f"alpha*beta  = {F49.conv_a_str(F49.mult(alpha, beta))}")
print(f"alpha^-1    = {F49.conv_a_str(F49.inv_mult(alpha))}")
print(f"alpha^48    = {F49.conv_a_str(F49.pot(alpha, 48))}   (debe ser 1)")

# Elementos aleatorios
print("\n5 elementos aleatorios de F_49:")
for _ in range(5):
    x = F49.aleatorio()
    print(f"  int={F49.conv_a_int(x):3d}  ->  {F49.conv_a_str(x)}")


# ---------------------------------------------------------------
# 3. Cuadrado latino en F_11
# ---------------------------------------------------------------
separador("3) Cuadrado latino en F_11 con a = 3")
F11 = cuerpo_fp(11)
cl = F11.cuadrado_latino(3)
for fila in cl:
    print(fila)