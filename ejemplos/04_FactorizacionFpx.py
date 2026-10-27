"""
Ejemplo 4: Factorización Cantor–Zassenhaus en F_p[x].

- Test de irreducibilidad (Rabin).
- Square-free, distinct-degree y equal-degree.
- Factorización completa con multiplicidades.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from pyalgcomputacional.cuerpo_fp import cuerpo_fp
from pyalgcomputacional.anillo_fp_x import anillo_fp_x
from pyalgcomputacional.factorizacion import (
    sqfree_fact_fpx, didegr_fact_fpx, eqdegr_fact_fpx,
    multiplicidad_fpx, fact_fpx,
)


def separador(titulo):
    print("\n" + "=" * 60)
    print(titulo)
    print("=" * 60)


F = cuerpo_fp(101)
Fx = anillo_fp_x(F, var='x')


# ---------------------------------------------------------------
# 1. Test de irreducibilidad
# ---------------------------------------------------------------
separador("1) Test de irreducibilidad de Rabin en F_101[x]")

candidatos = [
    "x^2 + 1",
    "x^2 + 2",           # -2 no es cuadrado mod 101 => irreducible
    "x^3 + x + 1",
    "x^4 + x + 1",
    "x^5 + 2*x + 3",
    "x^6 + x + 1",
]
for s in candidatos:
    f = Fx.elem_de_str(s)
    print(f"  {s:20s}  irreducible = {Fx.es_irreducible(f)}")


# ---------------------------------------------------------------
# 2. Square-free factorization
# ---------------------------------------------------------------
separador("2) Square-free factorization")

# f = (x+1)^2 * (x^2+1) * (x+3)
f = Fx.elem_de_str("x^5 + 5*x^4 + 10*x^3 + 10*x^2 + 5*x + 1")   # (x+1)^5
g = sqfree_fact_fpx(Fx, f)
print(f"f  = {Fx.conv_a_str(f)}")
print(f"sqfree(f) = {Fx.conv_a_str(g)}")


# ---------------------------------------------------------------
# 3. Distinct-degree factorization
# ---------------------------------------------------------------
separador("3) Distinct-degree factorization")

# Producto de irreducibles de grados 1, 2, 3
f1 = Fx.elem_de_str("x + 1")               # grado 1
f2 = Fx.elem_de_str("x^2 + 2")             # grado 2 irreducible
f3 = Fx.elem_de_str("x^3 + x + 1")         # grado 3 irreducible
g = Fx.mult(Fx.mult(f1, f2), f3)
print(f"g = f1*f2*f3 = {Fx.conv_a_str(g)}")

hlist = didegr_fact_fpx(Fx, g)
for i, h in enumerate(hlist):
    print(f"  bloque grado {i+1}: {Fx.conv_a_str(h)}")


# ---------------------------------------------------------------
# 4. Equal-degree factorization
# ---------------------------------------------------------------
separador("4) Equal-degree factorization (grado 2)")

h1 = Fx.elem_de_str("x^2 + 2")
h2 = Fx.elem_de_str("x^2 + 3")
h3 = Fx.elem_de_str("x^2 + x + 1")
h = Fx.mult(Fx.mult(h1, h2), h3)

print(f"h = {Fx.conv_a_str(h)}")
factores = eqdegr_fact_fpx(Fx, 2, h)
for u in factores:
    print(f"  factor irreducible grado 2: {Fx.conv_a_str(u)}")


# ---------------------------------------------------------------
# 5. Factorización completa con multiplicidades
# ---------------------------------------------------------------
separador("5) Factorización completa de Cantor–Zassenhaus")

p1 = Fx.elem_de_str("x + 1")
p2 = Fx.elem_de_str("x^2 + 2")
p3 = Fx.elem_de_str("x^3 + x + 1")

p1_cubo = Fx.mult(Fx.mult(p1, p1), p1)
p2_cuad = Fx.mult(p2, p2)
f = Fx.mult(Fx.mult(p1_cubo, p2_cuad), p3)

print(f"f = {Fx.conv_a_str(f)}")
print(f"grado(f) = {Fx.grado(f)}")

fact = fact_fpx(Fx, f)
print("\nFactorización:")
for u, e in fact:
    print(f"  ({Fx.conv_a_str(u)})^{e}")

prod = Fx.uno()
for u, e in fact:
    for _ in range(e):
        prod = Fx.mult(prod, u)
lead = f[-1]
prod = Fx.mult_por_escalar(prod, lead)
print(f"\n¿Reconstruye f (salvo constante)? {Fx.es_igual(prod, f)}")