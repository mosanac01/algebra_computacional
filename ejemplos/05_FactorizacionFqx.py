"""
Ejemplo 5: Factorización Cantor–Zassenhaus en F_q[x].

- Construcción de F_25 = F_5[a]/<a^2+2>.
- Factorización de polinomios con coeficientes en F_25.
- Comparación con factorización sobre F_5.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from pyalgcomputacional.cuerpo_fp import cuerpo_fp
from pyalgcomputacional.cuerpo_fq import cuerpo_fq
from pyalgcomputacional.anillo_fp_x import anillo_fp_x
from pyalgcomputacional.anillo_fq_x import anillo_fq_x
from pyalgcomputacional.factorizacion import (
    sqfree_fact_fqx, didegr_fact_fqx, fact_fqx,
)


def separador(titulo):
    print("\n" + "=" * 60)
    print(titulo)
    print("=" * 60)


# ---------------------------------------------------------------
# Construcción de F_25
# ---------------------------------------------------------------
F5 = cuerpo_fp(5)
F5a = anillo_fp_x(F5, var='a')

# a^2 + 2 es irreducible sobre F_5  (x^2 = -2 = 3, 3 no es cuadrado mod 5)
g = F5a.elem_de_str("a^2 + 2")
F25 = cuerpo_fq(F5, g, var='a')
print(f"Cuerpo base: {F25}")
print(f"q = {F25.q}, p = {F25.p}, n = {F25.n}")

F25x = anillo_fq_x(F25, var='x')
print(f"Anillo:      {F25x}")


# ---------------------------------------------------------------
# 1. Test de irreducibilidad en F_25[x]
# ---------------------------------------------------------------
separador("1) Irreducibilidad en F_25[x]")

candidatos = [
    "x^2 + 1",
    "x^2 + a",
    "x^2 + a*x + 1",
    "x^3 + x + 1",
    "x^3 + a*x + 1",
]
for s in candidatos:
    f = F25x.elem_de_str(s)
    print(f"  {s:20s}  irreducible = {F25x.es_irreducible(f)}")


# ---------------------------------------------------------------
# 2. Square-free factorization en F_25[x]
# ---------------------------------------------------------------
separador("2) Square-free factorization en F_25[x]")

# f = (x + a)^2 * (x^2 + a)
p1 = F25x.elem_de_str("x + a")
p2 = F25x.elem_de_str("x^2 + a")

# Multiplicamos: p1 * p1 * p2
p1_cuad = F25x.mult(p1, p1)
f = F25x.mult(p1_cuad, p2)

print(f"f = {F25x.conv_a_str(f)}")

g = sqfree_fact_fqx(F25x, f)
print(f"sqfree(f) = {F25x.conv_a_str(g)}")


# ---------------------------------------------------------------
# 3. Factorización completa en F_25[x]
# ---------------------------------------------------------------
separador("3) Factorización completa en F_25[x]")

# f = (x + a)^3 * (x^2 + a) * (x^2 + a*x + 1)
p1 = F25x.elem_de_str("x + a")
p2 = F25x.elem_de_str("x^2 + a")
p3 = F25x.elem_de_str("x^2 + a*x + 1")

# Multiplicamos: p1^3 * p2 * p3
p1_cubo = F25x.mult(F25x.mult(p1, p1), p1)
f = F25x.mult(F25x.mult(p1_cubo, p2), p3)

print(f"f = {F25x.conv_a_str(f)}")
print(f"grado(f) = {F25x.grado(f)}")

fact = fact_fqx(F25x, f)
print("\nFactorización:")
for u, e in fact:
    print(f"  ({F25x.conv_a_str(u)})^{e}")

# Comprobación
prod = F25x.uno()
for u, e in fact:
    for _ in range(e):
        prod = F25x.mult(prod, u)
lead = f[-1]
prod = F25x.mult_por_escalar(prod, lead)
print(f"\n¿Reconstruye f (salvo constante)? {F25x.es_igual(prod, f)}")