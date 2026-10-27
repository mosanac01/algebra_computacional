"""
Ejemplo 2: Anillos de polinomios F_p[x] y F_q[x].

- Construcción desde strings, tuplas o enteros.
- Operaciones aritméticas.
- GCD y gcd extendido.
- Inverso modular y potencias modulares.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from pyalgcomputacional.cuerpo_fp import cuerpo_fp
from pyalgcomputacional.cuerpo_fq import cuerpo_fq
from pyalgcomputacional.anillo_fp_x import anillo_fp_x
from pyalgcomputacional.anillo_fq_x import anillo_fq_x


def separador(titulo):
    print("\n" + "=" * 60)
    print(titulo)
    print("=" * 60)


# ---------------------------------------------------------------
# 1. Anillo F_5[x]
# ---------------------------------------------------------------
separador("1) Anillo F_5[x]")
F5 = cuerpo_fp(5)
F5x = anillo_fp_x(F5, var='x')
print(f"Anillo: {F5x}")

f = F5x.elem_de_str("2*x^3 + 3*x + 4")
g = F5x.elem_de_str("x^2 + 1")
print(f"f = {F5x.conv_a_str(f)}")
print(f"g = {F5x.conv_a_str(g)}")

print(f"f + g       = {F5x.conv_a_str(F5x.suma(f, g))}")
print(f"f * g       = {F5x.conv_a_str(F5x.mult(f, g))}")
print(f"-f          = {F5x.conv_a_str(F5x.inv_adit(f))}")
print(f"grado(f)    = {F5x.grado(f)}")

q, r = F5x.divmod(f, g)
print(f"f = g*q + r con q = {F5x.conv_a_str(q)}, r = {F5x.conv_a_str(r)}")

print(f"gcd(f, g)   = {F5x.conv_a_str(F5x.gcd(f, g))}")


# ---------------------------------------------------------------
# 2. gcd extendido y inverso modular
# ---------------------------------------------------------------
separador("2) gcd extendido e inverso modular en F_5[x]")

# f = x^2 + 1 y h = x^2 + x + 1 son coprimos en F_5[x]
h = F5x.elem_de_str("x^2 + x + 1")
d, s, t = F5x.gcd_ext(f, h)
print(f"gcd_ext(f, h) = (d, s, t) con d = s*f + t*h")
print(f"  d = {F5x.conv_a_str(d)}")
print(f"  s = {F5x.conv_a_str(s)}")
print(f"  t = {F5x.conv_a_str(t)}")

inv = F5x.inv_mod(f, h)
print(f"\ninverso de f mod h = {F5x.conv_a_str(inv)}")
comprobacion = F5x.mod(F5x.mult(f, inv), h)
print(f"f * inv mod h      = {F5x.conv_a_str(comprobacion)}   (debe ser 1)")

# Potencia modular
pot = F5x.pot_mod(f, 123, h)
print(f"f^123 mod h        = {F5x.conv_a_str(pot)}")


# ---------------------------------------------------------------
# 3. Anillo F_9[y] con F_9 = F_3[a]/<a^2+1>
# ---------------------------------------------------------------
separador("3) Anillo F_9[y] con F_9 = F_3[a]/<a^2+1>")

F3 = cuerpo_fp(3)
F3a = anillo_fp_x(F3, var='a')
g9 = F3a.elem_de_str("a^2 + 1")
F9 = cuerpo_fq(F3, g9, var='a')
print(f"Cuerpo base: {F9}")

F9y = anillo_fq_x(F9, var='y')
print(f"Anillo:      {F9y}")

# Polinomio con coeficientes en F_9
p1 = F9y.elem_de_str("(a + 1)*y^2 + (2a)*y + a + 2")
p2 = F9y.elem_de_str("y + 1")
print(f"p1 = {F9y.conv_a_str(p1)}")
print(f"p2 = {F9y.conv_a_str(p2)}")
print(f"p1 + p2 = {F9y.conv_a_str(F9y.suma(p1, p2))}")
print(f"p1 * p2 = {F9y.conv_a_str(F9y.mult(p1, p2))}")

q, r = F9y.divmod(p1, p2)
print(f"p1 = p2*q + r con q = {F9y.conv_a_str(q)}, r = {F9y.conv_a_str(r)}")
print(f"gcd(p1, p2) = {F9y.conv_a_str(F9y.gcd(p1, p2))}")