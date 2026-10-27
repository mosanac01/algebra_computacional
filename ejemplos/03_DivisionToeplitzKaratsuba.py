"""
Ejemplo 3: Multiplicación Karatsuba y división polinómica vía Toeplitz.

- Comparación de Karatsuba vs multiplicación clásica.
- División rápida usando matrices de Toeplitz.
- Inversión de matrices de Toeplitz triangulares.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

import time
from pyalgcomputacional.cuerpo_fp import cuerpo_fp
from pyalgcomputacional.anillo_fp_x import anillo_fp_x
from pyalgcomputacional.karatsuba import fp_x_mult_karatsuba
from pyalgcomputacional.division import fp_x_divmod
from pyalgcomputacional.toeplitz import (
    fp_toep_inf_vec, fp_toep_sup_vec,
    fp_toep_inf_inv, fp_toep_sup_inv,
)


def separador(titulo):
    print("\n" + "=" * 60)
    print(titulo)
    print("=" * 60)


F = cuerpo_fp(65537)     # primo de Fermat: FFT natural, buen tamaño
Fx = anillo_fp_x(F, var='x')


# ---------------------------------------------------------------
# 1. Karatsuba vs multiplicación clásica
# ---------------------------------------------------------------
separador("1) Karatsuba vs multiplicación clásica en F_65537[x]")

n = 500
# Coeficientes pseudoaleatorios deterministas
import random
random.seed(42)
coefs_f = tuple(random.randint(0, F.p - 1) for _ in range(n + 1))
coefs_g = tuple(random.randint(0, F.p - 1) for _ in range(n + 1))

f = Fx.elem_de_tuple(coefs_f)
g = Fx.elem_de_tuple(coefs_g)

t0 = time.perf_counter()
h_clas = Fx.mult(f, g)
t_clas = time.perf_counter() - t0

t0 = time.perf_counter()
h_kar = fp_x_mult_karatsuba(Fx, f, g)
t_kar = time.perf_counter() - t0

print(f"grado(f) = {Fx.grado(f)}, grado(g) = {Fx.grado(g)}")
print(f"Tiempo clásico   : {t_clas*1000:8.2f} ms")
print(f"Tiempo Karatsuba : {t_kar*1000:8.2f} ms")
print(f"¿Resultados iguales? {Fx.es_igual(h_clas, h_kar)}")


# ---------------------------------------------------------------
# 2. División rápida vía Toeplitz
# ---------------------------------------------------------------
separador("2) División polinómica vía Toeplitz")

f = Fx.elem_de_str("x^8 + 3*x^6 + 2*x^4 + 7*x^2 + 1")
g = Fx.elem_de_str("x^3 + 2*x + 5")
print(f"f = {Fx.conv_a_str(f)}")
print(f"g = {Fx.conv_a_str(g)}")

q_rap, r_rap = fp_x_divmod(Fx, f, g)
q_clas, r_clas = Fx.divmod(f, g)

print(f"q (Toeplitz) = {Fx.conv_a_str(q_rap)}")
print(f"r (Toeplitz) = {Fx.conv_a_str(r_rap)}")
print(f"¿Coincide con divmod clásico? "
      f"{Fx.es_igual(q_rap, q_clas) and Fx.es_igual(r_rap, r_clas)}")


# ---------------------------------------------------------------
# 3. Inversión de Toeplitz triangular
# ---------------------------------------------------------------
separador("3) Inversión de matrices de Toeplitz triangular")

a_inf = (2, 3, 5, 7, 11)   # primera columna de T inferior
inv_inf = fp_toep_inf_inv(F, len(a_inf), a_inf)
print(f"Primera columna T_inf     : {a_inf}")
print(f"Primera columna T_inf^-1  : {inv_inf}")

# Comprobamos T * T^-1 = I  (en los primeros n coeficientes)
n = len(a_inf)
col_0 = (1,) + (0,) * (n - 1)
prod = fp_toep_inf_vec(F, n, a_inf, inv_inf)
print(f"T_inf * T_inf^-1 (col 1)  : {prod}")
print(f"¿Es la identidad?          {prod == col_0}")

a_sup = (2, 3, 5, 7, 11)
inv_sup = fp_toep_sup_inv(F, len(a_sup), a_sup)
print(f"\nPrimera fila T_sup        : {a_sup}")
print(f"Primera fila T_sup^-1     : {inv_sup}")

prod_sup = fp_toep_inf_vec(F, n, a_sup, inv_sup)
print(f"T_sup * T_sup^-1          : {prod_sup}")
print(f"¿Es la identidad?          {prod_sup == col_0}")