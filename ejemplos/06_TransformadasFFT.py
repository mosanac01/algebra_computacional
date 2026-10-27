"""
Ejemplo 6: FFT (Cooley–Tukey) sobre F_p y F_q.

- Transformada discreta de Fourier sobre F_65537 (raíz primitiva de orden 2^k).
- Inversa y comprobación de la identidad.
- Convolución rápida mediante FFT.
- FFT sobre una extensión F_q.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from pyalgcomputacional.cuerpo_fp import cuerpo_fp
from pyalgcomputacional.cuerpo_fq import cuerpo_fq
from pyalgcomputacional.anillo_fp_x import anillo_fp_x
from pyalgcomputacional.transformadas import fp_fft, fp_ifft, fq_fft, fq_ifft


def separador(titulo):
    print("\n" + "=" * 60)
    print(titulo)
    print("=" * 60)


# ---------------------------------------------------------------
# 1. FFT sobre F_65537
# ---------------------------------------------------------------
separador("1) FFT sobre F_65537")

F = cuerpo_fp(65537)          # 65536 = 2^16 divide a p-1
n = 8
k = 3                         # n = 2^3

# 3 es raíz primitiva de F_65537; g = 3^(65536/n) tiene orden n
g = F.pot(3, (F.p - 1) // n)
print(f"n = {n}, g = {g}")
print(f"g^n = {F.pot(g, n)}  (debe ser 1)")
print(f"g^(n/2) = {F.pot(g, n // 2)}  (debe ser -1)")

a = tuple(range(1, n + 1))
print(f"a = {a}")

A = fp_fft(F, g, k, a)
print(f"FFT(a)  = {A}")

a_rec = fp_ifft(F, g, k, A)
print(f"IFFT(A) = {a_rec}")
print(f"¿Recupera a? {a_rec == a}")


# ---------------------------------------------------------------
# 2. Convolución vía FFT
# ---------------------------------------------------------------
separador("2) Convolución rápida vía FFT")

# Elegimos n como potencia de 2 mayor o igual que len(f)+len(g)-1
f = (1, 2, 3, 4)
g_poly = (5, 6, 7, 0)
n_conv = 8
k_conv = 3

# Rellenamos con ceros
f_ext = f + (0,) * (n_conv - len(f))
g_ext = g_poly + (0,) * (n_conv - len(g_poly))

g_root = F.pot(3, (F.p - 1) // n_conv)

Ff = fp_fft(F, g_root, k_conv, f_ext)
Gg = fp_fft(F, g_root, k_conv, g_ext)
H  = tuple(F.mult(Ff[i], Gg[i]) for i in range(n_conv))
h  = fp_ifft(F, g_root, k_conv, H)

print(f"f = {f}")
print(f"g = {g_poly}")
print(f"Convolución vía FFT = {h}")

# Convolución clásica
h_clas = [0] * (len(f) + len(g_poly) - 1)
for i, fi in enumerate(f):
    for j, gj in enumerate(g_poly):
        h_clas[i + j] = F.suma(h_clas[i + j], F.mult(fi, gj))
print(f"Convolución clásica = {tuple(h_clas)}")
print(f"¿Coinciden? {tuple(h[:len(h_clas)]) == tuple(h_clas)}")


# ---------------------------------------------------------------
# 3. FFT sobre F_17^2  (extensión de F_17)
# ---------------------------------------------------------------
separador("3) FFT sobre una extensión F_q")

# Construimos F_289 = F_17[a]/<a^2+1>  (17 ≡ 1 mod 4 => -1 es cuadrado,
# así que a^2+1 no sirve; usamos a^2+3 que sí es irreducible)
F17 = cuerpo_fp(17)
F17a = anillo_fp_x(F17, var='a')
g_irr = F17a.elem_de_str("a^2 + 3")
F289 = cuerpo_fq(F17, g_irr, var='a')
print(f"Cuerpo: {F289}")

n_q = 8
k_q = 3

# Necesitamos una raíz de orden n_q en F_289*; el orden del grupo es 288 = 2^5 * 9
# Tomamos un generador aleatorio y elevamos a 288/n_q
import random
random.seed(0)
generador = None
for _ in range(100):
    cand = F289.aleatorio()
    if F289.es_cero(cand):
        continue
    if F289.pot(cand, 288 // 2) != F289.uno() and F289.pot(cand, 288 // 3) != F289.uno():
        generador = cand
        break

if generador is None:
    print("No se encontró generador, saltando ejemplo.")
else:
    print(f"Generador de F_289*: {F289.conv_a_str(generador)}")
    g_q = F289.pot(generador, 288 // n_q)
    print(f"Raíz de orden {n_q}: {F289.conv_a_str(g_q)}")
    print(f"g_q^{n_q} = {F289.conv_a_str(F289.pot(g_q, n_q))}")

    # Vector de elementos de F_q
    a_q = tuple(F289.elem_de_int(i) for i in range(n_q))
    print(f"a = {[F289.conv_a_str(x) for x in a_q]}")

    A_q = fq_fft(F289, g_q, k_q, a_q)
    print(f"FFT(a) = {[F289.conv_a_str(x) for x in A_q]}")

    a_rec_q = fq_ifft(F289, g_q, k_q, A_q)
    print(f"¿Recupera a? {all(F289.es_igual(x, y) for x, y in zip(a_q, a_rec_q))}")