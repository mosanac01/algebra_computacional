# pyalgcomputacional

Librería de álgebra computacional en Python que implementa:

- Cuerpos finitos `F_p` y extensiones `F_q = F_p[a]/⟨g(a)⟩`.
- Anillos de polinomios `F_p[x]` y `F_q[x]`.
- División polinómica rápida mediante **matrices de Toeplitz**.
- Multiplicación **Karatsuba**.
- Factorización **Cantor–Zassenhaus** en `F_p[x]` y `F_q[x]`.
- **FFT** (Cooley–Tukey) sobre `F_p` y `F_q`.

## Instalación

```bash
pip install -r requirements.txt