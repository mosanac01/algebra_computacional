# pyalgcomputacional 🧮

Una librería modular de álgebra computacional en Python. El paquete implementa estructuras algebraicas desde cero e incorpora implementaciones de algoritmos de multiplicación y división polinómica rápida, así como transformadas de Fourier sobre cuerpos finitos, con el objetivo de mejorar la eficiencia de determinadas operaciones algebraicas, evitando depender de librerías externas para la lógica matemática central.

Este proyecto se desarrolló individualmente para la asignatura de Álgebra Computacional del Grado en Matemáticas.

El profesor proporcionó los esqueletos iniciales, que definían las clases y las interfaces de las funciones. A partir de ellos y de los contenidos teóricos de la asignatura, implementé las estructuras algebraicas y los algoritmos, tomando decisiones sobre la representación interna de los elementos y su implementación.

Posteriormente, reorganicé el código original en un paquete modular para facilitar su mantenimiento, comprensión y reutilización.


## 🚀 Características Principales

* **Estructuras Algebraicas Básicas:**
  * Cuerpos finitos $F_p$ (primos) y extensiones $F_q = F_p[\alpha]/\langle g(\alpha)\rangle$.
  * Anillos de polinomios $F_p[x]$ y $F_q[x]$.
* **Algoritmos Rápidos (Divide y Vencerás):**
  * **Multiplicación polinómica:** Algoritmo de **Karatsuba**.
  * **División rápida polinómica:** Inversión mediante **matrices de Toeplitz**.
  * **Transformadas:** FFT e IFFT utilizando el algoritmo de **Cooley-Tukey** sobre $F_p$ y $F_q$.
* **Factorización Avanzada:**
  * Factorización completa de polinomios usando el algoritmo de **Cantor-Zassenhaus** (Square-free, Distinct-degree y Equal-degree factorization).


## 📁 Arquitectura del Proyecto

```text
algebra_computacional/
├── pyalgcomputacional/
│   ├── __init__.py         # Inicializador del paquete Python
│   ├── anillo_fp_x.py      # Lógica de polinomios sobre Fp
│   ├── anillo_fq_x.py      # Lógica de polinomios sobre extensiones Fq
│   ├── cuerpo_fp.py        # Cuerpos finitos base
│   ├── cuerpo_fq.py        # Extensiones de cuerpos finitos
│   ├── division.py         # División rápida (Toeplitz)
│   ├── factorizacion.py    # Algoritmo Cantor-Zassenhaus (Puro)
│   ├── karatsuba.py        # Multiplicación rápida (Puro)
│   ├── toeplitz.py         # Álgebra lineal matricial sobre cuerpos finitos
│   └── transformadas.py    # FFT / Cooley-Tukey (Puro)
├── ejemplos/
│   ├── 01_CuerposFinitos.py
│   ├── 02_AnillosPolinomios.py
│   ├── 03_DivisionToeplitzKaratsuba.py
│   ├── 04_FactorizacionFpx.py
│   ├── 05_FactorizacionFqx.py
│   └── 06_TransformadasFFT.py
├── .gitignore              # Reglas de exclusión para Git
├── requirements.txt
└── README.md
```

## ⚙️ Instalación y Requisitos

El proyecto requiere **Python 3.10+**.

1. Clona el repositorio:
```bash
git clone https://github.com/mosanac01/algebra_computacional.git
cd algebra_computacional
```

2. Instala las dependencias (requiere `sympy>=1.12` exclusivamente para tests de primalidad iniciales):
```bash
pip install -r requirements.txt
```

## 💻 Ejemplo de Uso (Quick Start)

Ejemplo de cómo multiplicar y factorizar polinomios usando la librería:

```python
from pyalgcomputacional.cuerpo_fp import cuerpo_fp
from pyalgcomputacional.anillo_fp_x import anillo_fp_x
from pyalgcomputacional.factorizacion import fact_fpx

# 1. Definir el cuerpo F_11 y el anillo de polinomios F_11[x]
fp = cuerpo_fp(11)
fpx = anillo_fp_x(fp, var='x')

# 2. Fabricar polinomios
f = fpx.elem_de_str("x^2 + 5*x + 6")
print(f"Polinomio f(x) = {fpx.conv_a_str(f)}")

# 3. Factorizar usando Cantor-Zassenhaus
factores = fact_fpx(fpx, f)
print("Factores irreducibles:")
for base, exponente in factores:
    print(f" -> ({fpx.conv_a_str(base)})^{exponente}")
```

## 🧪 Ejemplos

El directorio `ejemplos/` contiene seis scripts que demuestran el funcionamiento de las estructuras y algoritmos.
