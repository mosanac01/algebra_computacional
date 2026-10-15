import random
from sympy import isprime

class cuerpo_fp:

    def __init__(self, p): # construye el cuerpo de p elementos Fp = Z/pZ
        if not isinstance(p, int):
            raise TypeError("p debe ser un entero")
        if p <= 1 or not isprime(p):
            raise ValueError("p debe ser un número primo")
        self.p = p

    def __str__(self):
        return f"F_{self.p}"
    
    
    # Elementos neutros
        
    def cero(self):        # devuelve el elemento 0
        return 0
    
    def uno(self):         # devuelve el elemento 1
        return 1
    

    # Las siguientes funciones fabrican elementos del cuerpo F_p
    
    def elem_de_int(self, n): # fabrica el elemento dado por la clase de n
        if not isinstance(n, int):
            raise TypeError("n debe ser un entero")
        return n % self.p
    
    def elem_de_str(self, s): # fabrica el elemento a partir de un string (parser)
        if not isinstance(s, str):
            raise TypeError("s debe ser un string")
        try:
            return int(s) % self.p
        except ValueError:
            raise ValueError("s debe representar un entero válido")


    # Las siguientes funciones convierten un elemento del cuerpo F_p a otro formato
    # Asumo que para las siguientes funciones los valores que acepta son elementos de F_p
    # es decir, ya son enteros entre 0 y p-1 (aun asi hare una comprobacion minima)

    def conv_a_int(self, a):  # devuelve un entero entre 0 y p-1
        return int(a) % self.p
    
    def conv_a_str(self, a):  # pretty-printer
        return str(int(a) % self.p)


    # Las siguientes operaciones asumen que a y b son enteros

    def suma(self, a, b):     # a+b
        return (a + b) % self.p
    
    def inv_adit(self, a):    # -a
        return (-a) % self.p
    
    def mult(self, a, b):     # a*b
        return (a * b) % self.p

    def inv_mult(self, a):
        a = a % self.p
        if a == 0:
            raise ZeroDivisionError("No se puede invertir el elemento cero")
        return pow(a, -1, self.p)

    def pot(self, a, k):      # a^k (k entero)
        if k < 0:
            a = self.inv_mult(a) # lo reemplazo por su inverso
            k = -k
        # Válido para cualquier entero k (incluso k=0)
        return pow(a, k, self.p)


    # Comparaciones

    def es_cero(self, a):     # a == 0
        return (a % self.p) == 0
    
    def es_uno(self, a):      # a == 1
        return (a % self.p) == 1
    
    def es_igual(self, a, b): # a == b
        return (a % self.p) == (b % self.p)


    # Aleatoriedad

    def aleatorio(self):      # fabrica un elemento aleatorio con prob uniforme
        return random.randint(0, self.p - 1)


    # Tablas de operaciones

    def tabla_suma(self):     # devuelve la matriz de pxp (lista de listas) de la suma
        return [[self.suma(i, j) for j in range(self.p)] for i in range(self.p)]

    def tabla_mult(self):     # devuelve la matriz de pxp (lista de listas) de la mult
        return [[self.mult(i, j) for j in range(self.p)] for i in range(self.p)]    

    def tabla_inv_adit(self): # devuelve una lista de long p con los inv_adit
        return [self.inv_adit(i) for i in range(self.p)]

    def tabla_inv_mult(self): # devuelve una lista de long p con los inv_mult (en el índice 0 pone un '*')
        return ['*'] + [self.inv_mult(i) for i in range(1, self.p)]

    
    # Cuadrado latino

    def cuadrado_latino(self, a): # cuadrado latino a*i+j (con a != 0)
        a = a % self.p
        if a == 0:
            raise ValueError("a debe ser distinto de 0")
        return [[self.suma(self.mult(a, i), j) for j in range(self.p)] for i in range(self.p)]
