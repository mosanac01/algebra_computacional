from .cuerpo_fp import cuerpo_fp
from .anillo_fp_x import anillo_fp_x

class cuerpo_fq:
    def __init__(self, fp, g, var='a'): # construye el cuerpo Fp[var]/<g(var)> donde g es objeto fabricado por fp
        if not isinstance(fp, cuerpo_fp):
            raise TypeError("fp debe ser una instancia de cuerpo_fp")
        
        # Crea anillo Fp[x]
        self.anillo = anillo_fp_x(fp, var)

        # Comprobar que g(x) es un polinomio válido
        if self.anillo.es_cero(g):
            raise ValueError("El polinomio g no puede ser el polinomio nulo")
        if not self.anillo.es_irreducible(g):
            raise ValueError("El polinomio g debe ser irreducible en Fp[x]")
        
        # Mas atributos
        self.fp = cuerpo_fp(fp.p)   # Copia segura
        self.p  = self.fp.p
        self.g  = self.anillo.hacer_monico(g)
        self.n  = self.anillo.grado(self.g)
        self.q  = self.p ** self.n
        self.var = var

    def __str__(self):
        return f"F_{self.p}^{self.n} = F_{self.p}[{self.var}]/<{self.anillo.conv_a_str(self.g)}>"


    # Elementos neutros

    def cero(self):                # 0
        return self.anillo.cero()

    def uno(self):                 # 1
        return self.anillo.uno()


    # Auxiliar

    def r_mod_g(self, a):   # Reduce el polinomio a modulo g(x) en Fp[x]
        return self.anillo.mod(a, self.g)

    
    # Las siguientes funciones fabrican elementos del cuerpo F_q

    def elem_de_tuple(self, a):    # fabrica elemento a partir de tupla de coeficientes
        return self.r_mod_g(self.anillo.elem_de_tuple(a))

    def elem_de_int(self, a):      # fabrica elemento a partir de entero
        return self.r_mod_g(self.anillo.elem_de_int(a))

    def elem_de_str(self, s):      # fabrica elemento parseando string
        return self.r_mod_g(self.anillo.elem_de_str(s))


    # Las siguientes funciones conviernten un elemento de F_q a otro formato

    def conv_a_tuple(self, a):     # devuelve tupla de coeficientes sin ceros "extra"
        # return self.anillo.conv_a_tuple(a) (esta es mas optima, pero la otra asegura 
        # que el polinomio esta reducido modulo g(x))
        return self.anillo.conv_a_tuple(self.r_mod_g(a))

    def conv_a_int(self, a):       # devuelve el entero correspondiente
        return self.anillo.conv_a_int(self.r_mod_g(a))

    def conv_a_str(self, a):       # pretty-printer
        return self.anillo.conv_a_str(self.r_mod_g(a))


    # Comparaciones

    def es_cero(self, a):          # a == 0
        return self.anillo.es_cero(self.r_mod_g(a))

    def es_uno(self, a):           # a == 1
        return self.anillo.es_uno(self.r_mod_g(a))

    def es_igual(self, a, b):      # a == b
        return self.anillo.es_igual(self.r_mod_g(a), self.r_mod_g(b))

    
    # Las siguientes operaciones asumen a y b son polinomios de F_q

    def suma(self, a, b):          # a+b
        return self.r_mod_g(self.anillo.suma(a, b))

    def inv_adit(self, a):         # -a
        return self.r_mod_g(self.anillo.inv_adit(a))

    def mult(self, a, b):          # a*b
        return self.r_mod_g(self.anillo.mult(a, b))

    def pot(self, a, k):           # a^k (k entero)
        return self.anillo.pot_mod(a, k, self.g)

    def inv_mult(self, a):
        if self.es_cero(a):
            raise ZeroDivisionError("No se puede invertir el elemento nulo")
        return self.anillo.inv_mod(a, self.g)
    
    
    # Aleatoriedad

    def aleatorio(self):           # devuelve un elemento aleatorio con prob uniforme
        return self.elem_de_tuple(tuple(self.fp.aleatorio() for _ in range(self.n)))
    

    # Tablas de operaciones

    def tabla_suma(self):          # matriz de qxq correspondiente a la suma (con la notación int)
        return [[self.conv_a_int(self.suma(self.elem_de_int(i), self.elem_de_int(j))) for j in range(self.q)] for i in range(self.q)]
    
    def tabla_mult(self):          # matriz de qxq correspondiente a la mult (con la notación int)
        return [[self.conv_a_int(self.mult(self.elem_de_int(i), self.elem_de_int(j))) for j in range(self.q)] for i in range(self.q)]

    def tabla_inv_adit(self):      # lista de inv_adit (con la notación int)
        return [self.conv_a_int(self.inv_adit(self.elem_de_int(j))) for j in range(self.q)]

    def tabla_inv_mult(self):      # lista de inv_mult (con la notación int)
        return ["*"] + [self.conv_a_int(self.inv_mult(self.elem_de_int(j))) for j in range(1, self.q)]


    # Cuadrado latino

    def cuadrado_latino(self, a):  # cuadrado latino para a != 0 (con notación int)
        a = self.r_mod_g(a)
        if self.es_cero(a):
            raise ValueError("a debe ser no nulo")
        return [[self.conv_a_int(self.suma(self.mult(a, self.elem_de_int(i)), self.elem_de_int(j))) for j in range(self.q)] for i in range(self.q)]
