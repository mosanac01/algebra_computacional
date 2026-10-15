from .cuerpo_fp import cuerpo_fp
from .utils import factorizar_entero

class anillo_fp_x:

    def __init__(self, fp, var='x'): # construye el anillo Fp[var]
        if not isinstance(fp, cuerpo_fp):
            raise TypeError("fp debe ser una instancia de cuerpo_fp")
        self.fp  = cuerpo_fp(fp.p)   # copia el cuerpo fp para evitar modificaciones indeseadas
        self.p   = self.fp.p         # toma p de la copia fp
        self.var = var               # nombre de la variable (string)

    def __str__(self):
        return f"F_{self.p}[{self.var}]"


    # Auxiliar
    
    def reducir(self, a):               # Quita ceros de cola, asegura tupla inmutable y coeficientes mod p.
        if not a: return (0,)
        a = [x % self.p for x in a]     # reduce coeficientes mod p
        while len(a) > 1 and a[-1] == 0:
            a.pop()
        return tuple(a)


    # Elementos neutros

    def cero(self):                 # 0
        return (0,)

    def uno(self):                  # 1
        return (1,)


    # Las siguientes funciones fabrican polinomios de F_p[x]

    def elem_de_tuple(self, a):     # fabrica un polinomio a partir de la tupla (a0, a1, ...)
        return self.reducir(a)
    
    def elem_de_int(self, a):       # fabrica un polinomio a partir de los dígitos de a en base p
        if not isinstance(a, int):
            raise TypeError("a debe ser un entero")
        if a < 0:
            raise ValueError("a debe ser un entero no negativo")
        if a == 0:
            return (0,)
        coef = []
        while a > 0:
            coef.append(a % self.p)
            a //= self.p            # divide entre p y actualiza a
        return self.reducir(coef)
    
    def elem_de_str(self, s):       # fabrica un polinomio a partir de un string (parser)
        """
        Toma un string donde los polinomios se escriben como sumas o restas de monomios del tipo "ax^b" 
        o "a*x^b" con a y b enteros es decir, ejemplos validos son "3x^2 + 2x + 1", "", "-x" "-8*x^3 + 3*x". 
        No admite expresiones como "2x3" que no tendrian sentido pues en todo caso faltaria ^ entre x y 3.
        Ni tampoco "3*4x^2", "3x*2", etc.
        """
        if not isinstance(s, str):          # asegura que s es un string
            raise TypeError("s debe ser un string")
        
        s = s.replace(' ', '')              # quitar espacios

        if '++' in s or '--' in s or '+-' in s or '-+' in s:
            raise ValueError("Signos consecutivos inválidos en la expresión")

        if s in ('0', ''):                  # polinomio nulo
            return (0,)
        
        if self.var not in s:               # caso sin variable, solo coeficiente
            try:
                s_aux = int(s) % self.p
                return (s_aux, )             # tupla de un elemento
            except ValueError:
                raise ValueError(f"No se pudo interpretar '{s}' como polinomio")
            
        # Caso general:
        # separar en sumandos
        if s[0] not in '+-':
            s = '+' + s    # asegura que el primer término tiene signo (si no tiene - pone +)

        sumandos = []
        i = 0
        while i < len(s):  # Escanea la cadena y corta en cada signo + o -
            j = i + 1
            while j < len(s) and s[j] not in '+-':
                j += 1
            sumandos.append(s[i:j])
            i = j

        # procesar cada sumando
        coef_dict = {}     # diccionario:  grado -> coeficiente

        for term in sumandos:
            if term == '+' or term == '-':
                raise ValueError("Término inválido en el polinomio")
            
            signo = 1 if term[0] == '+' else -1     # extraer signo
            term = term[1:]                         # quitar signo

            if '*' in term:
                # solo puede haber un '*' y debe estar justo antes de la variable
                if term.count('*') > 1:
                    raise ValueError(f"Uso invalido de '*' en el termino '{term}'")
                star_pos = term.index('*')
                var_pos = term.index(self.var)
                if star_pos != var_pos - 1:
                    raise ValueError(f"El '*' solo puede ir entre el coeficiente y la variable en '{term}'")
                term = term.replace('*', '')    # Elimina el '*'
            
            if self.var not in term:                # Caso solo coeficiente
                try:
                    c = int(term) * signo
                except ValueError:
                    raise ValueError(f"Coeficiente invalido en termino '{term}'")
                deg = 0

            else:                                   # Caso con variable
                var_pos = term.index(self.var)      # Posición donde empieza la variable
                caret_pos = term.find('^', var_pos + len(self.var)) # Buscar '^' despues de la variable, si no está -1

                if caret_pos != -1:     # Hay exponente
                    if var_pos > 0:
                        c_part = term[:var_pos]     # Parte del coeficiente
                        if c_part not in ('', '+', '-'):
                            try:
                                c = int(c_part) * signo # Coeficiente a número con signo (si es vacío o solo signo, vale 1 o -1)
                            except ValueError:
                                raise ValueError(f"Coeficiente invalido en termino '{term}'")
                        else:
                            c = signo
                        
                    else:
                        c = signo                     # Caso x^k sin coeficiente, vale 1 o -1
                    
                    deg_part = term[caret_pos+1:]     # Parte del exponente (después de ^)

                    if not deg_part:
                        raise ValueError(f"Exponente ausente en termino '{term}'")
                    if not deg_part.isdigit():
                        raise ValueError(f"Exponente invalido en termino '{term}'")

                    deg = int(deg_part)
                    
                else:                               # Caso x sin exponente
                    if var_pos > 0:
                        c_part = term[:var_pos]         # Parte del coeficiente
                        if c_part not in ('', '+', '-'):
                            try:
                                c = int(c_part) * signo # Coeficiente a número con signo
                            except ValueError:
                                raise ValueError(f"Coeficiente invalido en termino '{term}'")
                        else:
                            c = signo

                    else:
                        c = signo
                    deg = 1                         # Exponente vale 1
            
            # En este punto deg siempre es entero >= 0 (ya validado anteriormente)
            c = c % self.p          # Coeficiente mod p
            if deg in coef_dict:    # Acumular coeficientes de mismo grado y sumar mod p
                coef_dict[deg] = (coef_dict[deg] + c) % self.p
            else:
                coef_dict[deg] = c

        if not coef_dict:           # si quedó vacío (todos coef 0) -> polinomio nulo
            return (0,)
        
        max_deg = max(coef_dict.keys())     # grado máximo // coef_dict.keys() da una lista de todos lso grados
        coef = [0]*(max_deg + 1)            # lista de coeficientes inicializada a 0
        
        for d, c in coef_dict.items():      # coef_dict.items() devuelve pares (grado, coeficiente)
            coef[d] = c                     # coloca el coeficiente c en la posición d de la lista coef

        return self.reducir(coef)   


    # Las siguientes funciones convierten un polinomio de F_p[x] (asumo que ya viene como tupla o si no como lista) a otro formato    

    def conv_a_tuple(self, a):      # devuelve la tupla de coeficientes
        return self.reducir(a)

    def conv_a_int(self, a):        # devuelve el entero correspondiente al polinomio
        a = self.reducir(a)
        n = 0
        for i in range(len(a)-1, -1, -1):   # Recorre los coeficientes de mayor a menor grado
            n = n * self.p + a[i]           # suma coeficientes por p^i
        return n
    
    def conv_a_str2(self, a):         # pretty-printer
        a = self.reducir(a)
        terms = []
        for i in range(len(a)-1, -1, -1):  # Recorre los coeficientes de mayor a menor grado
            if a[i] == 0:
                continue
            if i == 0:
                terms.append(str(a[i]))
            elif i == 1:
                if a[i] == 1:
                    terms.append(self.var)
                else:
                    terms.append(f"{a[i]}{self.var}")
            else:
                if a[i] == 1:
                    terms.append(f"{self.var}^{i}")
                else:
                    terms.append(f"{a[i]}{self.var}^{i}")
        if not terms:
            return '0'
        return '+'.join(terms)
    
    def conv_a_str(self, a):         # pretty-printer 
        a = self.reducir(a)
        terms = []
        for i in range(len(a)-1, -1, -1):  # Recorre los coeficientes de mayor a menor grado
            if a[i] == 0:
                continue
            if i == 0:      # Termino independiente
                terms.append(str(a[i]))
            elif i == 1:    # Grado 1
                if a[i] == 1:
                    terms.append(self.var)
                else:
                    terms.append(f"{a[i]}*{self.var}")
            else:           # A partir de grado 2
                if a[i] == 1:
                    terms.append(f"{self.var}^{i}")
                else:
                    terms.append(f"{a[i]}*{self.var}^{i}")
        if not terms:
            return '0'
        return '+'.join(terms)


    # Comparaciones

    def es_cero(self, a):          # a == 0
        return self.reducir(a) == (0,)
    
    def es_uno(self, a):           # a == 1
        return self.reducir(a) == (1,)
    
    def es_igual(self, a, b):      # a == b
        return self.reducir(a) == self.reducir(b)


    # Las siguientes operaciones asumen que a y b son elementos de Fp[x]

    def suma(self, a, b):           # a+b
        a = self.reducir(a)
        b = self.reducir(b)
        na, nb = len(a), len(b)
        if na < nb:
            na, nb, a, b = nb, na, b, a
        c = tuple(self.fp.suma(a[i], b[i]) for i in range(nb)) + a[nb:]
        return self.reducir(c)
    
    def inv_adit(self, a):          # -a
        a = self.reducir(a)
        return self.reducir(tuple(self.fp.inv_adit(x) for x in a))
    
    def mult(self, a, b):           # a*b
        a = self.reducir(a)
        b = self.reducir(b)
        na, nb = len(a), len(b)
        if self.es_cero(a) or self.es_cero(b):
            return self.cero()
        c = [self.fp.cero()] * (na + nb - 1)
        for i in range(na):
            for j in range(nb):
                c[i + j] = self.fp.suma(c[i + j], self.fp.mult(a[i], b[j]))
        return self.reducir(c)

    def mult_por_escalar(self, a, e): # a*e (con e en Z/pZ)
        a = self.reducir(a)
        e = self.fp.elem_de_int(e)
        return self.reducir(tuple(self.fp.mult(x, e) for x in a))   

    def divmod(self, a, b):        # devuelve q,r tales que a=bq+r y deg(r)<deg(b)
        a = self.reducir(a)
        b = self.reducir(b)
        na, nb = len(a), len(b)

        if self.es_cero(b):
            raise ZeroDivisionError("No se puede dividir por el polinomio nulo")
        if na < nb:
            return self.cero(), a
        
        a = list(a)                # convertir para permitir mutabilidad
        b_lead_inv = self.fp.inv_mult(b[nb - 1])    # Inverso del coeficiente lider de b
        q = [self.fp.cero()]*(na - nb + 1)          # Coeficientes de q inicializado a cero
        
        for k in range(na - nb, -1, -1):
            q[k] = self.fp.mult(a[k + nb - 1], b_lead_inv)
            for j in range(nb):
                a[k + j] = self.fp.suma(a[k+j], self.fp.inv_adit(self.fp.mult(q[k], b[j])))

        r = self.reducir(a)                         # resto de lo que queda de a
        return self.reducir(q), r

    def div(self, a, b):           # q
        return self.divmod(a, b)[0]
    
    def mod(self, a, b):           # r
        return self.divmod(a, b)[1]
    
    def grado(self, a):            # deg(a)
        a = self.reducir(a)
        if self.es_cero(a):
            return -1
        return len(a) - 1
    

    # Auxiliar

    def hacer_monico(self, a): # Convierte un polinomio no nulo a su forma reducida monico
        a = self.reducir(a)

        if self.es_cero(a):
            raise ValueError("No se puede hacer monico el polinomio cero")
        
        lead_coef = a[-1]
        if lead_coef == 1:
            return a
        
        lead_inv = self.fp.inv_mult(lead_coef)
        return self.mult_por_escalar(a, lead_inv)

        
    # GCD

    def gcd(self, a, b):           # devuelve g = gcd(a,b) mónico
        a = self.reducir(a)
        b = self.reducir(b)
        if self.es_cero(a) and self.es_cero(b):
            raise ValueError("gcd(0, 0) no esta definido")
        
        if self.es_cero(a):
            return self.hacer_monico(b)      # Devuelve b mónico
        if self.es_cero(b):
            return self.hacer_monico(a)      # Devuelve a mónico
        
        return self.gcd(b, self.mod(a, b))   # Recursión y Algoritmo de Euclides 

    def gcd_ext(self, a, b):       # devuelve g,x,y tales que g=ax+by, g=gcd(a,b) mónico
        a = self.reducir(a)
        b = self.reducir(b)
        if self.es_cero(a) and self.es_cero(b):
            raise ValueError("gcd(0, 0) no esta definido")
        
        if self.es_cero(a):
            y = self.fp.inv_mult(b[-1])
            return self.mult_por_escalar(b, y), self.cero(), (y,)
        if self.es_cero(b):
            x = self.fp.inv_mult(a[-1])
            return self.mult_por_escalar(a, x), (x,), self.cero()
        
        q, r = self.divmod(a, b)              # r = a - bq
        g, x1, y1 = self.gcd_ext(b, r)        # g = bx1 + ry1 = bx1 + (a - bq)y1 = ay1 + b(x1 - qy1)
        x = y1
        y = self.suma(x1, self.inv_adit(self.mult(q, y1)))

        return g, x, y
    

    # Inverso y Potencia Modulo

    def inv_mod(self, a, b):       # devuelve x tal que ax = 1 mod b
        a = self.reducir(a)
        b = self.reducir(b)

        if self.es_cero(b):
            raise ZeroDivisionError("El modulo b no puede ser el polinomio nulo")
        #if self.es_cero(a):            // 0x \equiv 1 (mod 1) si tiene solución, cualquier entero es solución
        #   raise ZeroDivisionError("El polinomio a = 0 no tiene inverso modulo b")

        g, x, y = self.gcd_ext(a, b)    # g = a*x + b*y

        if not self.es_uno(g):          # Si gcd != 1 el inverso no existe
            raise ValueError("Polinomio no invertible modulo b")
        
        return self.mod(x, b)      # Forma canonica de x en Fp[x]/(b)

    def pot_mod(self, a, k, b):    # a^k mod b
        a = self.reducir(a)
        b = self.reducir(b)

        if self.es_cero(b):
            raise ZeroDivisionError("No se puede calcular potencia modulo el polinomio nulo")
        
        if k == 0:
            return self.uno()
        
        if k < 0:
            a = self.inv_mod(a, b)
            k = -k
        
        resultado = self.uno()
        base = self.mod(a, b)

        while k > 0:
            if k % 2 == 1:
                resultado = self.mod(self.mult(resultado, base), b)
            base = self.mod(self.mult(base, base), b)
            k //= 2

        return self.reducir(resultado)
    

    # Irreducible

    def es_irreducible(self, f):   # test de irreducibilidad de Rabin
        f = self.reducir(f)
        if self.es_cero(f):
            raise ValueError("La irreducibilidad del polinomo nulo no está definida")
        
        n = self.grado(f)
        if n <= 1:          # Grado 0 o 1
            return True
        
        p = self.p

        # Comprobar que x^(p^n) = x (mod f)
        x = (0, 1)                      # el polinomio x
        g = self.pot_mod(x, p**n, f)    # x^(p^n) (mod f)
        if not self.es_igual(g, x):
            return False
        
        # Para cada divisor primo d de n comprobar gcd(f, x^(p^(n/d)) - x) = 1
        factores = factorizar_entero(n)
        for(d, _) in factores:
            exp = p**(n//d)
            g = self.pot_mod(x, exp, f)
            g = self.suma(g, self.inv_adit(x))      # g = x^(p^(n/d)) - x
            h = self.gcd(f, g)
            if not self.es_uno(h):
                return False
            
        # En caso contrario (pasa las condiciones anteriores) es irreducible
        return True
