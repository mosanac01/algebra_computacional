from .cuerpo_fp import cuerpo_fp
from .cuerpo_fq import cuerpo_fq
from .utils import factorizar_entero

class anillo_fq_x:
    def __init__(self, fq, var='x'): # Fq[var], var debe ser distinta que la de fq
        if not isinstance(fq, cuerpo_fq):
            raise TypeError("fq debe ser una instancia de cuerpo_fq")
        if var == fq.var:
            raise ValueError("La variable del anillo debe ser distinta a la del cuerpo base")
        
        self.fq   = cuerpo_fq(fq.fp, fq.g, fq.var) # Copia segura F_q
        self.fp   = cuerpo_fp(fq.p)                # Cuerpo base
        self.g    = self.fq.g                           # Polinomio irreducible usado para generar F_q
        
        self.var    = var

        self.q = self.fq.q                 # Cardinal del cuerpo F_q
        self.p = self.fq.p                 # Caracteristica
        self.n = self.fq.n                 # Grado de extension del cuerpo base
    
    def __str__(self):
        return f"F_{self.q}[{self.var}]"
    

    # Auxiliar Reducir

    def reducir(self, a):
        if not a:
            return (self.fq.cero(),)
        
        try:
            a = tuple(a)
        except TypeError:
            raise TypeError("a debe ser lista o tupla (iterable)")
        
        a = tuple(self.fq.r_mod_g(ai) for ai in a)  # reduce coeficientes modulo g(x)

        # Elimina ceros al final
        i = len(a) - 1
        while i > 0 and a[i] == self.fq.cero():
            i -= 1
        return a[:i+1]


    # Elementos neutros

    def cero(self):                # 0 = ((0,),)
        return (self.fq.cero(),)

    def uno(self):                 # 1 = ((1,),)
        return (self.fq.uno(),)


    # Las siguientes funciones fabrican polinomios de F_q[x]:

    def elem_de_tuple(self, a):
        return self.reducir(a)

    def elem_de_int(self, a):
        """
        Fabrica un elemento del anillo F_q[x] a partir de un entero no negativo.
        Interpreta los 'dígitos' de a en base q como coeficientes en F_q.
        Ejemplo: si q = 25 y a = 51, entonces 51 = 1·25^0 + 2·25^1  ->  (1, 2)  ->  f(x) = 1 + 2x
        """
        if not isinstance(a, int):
            raise TypeError("a debe ser un entero")
        if a < 0:
            raise ValueError("a debe ser un entero no negativo")
        if a == 0:
            return self.cero()
        
        coef = []
        while a > 0:
            digito = a % self.q                         # resto en base q
            coef.append(self.fq.elem_de_int(digito))    # convertir a elemento de F_q
            a //= self.q                                # cociente entero

        while len(coef) > 1 and coef[-1] == self.fq.cero():     # eliminar ceros finales
            coef.pop()
        return tuple(coef)

    def elem_de_str(self, s):
        """
        Fabrica un polinomio de F_q[x] a partir de un string, análogo al caso F_p[x].
        Los coeficientes ahora se interpretan en F_q, y pueden ser expresiones
        en la variable del cuerpo base (por ejemplo, 'a+1', '2a', '3a^2+4').
        """
        if not isinstance(s, str):
            raise TypeError("s debe ser un string")

        s = s.replace(' ', '')

        if '++' in s or '--' in s or '+-' in s or '-+' in s:
            raise ValueError("Signos consecutivos inválidos en la expresión")

        if s in ('0', ''):
            return self.cero()

        if self.var not in s:       # No hay variable x, solo un coeficiente en F_q
            try:
                return (self.fq.elem_de_str(s),)
            except ValueError:
                raise ValueError(f"No se pudo interpretar '{s}' como un polinomio de Fq[x]")
            
        # Caso general
        # Separar los sumandos
        # Asegura que el primer término tenga signo explícito
        if s[0] not in '+-':
            s = '+' + s

        sumandos = []
        nivel_par = 0
        inicio = 0

        for i in range(1, len(s)):
            c = s[i]
            if c == '(':
                nivel_par += 1
            elif c == ')':
                nivel_par -= 1
                if nivel_par < 0:
                    raise ValueError("Parentesis desbalanceados")
            elif c in '+-' and nivel_par == 0:
                sumandos.append(s[inicio:i])
                inicio = i

        sumandos.append(s[inicio:]) # ultimo termino

        coef_dict = {}

        for term in sumandos:
            if term in ('+', '-'):
                raise ValueError("Término inválido en el polinomio")

            signo = 1 if term[0] == '+' else -1
            term = term[1:]

            if self.var not in term:        # Solo coeficiente
                coef_str = term
                coef = self.fq.elem_de_str(coef_str)
                if signo == -1:
                    coef = self.fq.inv_adit(coef)
                deg = 0
                if deg in coef_dict:
                    coef_dict[deg] = self.fq.suma(coef_dict[deg], coef)
                else:
                    coef_dict[deg] = coef
                continue    # Pasa al siguiente term del for

            if '*' in term:
                # Puede haber varios '*' dentro del coeficiente,
                # pero solo uno fuera de paréntesis y justo antes de la variable principal (x).
                nivel_par = 0
                star_pos_externo = -1
                for i, c in enumerate(term):
                    if c == '(':
                        nivel_par += 1
                    elif c == ')':
                        nivel_par -= 1
                        if nivel_par < 0:
                            raise ValueError("Paréntesis desbalanceados")
                    elif c == '*' and nivel_par == 0:
                        # '*' fuera de paréntesis
                        if star_pos_externo != -1:
                            raise ValueError(f"Uso inválido de '*' en el término '{term}' (más de un '*' fuera de paréntesis)")
                        star_pos_externo = i

                # Si hay uno fuera de paréntesis, debe estar justo antes de la variable principal
                if star_pos_externo != -1:
                    var_pos = term.find(self.var, star_pos_externo + 1)
                    if star_pos_externo != var_pos - 1:
                        raise ValueError(f"El '*' externo solo puede ir entre el coeficiente y la variable principal en '{term}'")

                # --- Detectar el coeficiente antes del '*' externo ---
                coef = term[:star_pos_externo]

                # --- Detectar variable y exponente ---
                var_part = term[var_pos:]

                # Ahora tratamos el coeficiente, si tiene paréntesis los eliminamos
                if coef.startswith('('):
                    coef = coef[1:-1]
                
                # Convertir el coeficiente a elemento de F_q
                coef = self.fq.elem_de_str(coef)
                
                # Aplicar signo
                if signo == -1:
                    coef = self.fq.inv_adit(coef)

                # Ahora trabajamos con la parte de la variable
                exp_pos = var_part.find('^')
                if exp_pos == -1:
                    deg = 1
                else:
                    exp_str = var_part[exp_pos + 1:]
                    if not exp_str.isdigit():
                        raise ValueError(f"Exponente inválido en termino '{term}'")
                    deg = int(exp_str)

                # --- Guardar en el diccionario ---
                if deg in coef_dict:
                    coef_dict[deg] = self.fq.suma(coef_dict[deg], coef)
                else:
                    coef_dict[deg] = coef
                continue    # Pasa al siguiente term del for

            if term == self.var or term.startswith(f"{self.var}^"):     # Solo variable x^j con coeficiente 1 implicito
                coef = self.fq.uno()
                if signo == -1:
                    coef = self.fq.inv_adit(coef)

                # Extraer el exponente si hay "^", si no, es 1
                caret_pos = term.find('^')
                if caret_pos != -1:
                    exp_str = term[caret_pos + 1:]
                    if not exp_str.isdigit():
                        raise ValueError(f"Exponente inválido en término '{term}'")
                    deg = int(exp_str)
                else:
                    deg = 1

                # Guardar en el diccionario
                if deg in coef_dict:
                    coef_dict[deg] = self.fq.suma(coef_dict[deg], coef)
                else:
                    coef_dict[deg] = coef
                continue

            if not term.startswith('('):        # Caso z a^j x^i con z un entero y a la variable de Fq que trae de Fp[a]
                var_pos = term.index(self.var)
                coef_str = term[:var_pos]

                coef = self.fq.elem_de_str(coef_str)
                if signo == -1:
                    coef = self.fq.inv_adit(coef)

                # Determinar exponente
                caret_pos = term.find('^', var_pos + len(self.var))
                if caret_pos != -1:
                    exp_str = term[caret_pos + 1:]
                    if not exp_str.isdigit():
                        raise ValueError(f"Exponente inválido en término '{term}'")
                    deg = int(exp_str)
                else:
                    deg = 1

                # Guardar
                if deg in coef_dict:
                    coef_dict[deg] = self.fq.suma(coef_dict[deg], coef)
                else:
                    coef_dict[deg] = coef
                continue

            if term.startswith('('):        # Caso (coef) x^i
                close_par = term.find(')')
                if close_par == -1:
                    raise ValueError(f"Coeficiente mal formado en '{term}' (falta cierre de paréntesis)")

                coef_str = term[1:close_par]
                coef = self.fq.elem_de_str(coef_str)
                if signo == -1:
                    coef = self.fq.inv_adit(coef)

                # Determinar exponente
                caret_pos = term.find('^', close_par + 1)
                if caret_pos != -1:
                    exp_str = term[caret_pos + 1:]
                    if not exp_str.isdigit():
                        raise ValueError(f"Exponente inválido en término '{term}'")
                    deg = int(exp_str)
                else:
                    deg = 1

                # Guardar
                if deg in coef_dict:
                    coef_dict[deg] = self.fq.suma(coef_dict[deg], coef)
                else:
                    coef_dict[deg] = coef
                continue

        if not coef_dict:
            return self.cero()

        max_deg = max(coef_dict.keys())
        coef = [self.fq.cero()] * (max_deg + 1)
        for d, c in coef_dict.items():
            coef[d] = c

        return self.reducir(coef)


    # Las siguientes funciones convierten un polinomio de F_q[x]

    def conv_a_tuple(self, a):
        return self.reducir(a)

    def conv_a_int(self, a):
        a = self.reducir(a)
        n = 0
        for i in range(len(a) - 1, -1, -1):
            coef_int = self.fq.conv_a_int(a[i])
            n = n * self.q + coef_int
        return n

    def conv_a_str(self, a):
        a = self.reducir(a)
        terms = []
        for i in range(len(a) - 1, -1, -1):
            if self.fq.es_cero(a[i]):
                continue

            coef_str = self.fq.conv_a_str(a[i])

            if i == 0:
                # término independiente
                term = coef_str

            elif i == 1:
                # término de grado 1
                if self.fq.es_uno(a[i]):
                    term = self.var
                else:
                    term = f"({coef_str})*{self.var}"

            else:
                # término de grado >= 2
                if self.fq.es_uno(a[i]):
                    term = f"{self.var}^{i}"
                else:
                    term = f"({coef_str})*{self.var}^{i}"

            terms.append(term)

        if not terms:
            return "0"
        return " + ".join(terms)


    # Comparaciones

    def es_cero(self, a):
        a = self.reducir(a)
        return len(a) == 1 and self.fq.es_cero(a[0])

    def es_uno(self, a):
        a = self.reducir(a)
        return len(a) == 1 and self.fq.es_uno(a[0])

    def es_igual(self, a, b):
        a = self.reducir(a)
        b = self.reducir(b)

        if len(a) != len(b):
            return False

        for ai, bi in zip(a, b):
            if not self.fq.es_igual(ai, bi):
                return False
        return True


    # Las siguientes operaciones asumen que a y b son elementos de Fq[x]

    def suma(self, a, b):
        a = self.reducir(a)
        b = self.reducir(b)
        na, nb = len(a), len(b)

        # aseguramos que 'a' sea el más largo
        if na < nb:
            na, nb, a, b = nb, na, b, a

        c = tuple(self.fq.suma(a[i], b[i]) for i in range(nb)) + a[nb:]
        return self.reducir(c)

    def inv_adit(self, a):
        a = self.reducir(a)
        return self.reducir(tuple(self.fq.inv_adit(x) for x in a))

    def mult(self, a, b):
        a = self.reducir(a)
        b = self.reducir(b)
        na, nb = len(a), len(b)

        if self.es_cero(a) or self.es_cero(b):
            return self.cero()

        # lista de coeficientes inicializada con ceros (de F_q)
        c = [self.fq.cero()] * (na + nb - 1)

        for i in range(na):
            for j in range(nb):
                prod = self.fq.mult(a[i], b[j])
                c[i + j] = self.fq.suma(c[i + j], prod)

        return self.reducir(c)

    def mult_por_escalar(self, a, e):
        a = self.reducir(a)
        e = self.fq.r_mod_g(e)  # asegurar que e sea elemento válido de F_q

        if self.fq.es_cero(e) or self.es_cero(a):
            return self.cero()

        return self.reducir(tuple(self.fq.mult(x, e) for x in a))

    def divmod(self, a, b):         # devuelve (q, r) tales que  a = b*q + r  y  grado(r) < grado(b)
        a = self.reducir(a)
        b = self.reducir(b)

        if self.es_cero(b):
            raise ZeroDivisionError("No se puede dividir por el polinomio nulo")

        na, nb = len(a), len(b)
        if na < nb:
            return self.cero(), a

        # convertir a lista para manipulación mutable
        a = list(a)
        # inverso del coeficiente líder de b en F_q
        b_lead_inv = self.fq.inv_mult(b[-1])

        q = [self.fq.cero()] * (na - nb + 1)

        # algoritmo clásico de división polinómica
        for k in range(na - nb, -1, -1):
            q[k] = self.fq.mult(a[k + nb - 1], b_lead_inv)
            for j in range(nb):
                a[k + j] = self.fq.suma(a[k + j], self.fq.inv_adit(self.fq.mult(q[k], b[j])))

        r = self.reducir(a)
        return self.reducir(q), r

    def div(self, a, b):
        return self.divmod(a, b)[0]

    def mod(self, a, b):
        return self.divmod(a, b)[1]

    def grado(self, a):
        a = self.reducir(a)
        if self.es_cero(a):
            return -1
        return len(a) - 1


    # GCD

    def hacer_monico(self, a):
        a = self.reducir(a)

        if self.es_cero(a):
            raise ValueError("No se puede hacer monico el polinomio cero")
        
        lead_coef = a[-1]
        if self.fq.es_uno(lead_coef):
            return a
        
        lead_inv = self.fq.inv_mult(lead_coef)
        return self.mult_por_escalar(a, lead_inv)

    def gcd(self, a, b):    # Devuelve g = gcd(a, b) mónico, calculado mediante el algoritmo de Euclides.
        a, b = self.reducir(a), self.reducir(b)

        if self.es_cero(a) and self.es_cero(b):
            raise ValueError("gcd(0, 0) no está definido")

        while not self.es_cero(b):
            _, r = self.divmod(a, b)
            a, b = b, r

        return self.hacer_monico(a)

    def gcd_ext(self, a, b):        # devuelve g,x,y tales que g=ax+by, g=gcd(a,b) mónico
        a = self.reducir(a)
        b = self.reducir(b)

        if self.es_cero(a) and self.es_cero(b):
            raise ValueError("gcd(0, 0) no está definido")

        # Casos base
        if self.es_cero(a):
            y = self.fq.inv_mult(b[-1])
            return self.mult_por_escalar(b, y), self.cero(), (y,)
        if self.es_cero(b):
            x = self.fq.inv_mult(a[-1])
            return self.mult_por_escalar(a, x), (x,), self.cero()

        # Paso recursivo
        q, r = self.divmod(a, b)
        g, x1, y1 = self.gcd_ext(b, r)

        # g = b*x1 + r*y1 = b*x1 + (a - b*q)*y1 = a*y1 + b*(x1 - q*y1)
        x = y1
        y = self.suma(x1, self.inv_adit(self.mult(q, y1)))

        return g, x, y


    # Inverso y Potencia Modulo

    def inv_mod(self, a, b):    #  Devuelve x tal que a·x ≡ 1 (mod b)
        a = self.reducir(a)
        b = self.reducir(b)

        if self.es_cero(b):
            raise ZeroDivisionError("El módulo b no puede ser el polinomio nulo")
        #if self.es_cero(a):            ax \equiv 1 (mod b) si tiene solucion
        #    raise ZeroDivisionError("El polinomio a = 0 no tiene inverso módulo b")

        g, x, y = self.gcd_ext(a, b)  # g = a·x + b·y

        if not self.es_uno(g):
            raise ValueError("Polinomio no invertible módulo b (gcd ≠ 1)")

        return self.mod(x, b)  

    def pot_mod(self, a, k, b):     # a^k mod b en F_q[x]
        a = self.reducir(a)
        b = self.reducir(b)

        if self.es_cero(b):
            raise ZeroDivisionError("No se puede calcular potencia módulo el polinomio nulo")

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

    def es_irreducible(self, f):        # test de irreducibilidad de Rabin
        f = self.reducir(f)
        if self.es_cero(f):
            raise ValueError("La irreducibilidad del polinomio nulo no está definida")

        n = self.grado(f)
        if n <= 1:
            return True  # grado 0 o 1: siempre irreducible

        q = self.q
        x = (self.fq.cero(), self.fq.uno())  # polinomio x

        # Comprobar que x^(q^n) ≡ x (mod f)
        g = self.pot_mod(x, q**n, f)
        if not self.es_igual(g, x):
            return False

        # Para cada divisor primo d de n, comprobar gcd(f, x^(q^(n/d)) - x) = 1
        factores = factorizar_entero(n)
        for (p, _) in factores:
            exp = q**(n // p)
            g = self.pot_mod(x, exp, f)
            g = self.suma(g, self.inv_adit(x))  # g = x^(q^(n/d)) - x
            h = self.gcd(f, g)
            if not self.es_uno(h):
                return False

        # Si pasa todas las pruebas → irreducible
        return True