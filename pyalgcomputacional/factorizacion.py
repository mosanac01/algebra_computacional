# ==========================================
# FACTORIZACIÓN CANTOR-ZASSENHAUS Fp[x]
# ==========================================

# square-free factorization
# input: fpx --> anillo_fp_x
# input: f --> polinomio fabricado por fpx (objeto opaco) no nulo
# output: g = producto de los factores irreducibles mónicos de f, es decir,
# si f = c * f1^e1 * f2^e2 * ... * fr^er con los fi irreducibles mónicos
# distintos entre si, ei >= 1, c en fp, entonces g = f1 * f2 * ... * fr
def sqfree_fact_fpx(fpx, f):
    f = fpx.reducir(f)  # aseguramos que esta en mod p y sin ceros al final

    # f = 0
    if fpx.es_cero(f):
        raise ValueError("f no puede ser el polinomio nulo")
    
    # Grado(f) <= 0 --> es constante
    if fpx.grado(f) <= 0:
        return fpx.uno()    # al hacer monico, devolvemos 1
    
    p   = fpx.p

    # Calculo la derivada
    # f = a0 + a1 * x + a2 * x^2 + ··· --> f' = a1 + 2 * a2 * x + ···
    # f = (a0, a1, a2, ..., an) --> f' = (a1, 2*a2, 3*a3, ..., n*an)
    deg = fpx.grado(f)

    if deg <= 0:
        df = fpx.cero()
    else:
        df_list = [(i*f[i]) % p for i in range(1, deg + 1)] # recorre los indices desde 1 hasta el grado
        df = fpx.reducir(tuple(df_list))

    # Caso 1. Derivada nula --> f(x) = g(x^p)
    if fpx.es_cero(df):
        g = tuple(f[i] for i in range(0, len(f), p))
        g = fpx.reducir(g)

        # Recursion: f(x) = g(x^p)
        return sqfree_fact_fpx(fpx, g)
    
    # Caso 2. Derivada no nula
    gcd = fpx.gcd(f, df)  # gcd(f, f')
    g = fpx.div(f, gcd)   # q = f/gcd(f, f') y g = q reducido => es square - free
                
    # hacer g monico
    if fpx.grado(g) > 0:
        g = fpx.hacer_monico(g)

    # eliminar de f todos los factores contenidos en g
    residuo = f
    while True:
        y = fpx.gcd(residuo, g)
        if fpx.grado(y) <= 0:
            break
        residuo = fpx.div(residuo, y)
    
    # Si hay residuo entonces es una potencia p-esima
    if fpx.grado(residuo) > 0:

        # extraer h(x) tal que residuo(x) = h(x^p)
        h = tuple(residuo[i] for i in range(0, len(residuo), p))
        h = fpx.reducir(h)

        # recursion
        g2 = sqfree_fact_fpx(fpx, h)

        # unir factores: lcm(g, g2) = g * g2 / gcd(g, g2)
        producto = fpx.mult(g, g2)
        interseccion = fpx.gcd(g, g2)

        return fpx.div(producto, interseccion)
    
    else:
        # no hay potencia p-esima -> devuelve g
        return g


# distinct-degree factorization
# input: fpx --> anillo_fp_x
# input: g --> polinomio de fpx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos cada uno con multiplicidad uno
# output: [h1, h2, ..., hr], donde hi = producto de los factores irreducibles
# mónicos de h de grado = i, el último hr debe ser no nulo y por supuesto
# g = h1 * h2 * ... * hr
def didegr_fact_fpx(fpx, g):
    g = fpx.hacer_monico(fpx.reducir(g))
    factores = []
    x = (0, 1)
    h = g
    d = 1

    # Algoritmo clásico de factorización por grados distintos
    while fpx.grado(h) >= 2 * d:
        x_pd = fpx.pot_mod(x, fpx.p**d, h)
        g_d = fpx.gcd(fpx.suma(x_pd, fpx.inv_adit(x)), h)
        if not fpx.es_uno(g_d):
            g_d = fpx.hacer_monico(g_d)
            factores.append((g_d, d))
            h = fpx.div(h, g_d)   # quitamos el bloque de grado d
            # repetimos con el mismo d hasta que g_d = 1
        else:
            d += 1

    # Lo que quede al final (grado > 0) es irreducible de grado >= d
    if not fpx.es_uno(h) and fpx.grado(h) > 0:
        h = fpx.hacer_monico(h)
        factores.append((h, fpx.grado(h)))

    # Si no se halló ningún factor (p.ej., g es irreducible), devolvemos [g]
    if not factores:
        return [g]

    # Construir la lista [h1, ..., hr], rellenando con 1 donde no haya factores
    maxd = max(d for _, d in factores)
    hlist = [fpx.uno() for _ in range(maxd)]
    for fac, d in factores:
        hlist[d - 1] = fpx.mult(hlist[d - 1], fac)

    return hlist

# equal-degree factorization
# input: fpx --> anillo_fp_x (supondremos p impar)
# input: r --> int
# input: h --> polinomio de fpx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos de grado r con multiplicidad uno
# output: [u1, ..., us], donde h = u1 * u2* ... * us y los ui son irreducibles
# mónicos de grado = r
def eqdegr_fact_fpx(fpx, r, h):
    # Auxiliares internas
    def random_poly(max_degree):
        while True:
            coefs = [fpx.fp.aleatorio() for _ in range(max_degree + 1)]
            poly = fpx.elem_de_tuple(tuple(coefs))
            if not fpx.es_cero(poly):
                return poly

    def find_linear_factors(f):
        factores = []
        h = f
        for i in range(fpx.p):
            x_minus_a = fpx.suma((0, 1), fpx.elem_de_tuple(((-i) % fpx.p,)))
            g = fpx.gcd(h, x_minus_a)
            if not fpx.es_uno(g):
                factores.append(fpx.hacer_monico(g))
                h = fpx.div(h, g)
                if fpx.es_uno(h):
                    break
        return factores

    # --- inicio del algoritmo principal ---
    h = fpx.hacer_monico(fpx.reducir(h))
    if fpx.grado(h) <= r:
        return [fpx.hacer_monico(h)]
    if r == 1:
        return find_linear_factors(h)

    factores = [h]
    max_iter = 50  # límite de seguridad

    for _ in range(max_iter):
        changed = False
        for fi in list(factores):
            if fpx.grado(fi) > r:
                # probar varios polinomios aleatorios
                for _ in range(10):
                    a = random_poly(fpx.grado(fi) - 1)
                    if fpx.es_cero(a):
                        continue

                    g = fpx.pot_mod(a, (fpx.p**r - 1)//2, fi)
                    g = fpx.suma(g, fpx.inv_adit(fpx.uno()))
                    h_gcd = fpx.gcd(fi, g)

                    if not fpx.es_uno(h_gcd) and not fpx.es_igual(h_gcd, fi):
                        h_gcd = fpx.hacer_monico(h_gcd)
                        coc = fpx.hacer_monico(fpx.div(fi, h_gcd))
                        factores.remove(fi)
                        factores.append(h_gcd)
                        factores.append(coc)
                        changed = True
                        break
        if not changed:
            break

    # refinamiento final: todos los factores deben ser de grado r
    final_factores = []
    for fi in factores:
        if fpx.grado(fi) == r:
            final_factores.append(fpx.hacer_monico(fi))
        elif fpx.grado(fi) > r:
            final_factores.extend(eqdegr_fact_fpx(fpx, r, fi))
        else:
            final_factores.append(fpx.hacer_monico(fi))

    return final_factores

# multiplicidad de factor irreducible mónico
# input: fpx --> anillo_fp_x
# input: f --> polinomio de fpx (objeto opaco) no nulo
# input: u --> polinomio irreducible mónico (objeto opaco) de grado >= 1
# output: multiplicidad de u como factor de f, es decir, el entero e >= 0
# mas grande tal que u^e | f
def multiplicidad_fpx(fpx, f, u):
    f = fpx.reducir(f)
    u = fpx.hacer_monico(fpx.reducir(u))
    e = 0
    while True:
        g = fpx.gcd(f, u)
        if fpx.grado(g) == 0:  # gcd = 1 => ya no divide
            break
        e += 1
        f = fpx.div(f, u)
    return e

# factorización de Cantor-Zassenhaus
# input: fpx --> anillo_fp_x (supondremos p impar)
# input: f --> polinomio de fpx (objeto opaco)
# output: [(f1,e1), ..., (fr,er)] donde f = c * f1^e1 * ... * fr^er es la
# factorización completa de f en irreducibles mónicos fi con multiplicidad
# ei >= 1 y los fi son distintos entre si y por supuesto c es el coeficiente
# principal de f
def fact_fpx(fpx, f):                     # mantener esta implementación
    g = sqfree_fact_fpx(fpx, f)
    h = didegr_fact_fpx(fpx, g)
    irreducibles = []
    for r in range(len(h)):
        if fpx.grado(h[r]) > 0:
            irreducibles += eqdegr_fact_fpx(fpx, r+1, h[r])
    factorizacion = []
    for u in irreducibles:
        e = multiplicidad_fpx(fpx, f, u)
        factorizacion += [(u,e)]
    return factorizacion



# ==========================================
# FACTORIZACIÓN CANTOR-ZASSENHAUS Fq[x]
# ==========================================

# square-free factorization
# input: fqx --> anillo_fq_x
# input: f --> polinomio fabricado por fqx (objeto opaco) no nulo
# output: g = producto de los factores irreducibles mónicos de f, es decir,
# si f = c * f1^e1 * f2^e2 * ... * fr^er con los fi irreducibles mónicos
# distintos entre si, ei >= 1, c en fq, entonces g = f1 * f2 * ... * fr
def sqfree_fact_fqx(fqx, f):
    f = fqx.reducir(f)  # aseguramos que esta en mod q y sin ceros al final

    # f = 0
    if fqx.es_cero(f):
        raise ValueError("f no puede ser el polinomio nulo")
    
    # Grado(f) <= 0 --> es constante
    if fqx.grado(f) <= 0:
        return fqx.uno()    # al hacer monico, devolvemos 1
    
    p   = fqx.p
    q   = fqx.q

    # calcula la derivada
    deg = fqx.grado(f)

    if deg <= 0:
        df = fqx.cero()
    else:
        df_coef = []
        for i in range(1, deg + 1):
            coef_i = f[i]
            k = i % p
            if k == 0:
                df_coef.append(fqx.fq.cero())   # coef nulo (equvale a (0,))
            else:
                k_fq = fqx.fq.elem_de_int(k)
                df_coef.append(fqx.fq.mult(coef_i, k_fq))
        df = fqx.reducir(tuple(df_coef))
    
    # Caso 1. Derivada nula --> f(x) = g(x^p)
    if fqx.es_cero(df):
        exp = q // p # raiz p-esima en Fq : a -> a^(q/p)

        g = tuple(fqx.fq.pot(f[i], exp) for i in range(0, len(f), p))
        g = fqx.reducir(g)

        # Recursion: f(x) = g(x^p)
        return sqfree_fact_fqx(fqx, g)
    
    # Caso 2. Derivada no nula
    gcd = fqx.gcd(f, df)  # gcd(f, f')
    g = fqx.div(f, gcd)   # q = f/gcd(f, f') y g = q reducido => es square - free

    # hacer g monico
    if fqx.grado(g) > 0:
        g = fqx.hacer_monico(g)

    # eliminar de f todos los factores contenidos en g
    residuo = f
    while True:
        y = fqx.gcd(residuo, g)
        if fqx.grado(y) <= 0:
            break
        residuo = fqx.div(residuo, y)

    # Si hay residuo entonces es una potencia p-esima
    if fqx.grado(residuo) > 0:
        exp = q // p # raiz p-esima en Fq : a -> a^(q/p)

        # extraer h(x) tal que residuo(x) = h(x^p)
        h = tuple(fqx.fq.pot(residuo[i], exp) for i in range(0, len(residuo), p))
        h = fqx.reducir(h)

        # recursion
        g2 = sqfree_fact_fqx(fqx, h)

        # unir factores: lcm(g, g2) = g * g2 / gcd(g, g2)
        prod = fqx.mult(g, g2)
        inter = fqx.gcd(g, g2)

        return fqx.div(prod, inter)
    
    else:
        # no hay potencia p-esima -> devuelve g
        return g
    

# distinct-degree factorization
# input: fqx --> anillo_fq_x
# input: g --> polinomio de fqx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos cada uno con multiplicidad uno
# output: [h1, h2, ..., hr], donde hi = producto de los factores irreducibles
# mónicos de h de grado = i, el último hr debe ser no nulo y por supuesto
# g = h1 * h2 * ... * hr
def didegr_fact_fqx(fqx, g):
    g = fqx.hacer_monico(fqx.reducir(g))
    factores = []
    x = (fqx.fq.cero(), fqx.fq.uno())  # polinomio x
    h = g
    d = 1

    # Algoritmo clásico: gcd(h, x^{q^d} - x)
    while fqx.grado(h) >= 2 * d:
        x_qd = fqx.pot_mod(x, fqx.q**d, h)
        g_d = fqx.gcd(fqx.suma(x_qd, fqx.inv_adit(x)), h)  # g_d = gcd(x^{q^d} - x, h)
        if not fqx.es_uno(g_d):
            g_d = fqx.hacer_monico(g_d)
            factores.append((g_d, d))
            h = fqx.div(h, g_d)      # quitamos el bloque de grado d
            # se repite con el mismo d hasta que g_d = 1
        else:
            d += 1

    # Lo que quede al final (grado > 0) es irreducible de grado >= d
    if not fqx.es_uno(h) and fqx.grado(h) > 0:
        h = fqx.hacer_monico(h)
        factores.append((h, fqx.grado(h)))

    # Si no se halló ningún factor (por ejemplo si g es irreducible), devolvemos [g]
    if not factores:
        return [g]

    # Construir la lista [h1, ..., hr], rellenando con 1 donde no haya factores
    maxd = max(d for _, d in factores)
    hlist = [fqx.uno() for _ in range(maxd)]
    for fac, d in factores:
        hlist[d - 1] = fqx.mult(hlist[d - 1], fac)

    return hlist

# equal-degree factorization
# input: fqx --> anillo_fq_x (supondremos q impar)
# input: r --> int
# input: h --> polinomio de fqx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos de grado r con multiplicidad uno
# output: [u1, ..., us], donde h = u1 * u2* ... * us y los ui son irreducibles
# mónicos de grado = r
def eqdegr_fact_fqx(fqx, r, h):
    # Auxiliares internas (adaptadas de tus métodos)
    def random_poly(max_degree):
        # polinomio aleatorio no nulo de grado <= max_degree
        while True:
            coefs = [fqx.fq.aleatorio() for _ in range(max_degree + 1)]
            poly = fqx.elem_de_tuple(tuple(coefs))
            if not fqx.es_cero(poly):
                return poly

    def find_linear_factors(f):
        # solo se usa cuando r = 1: buscamos raíces en F_q
        factores = []
        h_local = f
        for i in range(fqx.q):
            a = fqx.fq.elem_de_int(i)
            x_minus_a = fqx.suma(
                (fqx.fq.cero(), fqx.fq.uno()),      # x
                (fqx.fq.inv_adit(a),)               # -a
            )
            g = fqx.gcd(h_local, x_minus_a)
            if not fqx.es_uno(g):
                factores.append(fqx.hacer_monico(g))
                h_local = fqx.div(h_local, g)
                if fqx.es_uno(h_local):
                    break
        return factores

    # --- inicio del algoritmo principal ---
    h = fqx.hacer_monico(fqx.reducir(h))

    if fqx.grado(h) <= r:
        return [fqx.hacer_monico(h)]
    if r == 1:
        return find_linear_factors(h)

    factores = [h]
    max_iter = 50  # límite de seguridad

    for _ in range(max_iter):
        changed = False
        for fi in list(factores):
            if fqx.grado(fi) > r:
                # probar varios polinomios aleatorios
                for _ in range(10):
                    a = random_poly(fqx.grado(fi) - 1)
                    if fqx.es_cero(a):
                        continue

                    # g = a^{(q^r - 1)/2} - 1  (mod fi)
                    exp = (fqx.q**r - 1) // 2
                    g = fqx.pot_mod(a, exp, fi)
                    g = fqx.suma(g, fqx.inv_adit(fqx.uno()))
                    h_gcd = fqx.gcd(fi, g)

                    if (not fqx.es_uno(h_gcd)) and (not fqx.es_igual(h_gcd, fi)):
                        h_gcd = fqx.hacer_monico(h_gcd)
                        coc = fqx.hacer_monico(fqx.div(fi, h_gcd))
                        factores.remove(fi)
                        factores.append(h_gcd)
                        factores.append(coc)
                        changed = True
                        break
        if not changed:
            break

    # refinamiento final: todos los factores deben ser de grado r
    final_factores = []
    for fi in factores:
        if fqx.grado(fi) == r:
            final_factores.append(fqx.hacer_monico(fi))
        elif fqx.grado(fi) > r:
            final_factores.extend(eqdegr_fact_fqx(fqx, r, fi))
        else:
            final_factores.append(fqx.hacer_monico(fi))

    return final_factores

# multiplicidad de factor irreducible mónico
# input: fqx --> anillo_fq_x
# input: f --> polinomio de fqx (objeto opaco) no nulo
# input: u --> polinomio irreducible mónico (objeto opaco) de grado >= 1
# output: multiplicidad de u como factor de f, es decir, el entero e >= 0
# mas grande tal que u^e | f
def multiplicidad_fqx(fqx, f, u):
    f = fqx.reducir(f)
    u = fqx.hacer_monico(fqx.reducir(u))
    e = 0

    while True:
        g = fqx.gcd(f, u)
        # si gcd = 1 (grado 0), u ya no divide a f
        if fqx.grado(g) == 0:
            break
        e += 1
        f = fqx.div(f, u)

    return e

# factorización de Cantor-Zassenhaus
# input: fqx --> anillo_fq_x (supondremos q impar)
# input: f --> polinomio de fqx (objeto opaco)
# output: [(f1,e1), ..., (fr,er)] donde f = c * f1^e1 * ... * fr^er es la
# factorización completa de f en irreducibles mónicos fi con multiplicidad
# ei >= 1 y los fi son distintos entre si y por supuesto c es el coeficiente
# principal de f
def fact_fqx(fqx, f):                     # mantener esta implementación
    g = sqfree_fact_fqx(fqx, f)
    h = didegr_fact_fqx(fqx, g)
    irreducibles = []
    for r in range(len(h)):
        if fqx.grado(h[r]) > 0:
            irreducibles += eqdegr_fact_fqx(fqx, r+1, h[r])
    factorizacion = []
    for u in irreducibles:
        e = multiplicidad_fqx(fqx, f, u)
        factorizacion += [(u,e)]
    return factorizacion

