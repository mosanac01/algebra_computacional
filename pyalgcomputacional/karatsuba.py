# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio (objeto opaco creado por fpx)
# output: f*g calculado usando el método de Karatsuba
def fp_x_mult_karatsuba(fpx, f, g):
    f = fpx.reducir(f)
    g = fpx.reducir(g)

    # Caso base (alguno es cero)
    if fpx.es_cero(f) or fpx.es_cero(g):
        return fpx.cero()
    
    deg_f = fpx.grado(f)
    deg_g = fpx.grado(g)

    # Caso base (grado pequeno)
    if deg_f <= 1 or deg_g <= 1:
        return fpx.mult(f, g)   # multiplicacion normal
    
    # Grado maximo
    n = max(deg_f, deg_g)

    # Si n es pequeno usar multiplicacion normal
    if n <= 32:
        return fpx.mult(f, g)
    
    # Mitad del grado maximo
    m = n // 2

    # Dividir f
    f0 = fpx.reducir(f[:m])
    f1 = fpx.reducir(f[m:])

    # Dividir g
    g0 = fpx.reducir(g[:m])
    g1 = fpx.reducir(g[m:])
    
    # Recursion
    # z0 = f0*g0
    # z2 = f1*g1
    z0 = fp_x_mult_karatsuba(fpx, f0, g0)
    z2 = fp_x_mult_karatsuba(fpx, f1, g1)

    # z1 = (f0+f1)*(g0+g1) - z0 - z2
    f_sum = fpx.suma(f0, f1)
    g_sum = fpx.suma(g0, g1)

    z1 = fp_x_mult_karatsuba(fpx, f_sum, g_sum)
    z1 = fpx.suma(z1, fpx.inv_adit(z0))
    z1 = fpx.suma(z1, fpx.inv_adit(z2))

    # Combinar resultados z0 + z1*x^m + z2*x^(2m)
    resultado = [0] * (max(len(z0), m + len(z1), 2 * m + len(z2)))

    for i in range(len(z0)):
        resultado[i] = fpx.fp.suma(resultado[i], z0[i])

    for i in range(len(z1)):
        resultado[i + m] = fpx.fp.suma(resultado[i + m], z1[i])

    for i in range(len(z2)):
        resultado[i + 2 * m] = fpx.fp.suma(resultado[i + 2 * m], z2[i])

    # Reducir el resultado:
    return fpx.reducir(resultado)    



# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio (objeto opaco creado por fqx)
# output: f*g calculado usando el método de Karatsuba
def fq_x_mult_karatsuba(fqx, f, g):
    f = fqx.reducir(f)
    g = fqx.reducir(g)

    # Caso base
    if fqx.es_cero(f) or fqx.es_cero(g):
        return fqx.cero()
    
    deg_f = fqx.grado(f)
    deg_g = fqx.grado(g)

    # Caso base (grado pequeno)
    if deg_f <= 1 or deg_g <= 1:
        return fqx.mult(f, g)   # multiplicacion normal
    
    # Grado maximo
    n = max(deg_f, deg_g)

    # Si n es pequeno usar multiplicacion normal
    if n <= 32:
        return fqx.mult(f, g)
    
    # Mitad del grado maximo
    m = n // 2

    # Dividir f
    f0 = fqx.reducir(f[:m])
    f1 = fqx.reducir(f[m:])

    # Dividir g
    g0 = fqx.reducir(g[:m])
    g1 = fqx.reducir(g[m:])
    
    # Recursion
    # z0 = f0*g0
    # z2 = f1*g1
    z0 = fq_x_mult_karatsuba(fqx, f0, g0)
    z2 = fq_x_mult_karatsuba(fqx, f1, g1)

    # z1 = (f0+f1)*(g0+g1) - z0 - z2
    f_sum = fqx.suma(f0, f1)
    g_sum = fqx.suma(g0, g1)

    z1 = fq_x_mult_karatsuba(fqx, f_sum, g_sum)
    z1 = fqx.suma(z1, fqx.inv_adit(z0))
    z1 = fqx.suma(z1, fqx.inv_adit(z2))

    # Combinar resultados z0 + z1*x^m + z2*x^(2m)
    resultado = [fqx.fq.cero()] * (max(len(z0), m + len(z1), 2 * m + len(z2)))

    for i in range(len(z0)):
        resultado[i] = fqx.fq.suma(resultado[i], z0[i])

    for i in range(len(z1)):
        resultado[i + m] = fqx.fq.suma(resultado[i + m], z1[i])

    for i in range(len(z2)):
        resultado[i + 2 * m] = fqx.fq.suma(resultado[i + 2 * m], z2[i])

    # Reducir el resultado:
    return fqx.reducir(resultado)    

