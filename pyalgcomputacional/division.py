from .karatsuba import fp_x_mult_karatsuba, fq_x_mult_karatsuba
from .toeplitz import fp_toep_sup_inv, fp_toep_inf_vec, fq_toep_sup_inv, fq_toep_inf_vec


# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio no nulo (objeto opaco creado por fpx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fp_x_divmod(fpx, f, g):
    f = fpx.reducir(f)
    g = fpx.reducir(g)

    if fpx.es_cero(g):
        raise ValueError("No se puede dividir por cero")
    
    # Caso deg(f) < deg(g) => q = 0, r = f
    deg_f = fpx.grado(f)
    deg_g = fpx.grado(g)

    if deg_f < deg_g:
        return (fpx.cero(), f)
    
    # Caso deg(g) < deg(f)
    k = deg_f - deg_g + 1
    fp = fpx.fp

    # Dar la vuelta a los polinomios
    f_rev = f[::-1]
    g_rev = g[::-1]

    if len(g_rev) >= k:
        g_inv = tuple(g_rev[:k])
    else:
        g_inv = g_rev + (fp.cero(),) * (k - len(g_rev))

    if len(f_rev) >= k:
        f_inv = tuple(f_rev[:k])
    else:
        f_inv = f_rev + (fp.cero(),) * (k - len(f_rev))
    
    # Calcular la inversa de g
    inv_g_rev = fp_toep_sup_inv(fp, k, g_inv)

    # Calcular el cociente
    q_rev = fp_toep_inf_vec(fp, k, inv_g_rev, f_inv)

    # Dar la vuelta al cociente
    q = q_rev[::-1]
    q = fpx.reducir(q)

    # Calcular el resto r = f - g*q
    g_q = fp_x_mult_karatsuba(fpx, g, q)
    r = fpx.suma(f, fpx.inv_adit(g_q))

    return (q, r)



# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio no nulo (objeto opaco creado por fqx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fq_x_divmod(fqx, f, g):
    f = fqx.reducir(f)
    g = fqx.reducir(g)

    if fqx.es_cero(g):
        raise ValueError("No se puede dividir por cero")
    
    # Caso deg(f) < deg(g) => q = 0, r = f
    deg_f = fqx.grado(f)
    deg_g = fqx.grado(g)

    if deg_f < deg_g:
        return (fqx.cero(), f)
    
    # Caso deg(g) < deg(f)
    k = deg_f - deg_g + 1
    fq = fqx.fq

    # Dar la vuelta a los polinomios
    f_rev = f[::-1]
    g_rev = g[::-1]

    if len(g_rev) >= k:
        g_inv = tuple(g_rev[:k])
    else:
        g_inv = g_rev + (fq.cero(),) * (k - len(g_rev))

    if len(f_rev) >= k:
        f_inv = tuple(f_rev[:k])
    else:
        f_inv = f_rev + (fq.cero(),) * (k - len(f_rev))
    
    # Calcular la inversa de g
    inv_g_rev = fq_toep_sup_inv(fq, k, g_inv)

    # Calcular el cociente
    q_rev = fq_toep_inf_vec(fq, k, inv_g_rev, f_inv)

    # Dar la vuelta al cociente
    q = q_rev[::-1]
    q = fqx.reducir(q)

    # Calcular el resto r = f - g*q
    g_q = fq_x_mult_karatsuba(fqx, g, q)
    r = fqx.suma(f, fqx.inv_adit(g_q))

    return (q, r)

