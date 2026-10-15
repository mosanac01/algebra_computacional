from .anillo_fp_x import anillo_fp_x
from .anillo_fq_x import anillo_fq_x

# ==========================================
# MATRICES DE TOEPLITZ SOBRE Fp
# ==========================================

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_inf_vec(fp, n, a, b):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != n or len(b) != n:
        raise ValueError("Las longitudes de a y b deben ser iguales a n")
    
    # Construir Fp[x]
    fpx = cf.anillo_fp_x(fp)

    # Pasamos a y b a polinomios
    A = fpx.elem_de_tuple(a)
    B = fpx.elem_de_tuple(b)

    # C = A * B
    C = fp_x_mult_karatsuba(fpx, A, B)

    # Extraer los n primeros coeficientes de C
    resultado = list(C)
    while len(resultado) < n:
        resultado.append(fp.cero())

    return tuple(resultado[:n])


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_sup_vec(fp, n, a, b):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != n or len(b) != n:
        raise ValueError("Las longitudes de a y b deben ser iguales a n")
    
    # Construir Fp[x]
    fpx = cf.anillo_fp_x(fp)

    # Pasamos a y b a polinomios
    A = fpx.elem_de_tuple(a[::-1])  # hemos invertido a
    B = fpx.elem_de_tuple(b)

    # C = A * B
    C = fp_x_mult_karatsuba(fpx, A, B)

    # Extraer los n ultimos coeficientes de C
    resultado = list(C)
    while len(resultado) < 2*n - 1:
        resultado.append(fp.cero())

    resultado = resultado[n-1:2*n-1]   

    return tuple(resultado)

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fp (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_vec(fp, n, a, b):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != 2*n - 1:
        raise ValueError("a debe tener longitud 2*n - 1")
    if len(b) != n:
        raise ValueError("b debe tener longitud n")
    
    # Primera fila y primera columna
    fila = a[:n]
    # columna = (a[0],) + a[n:]
    columna = (fp.cero(),) + a[n:]

    # Superior
    v_sup = fp_toep_sup_vec(fp, n, fila, b)

    # Inferior
    v_inf = fp_toep_inf_vec(fp, n, columna, b)

    # Correccion diagonal
    a0 = a[0]
    resultado = []
    for i in range(n):
        valor = fp.suma(v_sup[i], v_inf[i])
        valor = fp.suma(valor, fp.inv_adit(fp.mult(a0, b[i])))
        resultado.append(valor)

    return tuple(resultado)


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fp_toep_inf_inv(fp, n, a):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != n:
        raise ValueError("a debe tener longitud n")
    if fp.es_cero(a[0]):
        raise ValueError("a[0] no debe ser cero")
    
    # Caso base n=1
    if n == 1:
        return (fp.inv_mult(a[0]),)
    
    # T = | T_0   0  |
    #     | T_1  T_0 |
    # Dividir a en dos mitades
    k = (n + 1) // 2
    m = n - k

    # Calculamos la primera columna de T_0^{-1}
    T_0 = a[:k]
    v = fp_toep_inf_inv(fp, k, T_0)

    # Caso m = 0
    if m == 0:
        return v
    
    # Ahora hay que recordar que T_1 puede no ser cuadrada
    # Calculamos el residuo y = T_1 * v
    # Como T_1 no es cuadrada, usamos la matriz T completa
    # que multiplicamos por v extendido con ceros
    v_ext = list(v)
    while len(v_ext) < n:
        v_ext.append(fp.cero())
    v_ext = tuple(v_ext)

    # Multiplicamos por T para obtener (T_0 * v, T_1 * v)^t = (e_1, y)^t
    prod = fp_toep_inf_vec(fp, n, a, v_ext)
    y = prod[k:]

    # Calculamos u = -T_0^{-1} * y
    # como m <= k rellenamos y con ceros
    y_ext = list(y)
    while len(y_ext) < k:
        y_ext.append(fp.cero())
    y_ext = tuple(y_ext)

    # Multiplicamos por T_0^{-1}
    w = fp_toep_inf_vec(fp, k, v, y_ext)

    # El resultado son los primeros m elementos de w cambiados de signo
    resultado = []
    for i in range(m):
        resultado.append(fp.inv_adit(w[i]))

    # Finalmente concatenamos con v:
    return v + tuple(resultado)
  
# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fp_toep_sup_inv(fp, n, a):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != n:
        raise ValueError("a debe tener longitud n")
    if fp.es_cero(a[0]):
        raise ValueError("a[0] no debe ser cero")
    
    # Caso base n=1
    if n == 1:
        return (fp.inv_mult(a[0]),)
    
    # T = | T_0  T_1 |
    #     |  0   T_0 |
    # Dividir a en dos mitades
    k = (n + 1) // 2
    m = n - k

    # Calculamos la primera fila de T_0^{-1}
    T_0 = a[:k]
    v = fp_toep_sup_inv(fp, k, T_0)

    # Caso m = 0
    if m == 0:
        return v
    
    # Ahora hay que recordar que T_1 puede no ser cuadrada
    # Calculamos el residuo y = v * T_1
    # Como T_1 no es cuadrada, usamos la matriz T completa
    # que multiplicamos por v extendido con ceros
    v_ext = list(v)
    while len(v_ext) < n:
        v_ext.append(fp.cero())

    v_ext = tuple(v_ext[::-1])   # Y le damos la vuelta

    # Multiplicamos por T (completa)
    prod = fp_toep_sup_vec(fp, n, a, v_ext)
    y = prod[::-1]  # Volvemos a darle la vuelta
    y = y[k:]    

    # Calculamos u = - y * T_0^{-1}
    # como m <= k rellenamos y con ceros
    y_ext = list(y)
    while len(y_ext) < k:
        y_ext.append(fp.cero())

    y_ext = tuple(y_ext[::-1])  # Y damos la vuelta

    # Multiplicamos por T_0^{-1}
    w = fp_toep_sup_vec(fp, k, v, y_ext)
    w = w[::-1]     # Volvemos a darle la vuelta

    # El resultado son los primeros m elementos de w cambiados de signo
    resultado = []
    for i in range(m):
        resultado.append(fp.inv_adit(w[i]))

    # Finalmente concatenamos con v:
    return v + tuple(resultado)


# ==========================================
# MATRICES DE TOEPLITZ SOBRE Fq
# ==========================================


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_inf_vec(fq, n, a, b):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != n or len(b) != n:
        raise ValueError("Las longitudes de a y b deben ser iguales a n")
    
    # Construir Fq[x]
    fqx = cf.anillo_fq_x(fq)

    # Pasamos a y b a polinomios
    A = fqx.elem_de_tuple(a)
    B = fqx.elem_de_tuple(b)

    # C = A * B
    C = fq_x_mult_karatsuba(fqx, A, B)

    # Extraer los n primeros coeficientes de C
    resultado = list(C)
    while len(resultado) < n:
        resultado.append(fq.cero())

    return tuple(resultado[:n])


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_sup_vec(fq, n, a, b):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != n or len(b) != n:
        raise ValueError("Las longitudes de a y b deben ser iguales a n")
    
    # Construir Fq[x]
    fqx = cf.anillo_fq_x(fq)

    # Pasamos a y b a polinomios
    A = fqx.elem_de_tuple(a[::-1])  # hemos invertido a
    B = fqx.elem_de_tuple(b)

    # C = A * B
    C = fq_x_mult_karatsuba(fqx, A, B)

    # Extraer los n ultimos coeficientes de C
    resultado = list(C)
    while len(resultado) < 2*n - 1:
        resultado.append(fq.cero())

    return tuple(resultado[n-1:2*n-1])


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fq (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_vec(fq, n, a, b):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != 2*n - 1:
        raise ValueError("a debe tener longitud 2*n - 1")
    if len(b) != n:
        raise ValueError("b debe tener longitud n")
    
    # Primera fila y primera columna
    fila = a[:n]
    # columna = (a[0],) + a[n:]
    columna = (fq.cero(),) + a[n:]

    # Superior
    v_sup = fq_toep_sup_vec(fq, n, fila, b)

    # Inferior
    v_inf = fq_toep_inf_vec(fq, n, columna, b)

    # Correccion diagonal
    a0 = a[0]
    resultado = []
    for i in range(n):
        valor = fq.suma(v_sup[i], v_inf[i])
        valor = fq.suma(valor, fq.inv_adit(fq.mult(a0, b[i])))
        resultado.append(valor)

    return tuple(resultado)


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fq_toep_inf_inv(fq, n, a):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != n:
        raise ValueError("a debe tener longitud n")
    if fq.es_cero(a[0]):
        raise ValueError("a[0] no debe ser cero")
    
    # Caso base n=1
    if n == 1:
        return (fq.inv_mult(a[0]),)

    # Dividir a en dos mitades
    k = (n + 1) // 2
    m = n - k

    # Calculamos la primera columna de T_0^{-1}
    T_0 = a[:k]
    v = fq_toep_inf_inv(fq, k, T_0)

    # Caso m = 0
    if m == 0:
        return v
    
    # Calculamos el residuo y = T_1 * v
    v_ext = list(v)
    while len(v_ext) < n:
        v_ext.append(fq.cero())
    v_ext = tuple(v_ext)

    # Multiplicamos por T para obtener (T_0 * v, T_1 * v)^t = (e_1, y)^t
    prod = fq_toep_inf_vec(fq, n, a, v_ext)
    y = prod[k:]

    # Calculamos u = -T_0^{-1} * y 
    # como m <= k rellenamos y con ceros
    y_ext = list(y)
    while len(y_ext) < k:
        y_ext.append(fq.cero())
    y_ext = tuple(y_ext)

    # Multiplicamos por T_0^{-1}
    w = fq_toep_inf_vec(fq, k, v, y_ext)

    # El resultado son los primeros m elementos de w cambiados de signo
    resultado = []
    for i in range(m):
        resultado.append(fq.inv_adit(w[i]))

    # Finalmente concatenamos con v:
    return v + tuple(resultado)
  

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fq_toep_sup_inv(fq, n, a):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n debe ser un entero mayor o igual a 1")
    if len(a) != n:
        raise ValueError("a debe tener longitud n")
    if fq.es_cero(a[0]):
        raise ValueError("a[0] no debe ser cero")
    
    # Caso base n=1
    if n == 1:
        return (fq.inv_mult(a[0]),)
    
    # Dividir a en dos mitades
    k = (n + 1) // 2
    m = n - k

    # Calculamos la primera fila de T_0^{-1}
    T_0 = a[:k]
    v = fq_toep_sup_inv(fq, k, T_0)

    # Caso m = 0
    if m == 0:
        return v
    
    # Ahora hay que recordar que T_1 puede no ser cuadrada
    # Calculamos el residuo y = v * T_1
    # Como T_1 no es cuadrada, usamos la matriz T completa
    # que multiplicamos por v extendido con ceros
    v_ext = list(v)
    while len(v_ext) < n:
        v_ext.append(fq.cero())

    v_ext = tuple(v_ext[::-1])   # Y le damos la vuelta

    # Multiplicamos por T (completa)
    prod = fq_toep_sup_vec(fq, n, a, v_ext)
    y = prod[::-1]  # Volvemos a darle la vuelta
    y = y[k:]    

    # Calculamos u = - y * T_0^{-1}
    # como m <= k rellenamos y con ceros
    y_ext = list(y)
    while len(y_ext) < k:
        y_ext.append(fq.cero())

    y_ext = tuple(y_ext[::-1])  # Y damos la vuelta

    # Multiplicamos por T_0^{-1}
    w = fq_toep_sup_vec(fq, k, v, y_ext)
    w = w[::-1]     # Volvemos a darle la vuelta

    # El resultado son los primeros m elementos de w cambiados de signo
    resultado = []
    for i in range(m):
        resultado.append(fq.inv_adit(w[i]))

    # Finalmente concatenamos con v:
    return v + tuple(resultado)

