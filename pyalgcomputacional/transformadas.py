# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# utilizar el algoritmo de Cooley-Tuckey
def fp_fft(fp, g, k, a):
    # Caso base: si k = 0, n = 1 => devuelve el propio elemento
    if k == 0:
        return a
    
    n = len(a)
    if n != 2**k:
        raise ValueError("La longitud de a debe ser igual a n = 2**k")
    
    m = n // 2

    # Dividir a en pares e impares
    a_par = a[0::2]
    a_impar = a[1::2]

    # Recursion
    # calcular g^2
    g2 = fp.mult(g, g)

    # Evaluamos en pares e impares
    y0 = fp_fft(fp, g2, k - 1, a_par)
    y1 = fp_fft(fp, g2, k - 1, a_impar)

    # Combinar resultados
    resultado = [fp.cero()] * n
    w = fp.uno()  # w = g^0
    for j in range(m):
        t = fp.mult(w, y1[j]) # t = w * y1[j]
        u = y0[j]
        
        resultado[j] = fp.suma(u, t)    # u + t
        resultado[j + m] = fp.suma(u, fp.inv_adit(t))   # u - t

        w = fp.mult(w, g)  # w = w * g

    return tuple(resultado)


# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fp_ifft(fp, g, k, a):
    # Inversa de g => g^-1
    inv_g = fp.inv_mult(g)

    # Calcular DFT_{n,g^-1}(a)
    dft = fp_fft(fp, inv_g, k, a)

    # Calcular n^(-1)
    # Como n = 2^k => n^{-1} = (2^k)^{-1} = (2^{-1})^k
    dos = fp.suma(fp.uno(), fp.uno())
    inv_dos = fp.inv_mult(dos)
    inv_n = fp.pot(inv_dos, k)

    # Multiplicar cada elemento por n^(-1)
    resultado = []
    for x in dft:
        resultado.append(fp.mult(x, inv_n))

    return tuple(resultado)



# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a q-1
# input: a -> tupla de longitud n de elementos de fq
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# utilizar el algoritmo de Cooley-Tuckey
def fq_fft(fq, g, k, a):
    # Caso base: si k = 0, n = 1 => devuelve el propio elemento
    if k == 0:
        return a
    
    n = len(a)

    if n != 2**k:
        raise ValueError("La longitud de a debe ser igual a n = 2**k")
    
    m = n // 2

    # Dividir a en pares e impares
    a_par = a[0::2]
    a_impar = a[1::2]

    # Recursion
    # calcular g^2
    g2 = fq.mult(g, g)

    # Evaluamos en pares e impares
    y0 = fq_fft(fq, g2, k - 1, a_par)
    y1 = fq_fft(fq, g2, k - 1, a_impar)

    # Combinar resultados
    resultado = [fq.cero()] * n
    w = fq.uno()  # w = g^0
    for j in range(m):
        t = fq.mult(w, y1[j]) # t = w * y1[j]
        u = y0[j]
        
        resultado[j] = fq.suma(u, t)    # u + t
        resultado[j + m] = fq.suma(u, fq.inv_adit(t))   # u - t

        w = fq.mult(w, g)  # w = w * g

    return tuple(resultado)


# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fq
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fq_ifft(fq, g, k, a):
    # Inversa de g => g^-1
    inv_g = fq.inv_mult(g)

    # Calcular DFT_{n,g^-1}(a)
    dft = fq_fft(fq, inv_g, k, a)

    # Calcular n^(-1)
    # Como n = 2^k => n^{-1} = (2^k)^{-1} = (2^{-1})^k
    dos = fq.suma(fq.uno(), fq.uno())
    inv_dos = fq.inv_mult(dos)
    inv_n = fq.pot(inv_dos, k)

    # Multiplicar cada elemento por n^(-1)
    resultado = []
    for x in dft:
        resultado.append(fq.mult(x, inv_n))

    return tuple(resultado)