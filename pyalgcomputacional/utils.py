def factorizar_entero(n):
    fact = []
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    if e != 0:
        fact.append((2, e))

    p = 3
    while p * p <= n:       # prueba con todos los primos p <= sqrt(n)
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e != 0:
            fact.append((p, e))
        p += 2              # Solo impares
    if n > 1:               # si queda un factor primo > sqrt(n)
        fact.append((n, 1))     
    return fact