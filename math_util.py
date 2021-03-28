def fib(n):
    a, b, c = 1, 1, 0
    for i in bin(n)[3:]:
        t = b * b
        a, b, c = a * a + t, b * (a + c), t + c * c
        if i == "1":
            a, b, c = a + b, a, b
    return b

# hare and tortoise cycle detection
def hat(f, x0):
    t, h = f(x0), f(f(x0))
    while  h != t:
        t, h = f(t), f(f(h))

    t, mu = x0, 0
    while h != t:
        t, h = f(t), f(h)
        mu += 1

    h, la = f(t), 1
    while h != t:
        h = f(h)
        la += 1

    return mu, la

def primesieve(n):
    sieve = [True] * (n + 1)
    i = 3
    while i**2 <= n:
        if sieve[i]:
            j = i**2
            while j <= n:
                sieve[j] = False
                j += 2 * i
        i += 2
    return [2] + [i for i in range(3, n+1, 2) if sieve[i]]

def d_factors(n):
    factors = []
    d = 2
    while d**2 <= n:
        if n % d == 0:
            n //= d
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n != 1:
        factors.append(n)
    return factors

def eulers_totient(n):
    factors = d_factors(n)
    num, den = 1, 1
    for f in factors:
        num, den = (f - 1) * num, f * den
    return n * num // den

def primitive_root(p, prime = True):
    phi = p - 1 if prime else eulers_totient(p)
    factors = d_factors(phi)
    for k in range(2, p + 1):
        primitive = True
        for f in factors:
            if pow(k, phi // f, p) == 1:
                primitive = False
                break
        if primitive:
            return k

class Edge:
    def __init__(self, cost, dest):
        self.cost = cost
        self.dest = dest

def prim(graph):
    import math
    cost = 0
    graphsz = len(graph)
    used = [False] * graphsz
    min_edge = [Edge(math.inf, -1) for i in range(graphsz)]
    min_edge[0].cost = 0
    for _ in graph:
        v = -1
        for i in range(graphsz):
            if not used[i] and (v == -1 or min_edge[i].cost < min_edge[v].cost):
                v = i

        if min_edge[v].cost == math.inf:
            return None

        used[v] = True
        cost += min_edge[v].cost

        for i in range(graphsz):
            if graph[v][i] < min_edge[i].cost and graph[v][i] != -1:
                min_edge[i] = Edge(graph[v][i], v)

    return cost

def nck(n, k):
    if n < k: return 0
    if n == k or k == 0: return 1
    res = k + 1 if k >= n - k else n - k + 1
    i, j = res + 1, 2
    for i in range(res + 1, n + 1):
        res *= i
        res //= j
        j += 1
    return res

def xgcd(a, b):
    x, x1, y, y1 = 0, 1, 1, 0
    while (a != 0):
        q, a, b = b // a, b % a, a
        y, y1 = y1, y - q * y1
        x, x1 = x1, x - q * x1
    return (b, x, y)

def fast_fib(n):
    if n < 2:
        return n
    fn, fnm = 1, 0
    for i in bin(n)[3:]:
        f2nm = fn**2 + fnm**2
        f2n = fn * (2 * fnm + fn)
        fn = f2n + f2nm if int(i) & 1 else f2n
        fnm = f2n if int(i) & 1 else f2nm
    return fn

def divisors(n):
    d = [1]
    for p, e in n:
        d += [div*p**k for k in range(1, e + 1) for div in d]
    return d

def sieve(n):
    primes = primesieve(int(n**0.5))
    factors = [[] for i in range(n + 1)]
    for i in range(2, n + 1):
        if (len(factors[i])):
            continue
        p = i
        e = 1
        while p <= n:
            for j in range(p, n + 1, p):
                if factors[j] and factors[j][-1][0] == i:
                    factors[j][-1] = (i, e)
                else:
                    factors[j].append((i, e))
            e += 1
            p = i**e

    return factors

def squarefree(n):
    sqfree = [1] * (n + 1)
    factors = [[] for i in range(n + 1)]
    for i in range(2, n + 1):
        for j in range(i**2, len(sqfree), i**2):
            sqfree[j] = 0
        if not len(factors[i]):
            for k in range(i, n + 1, i):
                factors[k].append(i)
    return [factors[i] for i in range(2, len(factors)) if len(factors[i]) and sqfree[i]]

def mu_sieve(n):
    import numpy as np
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    primes = primesieve(n)
    for p in primes:
        mu[p::p] *= -1
        mu[p**2::p**2] = 0
    return mu
