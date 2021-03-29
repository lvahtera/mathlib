# mathlib
field.h implements [Galois fields](https://en.wikipedia.org/wiki/Finite_field)  

misc_al_t.h implements several common algorithms that have no implementation in STL.  

mod_a_t.h implements several common functions used in [modular arithmetic](https://en.wikipedia.org/wiki/Modular_arithmetic) that have a non-trivial implementation.  

primitive_root.h is used to solve a mathematical problem of the [same name](https://en.wikipedia.org/wiki/Primitive_root_modulo_n#Finding_primitive_roots) that has its uses in e.g. cryptography  

tonellishanks.h implements the [Tonelli-Shanks algorithm](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm)  

primes_t.h contains several functions related to primes in one way or another:
  - Two differing implementations of the Sieve of Eratosthenes, which can be used to find all prime numbers up to a given limit.  
  - Several algorithms that can be used to test whether a number is a prime, including the [Miller-Rabin primality test](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test)  
  - Number factorisation  
  - [Euler's totient function](https://en.wikipedia.org/wiki/Euler%27s_totient_function)  
    
primecount.h can be used to rapidly count the number of prime numbers under a given limit without generating all of them based on the [Meissel-Lehmer algorithm](https://en.wikipedia.org/wiki/Meissel%E2%80%93Lehmer_algorithm)

math_util.py contains similar algorithms or simplified versions implemented in Python 3
