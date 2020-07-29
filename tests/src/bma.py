# An implementation of the Berlekamp Massey Algorithm for the binary field. Taken from Wikipedia.
# https://en.wikipedia.org/wiki/Berlekamp-Massey_algorithm#The_algorithm_for_the_binary_field

# Berlekamp–Massey algorithm for the binary field
def bma(s):
    # 1. Let s be the bits of the stream.
    s = list(map(int, s))
    N = len(s)

    # 2. initialize b and c arrays each of length N to be zeroes
    b, c = [0 for x in s], [0 for x in s]

    # 2. except b_0 <- 1, c_0 <-1
    b[0], c[0] = 1, 1

    # 3. assign L <- 0, m <- -1
    L, m = 0, -1

    # 4. For n=0 step 1 while n<N:
    n = 0
    while n < N:
        # * Let discrepancy d
        d = s[n]
        for i in range(1, L+1):
            d = d ^ (c[i] & s[n-i])

        # * if d=0, then c is already a polynomial which annihilates the portion of the stream from n-L to n.
        # * else
        if (d != 0):
            # - Let t be a copy of c.
            t = c.copy()

            for i in range(0, N-n+m):
                # - Set c_n-m <- (c_n-m ^ b_0), c_n-m+1 <- (c_n-m+1 ^ b_1), ... up to c_N-1 <- (c_N-1 ^ b_N-n-m-1)
                c[n-m+i] = c[n-m+i] ^ b[i]

            # - if L <= n/2, set L <- n + 1, set m <- n, and let b <- t; otherwise leave L, m and b alone.
            if (L <= n*0.5):
                L = n + 1 - L
                m = n
                b = t

        n += 1

    return L