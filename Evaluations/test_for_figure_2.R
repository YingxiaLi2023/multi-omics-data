# P-value:
1 - phyper(q=27-1, m=16*5, n=(31-16)*5, k=30)


# Adjusted p-value:
5*(1 - phyper(q=27-1, m=16*5, n=(31-16)*5, k=30))


# Expected number of combinations that feature a particular block
# under the null hyptohesis:

30*16*5/(31*5)
