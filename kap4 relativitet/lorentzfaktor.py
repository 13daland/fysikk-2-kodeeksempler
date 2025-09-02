c = 3*10**8         # Lyshastigheten c [m/s]
m_e = 9.109e-31     # Massen til et elektron m_e [kg]


def gammab(b: float):    # Gamma fra b.
    if b >= 1 or b < 0:
        raise ValueError(f"beta må være mellom 0 og 1. Du har oppgitt {b}")
    else:
        return 1/(1-b**2)**.5


def gammav(v: float, ls: float = c):        # Gamma fra v.
    return gammab(v/ls)


def print_gamma(b: float):           # Printer gamma som en funksjon av beta
    if b >= 1 or b < 0:
        raise ValueError(f"beta må være mellom 0 og 1. Du har oppgitt {b}")
    else:
        print(gammab(b))


def ek(m: float, b: float, ls: float = c):
    g = gammab(b)
    return m * ls ** 2 * (g - 1)

startval = 9999999990
endval = 10000000000

for i in range(startval, endval):
    beta = i/endval
    print(f"v={beta}c, gamma={gammab(beta)}, E_k={ek(75, beta)}")

print(0.9999999999*c-0.9999999998*c)
d1 = ek(75, 0.9999999999)-ek(75, 0.9999999998)
print(d1)
d2 = ek(75, 200/c)-ek(75, 100/c)
print(d2)
print(d1/d2)


