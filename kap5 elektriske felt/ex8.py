import numpy as np

# Konstanter
q1 = 2.000e-9               # ladning til kule 1, C
p1 = np.array([-0.005, 0])   # posisjon til kule 1, [x, y] m
q2 = -2.000e-9              # ladning til kule 2, C
p2 = np.array([0.005, 0])    # posisjon til kule 2, [x, y] m
k = 8.99e9

# Funksjon som beregner det elektriske feltet i posisjonen r
def e_felt(r):
    r1 = r - p1
    r1_norm = np.linalg.norm(r1)
    r1_enhet = r1 / r1_norm
    e1 = k * q1 / r1_norm ** 2 * r1_enhet

    r2 = r - p2
    r2_norm = np.linalg.norm(r2)
    r2_enhet = r2 / r2_norm
    e2 = k * q2 / r2_norm ** 2 * r2_enhet
    e = e1 + e2
    return e

# Beregning
punkt = np.array([-0.005, 0.01])
e = e_felt(punkt)

# Skriver resultat
print(f"Det elektriske feltet i punktet {punkt} er {e}.")
print(f"Det har en absoluttverdi på {np.linalg.norm(e)} V/m.")
print(f"Vinkelen mellom feltretningen og x-aksen er {np.arctan(e[1]/e[0])} radianer.")
print(np.arctan(e[1]/e[0])*180/3.14159)