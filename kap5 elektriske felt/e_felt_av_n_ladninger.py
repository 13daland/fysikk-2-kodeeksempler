import numpy as np
"""
# Konstanter
ladninger = (2.000e-9, 2.000e-9, -2e-9, -2e-9)
posisjoner = ((-0.05, 0), (0.05, 0), (0, 0.05), (0, -0.05))
"""
k = 8.99e9

# Funksjon som beregner det elektriske feltet i posisjonen r
def e_felt(r, q, p):
    e_liste = []

    if len(q) != len(p):
        print("Forskjellige antall posisjoner og ladninger!")
        return ValueError

    for i in range(len(q)):
        r1 = r - np.array(p[i])
        r1_norm = np.linalg.norm(r1)
        r1_enhet = r1 / r1_norm
        e1 = k * q[i] / r1_norm ** 2 * r1_enhet
        e_liste.append(e1)
    e_ut = np.array([0., 0.])
    for e_komponent in e_liste:
        # print(e_komponent)
        e_ut += e_komponent
    return e_ut

# Beregning

# punkt = np.array([0.1, 0.1])
# e = e_felt(punkt)

# Skriver resultat
# print(f"Det elektriske feltet i punktet {punkt} er {e}.")
# print(f"Det har en absoluttverdi på {np.linalg.norm(e)} V/m.")
# print(f"Vinkelen mellom feltretningen og x-aksen er {np.arctan(e[1] / e[0])} radianer.")