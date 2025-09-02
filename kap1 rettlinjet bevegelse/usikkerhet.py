

def gjennomsnitt(observasjoner):
    return sum(observasjoner)/len(observasjoner)


def absolutt_usikkerhet(observasjoner):
    return (max(observasjoner)-min(observasjoner))/2


def relativ_usikkerhet(observasjoner):
    return absolutt_usikkerhet(observasjoner)/gjennomsnitt(observasjoner)


x = [7.25, 7.47, 7.45, 7.26, 7.54, 7.4]

print("a):", gjennomsnitt(x))
print("b):", absolutt_usikkerhet(x))
print("c):", relativ_usikkerhet(x))


def sammensatt_usikkerhet(obs1, obs2, operasjon):
    x_avg = gjennomsnitt(obs1)
    x_abs = absolutt_usikkerhet(obs1)
    x_rel = relativ_usikkerhet(obs1)
    y_avg = gjennomsnitt(obs2)
    y_abs = absolutt_usikkerhet(obs2)
    y_rel = relativ_usikkerhet(obs2)

    if operasjon in ("+", "-"):
        if operasjon == "+":
            print("Du adderer størrelsene.")
            sum_avg = x_avg + y_avg
            print("Gjennomsnitt:", sum_avg)

        else:
            print("Du subtraherer størrelsene.")
            sum_avg = x_avg + y_avg
            print("Gjennomsnitt:", sum_avg)
        sum_abs = y_abs + x_abs
        print(f"Absolutt usikkerhet: {sum_abs:.2}")
        sum_rel = sum_abs / sum_avg
        print(f"Relativ usikkerhet: {sum_rel:.2}")

    elif operasjon in ("*", "/"):
        if operasjon == "*":
            print("Du multipliserer størrelsene.")
            sum_avg = x_avg * y_avg
            print("Gjennomsnitt:", sum_avg)
            minval = min(obs1) * min(obs2)
            maxval = max(obs1) * max(obs2)
        else:
            print("Du dividerer størrelsene.")
            sum_avg = x_avg / y_avg
            print("Gjennomsnitt:", sum_avg)
            minval = min(obs1) / max(obs2)
            maxval = max(obs1) / min(obs2)

        sum_abs = (maxval - minval)/2
        print(f"Absolutt usikkerhet: {sum_abs:.2}")
        sum_rel = x_rel + y_rel
        print(f"Relativ usikkerhet: {sum_rel:.2}")


sammensatt_usikkerhet((3.75, 3.85), (2.15, 2.25), "*")

