import matplotlib.pyplot as plt
import os

# Sett do_plot til True om du ønsker plot
do_plot = True
# Sett plot_distance til True om du ønsker å tegne avstand som funksjon av tid.
plot_distance = True
# Sett plot_velocity til True om du ønsker å tegne hastighet som funksjon av tid.
plot_velocity = True
# Sett plot_acceleration til True om du ønsker å tegne akselerasjon som funksjon av tid.
plot_acceleration = True

# Sett focus_distance lik False for å undersøke nærmere på aks og fart.
# Sett focus_distance lik True for å undersøke nærmere på avstand.
focus_distance = False

# Sett save_image til True for å lagre grafen du tegner (OBS kan bare gjøres lokalt.)
save_image = True
decide_name = False

# Konstanter
m = 6.64e-27   # massen av alfapartikkelen, kg
q1 = 3.20e-19  # ladningen til alfapartikkelen, C
q2 = 1.26e-17  # ladningen til gullkjernen, C
k = 8.99e9     # konstant i coulombs lov

# Startverdier
r = -1.0e-10    # startavstand fra kjernen, m
v = 1.5e7       # startfart, m/s
t = 0           # starttid, s
min_v = 0      # sluttfart (kan også endres til å stoppe simuleringen basert på andre parametre.)


# Regner ut akselerasjon
def akselerasjon(radius):
    return -k*q1*q2/(m*radius**2)


# Simulering av bevegelsen
dt = 1.0e-23  # tidssteg, s

a_liste = []
r_liste = []
v_liste = []
t_liste = []

while v > min_v:  # gjenta så lenge farten er positiv
    a = akselerasjon(r)     # beregning av akselerasjon, m/s^2
    v = v + a*dt            # beregning av ny fart
    r = r + v*dt            # beregning av ny posisjon
    t = t + dt              # beregning av ny tid
    a_liste.append(a)
    v_liste.append(v)
    r_liste.append(r)
    t_liste.append(t)

print("Alfapartikkelen snur når r =", r, "m")


if do_plot:
    n_plots = sum([1 for a in [plot_distance, plot_velocity, plot_acceleration] if a])
    fig, axs = plt.subplots(n_plots, sharex="col", )
    n_plots_done = 0
    axs[-1].set_xlabel("$t$ / s")

    if plot_distance:   # legger inn graf med avstand som funksjon av tid
        axs[n_plots_done].plot(t_liste, r_liste)
        axs[n_plots_done].set_title("Avstand $r$")
        axs[n_plots_done].set_ylabel("$r$ / m")
        n_plots_done += 1
    else:               # Tvinger fokus på aks + vel om ikke r vises
        focus_distance = False

    if plot_velocity:   # tegner v(t)
        axs[n_plots_done].plot(t_liste, v_liste)
        axs[n_plots_done].set_title("Hastighet $v$")
        axs[n_plots_done].set_ylabel("$v$ / (m/s)")
        n_plots_done += 1

    if plot_acceleration:   # tegner a(t)
        axs[n_plots_done].plot(t_liste, a_liste)
        axs[n_plots_done].set_title("Akselerasjon, $a$")
        axs[n_plots_done].set_ylabel("$a$ / (m/s$^2$)")

    if not focus_distance:  # Snevrer inn grafen så de interessante bitene blir mer åpenbare
        axs[-1].set_xlim(6.64e-18, 6.71e-18)

    fig.tight_layout()      # pynt for å fjerne tomrom mellom koordinatsystemer

    if save_image:          # Lagrer filen på en slik måte at den ikke overskriver andre filer
        list_of_graph_files = [file for file in os.listdir() if (file[:5] == "alpha" and file[-4:] == ".png")]
        print(list_of_graph_files)
        num_files = len(list_of_graph_files)
        new_file_name = f"alpha{num_files}.png"
        print(new_file_name)
        print(new_file_name in list_of_graph_files)
        n_tries = 0
        if decide_name:
            new_file_name =input("Skriv filnavnet du vil bruke > ")+".png"
            while new_file_name in list_of_graph_files:
                new_file_name = input("Det var visst opptatt.\nPrøv på nytt >")+".png"
        else:
            while new_file_name in list_of_graph_files:
                new_file_name = f"alpha{num_files}_{n_tries}.png"
                n_tries += 1

        fig.savefig(new_file_name, bbox_inches="tight", dpi=500)


    fig.show()