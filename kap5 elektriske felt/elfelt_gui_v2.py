import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LogNorm

# =================================================
# FYSISKE KONSTANTER
# =================================================

k = 8.9875517923e9
eps = 1e-12

# =================================================
# NUMERISKE PARAMETRE
# =================================================

N = 250
STEP = 6
SCALE = 50
R_CUT = 0.01
MAX_CHARGES = 8

# =================================================
# GUI-APP
# =================================================

class FieldApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Elektrisk felt og potensial – interaktiv utforsking")
        self.root.geometry("1300x850")

        # Halv bredde av visningsområdet
        self.half_size = 0.2

        # Ladninger
        self.charges = []

        # Midlertidig q når vi skal plassere ny ladning
        self.pending_q = None

        # Mode: 1 = vektorfelt, 2 = |E|-konturer, 3 = potensial + ekvipotensialer
        self.mode = 1

        # Fargeskala lages bare én gang
        self.colorbar_created = False

        # Dragging
        self.dragging_charge_index = None
        self.dragging_pointP = False

        # Punkt P
        self.Px = 0.05
        self.Py = 0.05

        # -------------------------
        # VENSTRE PANEL
        # -------------------------
        control_frame = tk.Frame(root)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        # Område
        tk.Label(control_frame, text="Visningsområde").pack(pady=(5, 0))
        tk.Label(control_frame, text="Halv bredde (m):").pack()
        self.size_entry = tk.Entry(control_frame, width=10)
        self.size_entry.insert(0, str(self.half_size))
        self.size_entry.pack()
        tk.Button(control_frame, text="Oppdater område", command=self.update_size).pack(pady=3)

        # Punkt P
        tk.Label(control_frame, text="Målepunkt P").pack(pady=(10, 0))
        tk.Label(control_frame, text="Px:").pack()
        self.Px_entry = tk.Entry(control_frame, width=10)
        self.Px_entry.insert(0, str(self.Px))
        self.Px_entry.pack()

        tk.Label(control_frame, text="Py:").pack()
        self.Py_entry = tk.Entry(control_frame, width=10)
        self.Py_entry.insert(0, str(self.Py))
        self.Py_entry.pack()

        tk.Button(control_frame, text="Sett P", command=self.set_P_from_entries).pack(pady=3)

        # Felt i P
        tk.Label(control_frame, text="E i punkt P").pack(pady=(10, 0))
        self.Ex_label = tk.Label(control_frame, text="Ex = ")
        self.Ex_label.pack()
        self.Ey_label = tk.Label(control_frame, text="Ey = ")
        self.Ey_label.pack()
        self.E_label = tk.Label(control_frame, text="|E| = ")
        self.E_label.pack()

        # Forhåndsoppsett
        tk.Label(control_frame, text="Forhåndsoppsett").pack(pady=(10, 0))
        self.preset_var = tk.StringVar()
        self.preset_var.set("Tom")
        presets = ["Tom", "Dipol", "To like", "Kvadrupol", "Tilfeldig"]
        tk.OptionMenu(control_frame, self.preset_var, *presets).pack()
        tk.Button(control_frame, text="Last oppsett", command=self.load_preset).pack(pady=3)

        # Ladninger
        tk.Label(control_frame, text="Ladninger (maks 8)").pack(pady=(10, 0))
        self.charge_frame = tk.Frame(control_frame)
        self.charge_frame.pack()

        tk.Label(control_frame, text="Ny ladning q (C):").pack()
        self.new_q_entry = tk.Entry(control_frame)
        self.new_q_entry.insert(0, "1e-9")
        self.new_q_entry.pack()

        tk.Button(control_frame, text="Plasser ny ladning", command=self.start_place_charge).pack(pady=2)
        tk.Button(control_frame, text="Tilfeldige posisjoner", command=self.randomize_positions).pack(pady=4)

        # MODE
        tk.Label(control_frame, text="Visning").pack(pady=(10, 0))
        tk.Button(control_frame, text="MODE 1: Vektorfelt", command=lambda: self.set_mode(1)).pack()
        tk.Button(control_frame, text="MODE 2: |E|-konturer", command=lambda: self.set_mode(2)).pack()
        tk.Button(control_frame, text="MODE 3: Potensial + ekvipotensialer", command=lambda: self.set_mode(3)).pack()

        tk.Button(control_frame, text="Oppdater figur", command=self.update_plot).pack(pady=8)

        # -------------------------
        # HØYRE PANEL (FIGUR)
        # -------------------------
        fig_frame = tk.Frame(root)
        fig_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.fig, self.ax = plt.subplots(figsize=(7, 7))
        self.canvas = FigureCanvasTkAgg(self.fig, master=fig_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.canvas.mpl_connect("button_press_event", self.on_mouse_press)
        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.canvas.mpl_connect("button_release_event", self.on_mouse_release)

        self.update_plot()

    # =================================================
    # GUI-FUNKSJONER
    # =================================================

    def set_mode(self, m):
        self.mode = m
        self.update_plot()

    def update_size(self):
        try:
            val = float(self.size_entry.get())
            if val <= 0:
                raise ValueError
            self.half_size = val
            self.update_plot()
        except ValueError:
            messagebox.showerror("Feil", "Ugyldig størrelse.")

    def set_P_from_entries(self):
        try:
            self.Px = float(self.Px_entry.get())
            self.Py = float(self.Py_entry.get())
            self.update_plot()
        except ValueError:
            messagebox.showerror("Feil", "Ugyldig P-koordinat.")

    def start_place_charge(self):
        if len(self.charges) >= MAX_CHARGES:
            messagebox.showerror("Feil", "Maks 8 ladninger.")
            return
        try:
            q = float(self.new_q_entry.get())
        except ValueError:
            messagebox.showerror("Feil", "Ugyldig q.")
            return
        self.pending_q = q

    # =================================================
    # MUS
    # =================================================

    def on_mouse_press(self, event):
        if event.inaxes != self.ax:
            return

        x = event.xdata
        y = event.ydata

        # P
        if np.hypot(x - self.Px, y - self.Py) < 0.02 * self.half_size:
            self.dragging_pointP = True
            return

        # Ladninger
        active = self.get_active_charges()
        for i, (q, xq, yq) in enumerate(active):
            if np.hypot(x - xq, y - yq) < 0.02 * self.half_size:
                self.dragging_charge_index = i
                return

        # Ny ladning
        if self.pending_q is not None:
            q = self.pending_q
            self.pending_q = None
            self.add_charge(q, x, y)
            return

    def on_mouse_move(self, event):
        if event.inaxes != self.ax:
            return

        x = event.xdata
        y = event.ydata

        if self.dragging_pointP:
            self.Px = x
            self.Py = y
            self.sync_P_entries()
            self.update_plot()
            return

        if self.dragging_charge_index is not None:
            # Finn tilsvarende aktive ladning i full liste
            active_indices = [i for i, c in enumerate(self.charges) if c["active_var"].get() == 1]
            if self.dragging_charge_index < len(active_indices):
                idx = active_indices[self.dragging_charge_index]
                c = self.charges[idx]
                c["x_entry"].delete(0, tk.END)
                c["x_entry"].insert(0, f"{x:.4f}")
                c["y_entry"].delete(0, tk.END)
                c["y_entry"].insert(0, f"{y:.4f}")
                self.update_plot()

    def on_mouse_release(self, event):
        self.dragging_charge_index = None
        self.dragging_pointP = False

    # =================================================
    # LADNINGER
    # =================================================

    def add_charge(self, q, x, y):
        row = tk.Frame(self.charge_frame)
        row.pack(fill=tk.X, pady=2)

        active_var = tk.IntVar(value=1)
        tk.Checkbutton(row, variable=active_var).pack(side=tk.LEFT)

        q_entry = tk.Entry(row, width=10)
        q_entry.insert(0, str(q))
        q_entry.pack(side=tk.LEFT)

        x_entry = tk.Entry(row, width=8)
        x_entry.insert(0, f"{x:.4f}")
        x_entry.pack(side=tk.LEFT)

        y_entry = tk.Entry(row, width=8)
        y_entry.insert(0, f"{y:.4f}")
        y_entry.pack(side=tk.LEFT)

        btn = tk.Button(row, text="Slett")
        btn.pack(side=tk.LEFT)

        charge = {
            "active_var": active_var,
            "q_entry": q_entry,
            "x_entry": x_entry,
            "y_entry": y_entry,
            "frame": row
        }

        btn.config(command=lambda c=charge: self.delete_charge(c))

        self.charges.append(charge)
        self.update_plot()

    def delete_charge(self, charge):
        charge["frame"].destroy()
        self.charges.remove(charge)
        self.update_plot()

    def get_active_charges(self):
        out = []
        for c in self.charges:
            if c["active_var"].get() != 1:
                continue
            try:
                q = float(c["q_entry"].get())
                x = float(c["x_entry"].get())
                y = float(c["y_entry"].get())
                out.append((q, x, y))
            except ValueError:
                pass
        return out

    def randomize_positions(self):
        for c in self.charges:
            x = np.random.uniform(-self.half_size, self.half_size)
            y = np.random.uniform(-self.half_size, self.half_size)
            c["x_entry"].delete(0, tk.END)
            c["x_entry"].insert(0, f"{x:.4f}")
            c["y_entry"].delete(0, tk.END)
            c["y_entry"].insert(0, f"{y:.4f}")
        self.update_plot()

    def clear_charges(self):
        for c in self.charges:
            c["frame"].destroy()
        self.charges = []

    def load_preset(self):
        self.clear_charges()
        L = self.half_size * 0.5
        q = 1e-9

        p = self.preset_var.get()

        if p == "Tom":
            pass

        elif p == "Dipol":
            self.add_charge(+q, -L, 0)
            self.add_charge(-q, +L, 0)

        elif p == "To like":
            self.add_charge(+q, -L, 0)
            self.add_charge(+q, +L, 0)

        elif p == "Kvadrupol":
            self.add_charge(+q, -L, -L)
            self.add_charge(-q, +L, -L)
            self.add_charge(-q, -L, +L)
            self.add_charge(+q, +L, +L)

        elif p == "Tilfeldig":
            n = np.random.randint(2, 5)
            for _ in range(n):
                x = np.random.uniform(-self.half_size, self.half_size)
                y = np.random.uniform(-self.half_size, self.half_size)
                qv = np.random.choice([-1, 1]) * q
                self.add_charge(qv, x, y)

        self.update_plot()

    # =================================================
    # FELT OG POTENSIAL
    # =================================================

    def compute_field_at(self, x, y, charges):
        Ex = 0.0
        Ey = 0.0
        for q, xq, yq in charges:
            Rx = x - xq
            Ry = y - yq
            R2 = Rx*Rx + Ry*Ry + eps
            R = np.sqrt(R2)
            Ex += k * q * Rx / (R2 * R)
            Ey += k * q * Ry / (R2 * R)
        return Ex, Ey, np.sqrt(Ex*Ex + Ey*Ey)

    def compute_potential(self, X, Y, charges):
        V = np.zeros_like(X)
        for q, xq, yq in charges:
            R = np.sqrt((X - xq)**2 + (Y - yq)**2 + eps)
            V += k * q / R
        return V

    def sync_P_entries(self):
        self.Px_entry.delete(0, tk.END)
        self.Px_entry.insert(0, f"{self.Px:.4f}")
        self.Py_entry.delete(0, tk.END)
        self.Py_entry.insert(0, f"{self.Py:.4f}")

    # =================================================
    # PLOTTING
    # =================================================

    def update_plot(self):
        self.ax.clear()

        L = self.half_size
        xmin, xmax = -L, L
        ymin, ymax = -L, L

        charges = self.get_active_charges()

        x = np.linspace(xmin, xmax, N)
        y = np.linspace(ymin, ymax, N)
        X, Y = np.meshgrid(x, y)

        Ex = np.zeros_like(X)
        Ey = np.zeros_like(Y)

        for q, xq, yq in charges:
            Rx = X - xq
            Ry = Y - yq
            R2 = Rx**2 + Ry**2 + eps
            R = np.sqrt(R2)
            Ex += k * q * Rx / (R2 * R)
            Ey += k * q * Ry / (R2 * R)

        E = np.sqrt(Ex**2 + Ey**2)

        if self.mode == 1:
            # Vektorfelt
            Xs = X[::STEP, ::STEP]
            Ys = Y[::STEP, ::STEP]
            Exs = Ex[::STEP, ::STEP]
            Eys = Ey[::STEP, ::STEP]
            Es = E[::STEP, ::STEP]

            mask = np.ones_like(Xs, dtype=bool)
            for q, xq, yq in charges:
                R = np.sqrt((Xs - xq)**2 + (Ys - yq)**2)
                mask &= (R > R_CUT)

            Xp = Xs[mask]
            Yp = Ys[mask]
            Exp = Exs[mask]
            Eyp = Eys[mask]
            Ep = Es[mask]

            if len(Ep) > 0 and np.isfinite(Ep).any():
                Emin = Ep.min()
                Emax = Ep.max()
                norm = LogNorm(vmin=Emin, vmax=Emax) if (Emin > 0 and Emax > Emin) else None

                Exn = Exp / (Ep + eps)
                Eyn = Eyp / (Ep + eps)

                # Flytt pilene slik at de er sentrert i feltpunktet
                Lp = 1.0 / SCALE
                Xstart = Xp - 0.5 * Lp * Exn
                Ystart = Yp - 0.5 * Lp * Eyn

                qv = self.ax.quiver(
                    Xstart, Ystart,
                    Exn, Eyn, Ep,
                    cmap="viridis",
                    norm=norm,
                    scale=SCALE,
                    scale_units="xy",
                    pivot="tail"
                )

                if norm is not None and not self.colorbar_created:
                    self.fig.colorbar(qv, ax=self.ax, label="|E| (log)")
                    self.colorbar_created = True

            self.ax.set_title("Elektrisk felt – vektorfelt")

        elif self.mode == 2:
            # |E|-konturer
            if len(charges) > 0:
                levels = np.logspace(np.log10(E.min() + 1e-3), np.log10(E.max()), 30)
                self.ax.contour(X, Y, E, levels=levels, cmap="viridis")
            self.ax.set_title("|E|-konturlinjer")

        else:
            # Potensial + ekvipotensialer
            if len(charges) > 0:
                V = self.compute_potential(X, Y, charges)
                levels = np.linspace(np.min(V), np.max(V), 25)
                self.ax.contour(X, Y, V, levels=levels, cmap="coolwarm")
            self.ax.set_title("Elektrisk potensial og ekvipotensialer")

        # Tegn ladninger (også inaktive svakt)
        for c in self.charges:
            try:
                q = float(c["q_entry"].get())
                xq = float(c["x_entry"].get())
                yq = float(c["y_entry"].get())
            except ValueError:
                continue

            active = c["active_var"].get() == 1

            if q > 0:
                style = "ro" if active else "r."
            else:
                style = "bo" if active else "b."

            self.ax.plot(xq, yq, style)

        # Tegn punkt P
        self.ax.plot(self.Px, self.Py, "r*", markersize=12)

        # Beregn felt i P
        ExP, EyP, EP = self.compute_field_at(self.Px, self.Py, charges)
        self.Ex_label.config(text=f"Ex = {ExP:.3e} V/m")
        self.Ey_label.config(text=f"Ey = {EyP:.3e} V/m")
        self.E_label.config(text=f"|E| = {EP:.3e} V/m")

        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        self.ax.set_aspect("equal")
        self.ax.set_xlabel("x (m)")
        self.ax.set_ylabel("y (m)")

        self.canvas.draw()

# =================================================
# START
# =================================================

root = tk.Tk()
app = FieldApp(root)
root.mainloop()
