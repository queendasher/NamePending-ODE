# plot_rc.py
import numpy as np
import matplotlib.pyplot as plt

# match your C++ parameters
R = 1.0
C = 0.01
omega = 100.0 * np.pi

# load: time, y0=Uc, y1=t
data = np.loadtxt("output_test_ode.txt", usecols=(0, 1, 2))
t = data[:, 0]
Uc = data[:, 1]

# current I = (U0 - Uc)/R, with U0(t)=cos(omega t)
U0 = np.cos(omega * t)
I = (U0 - Uc) / R

# --- Plot 1: Uc and Ic in the same plot (second y-axis for clarity) ---
fig, ax1 = plt.subplots()
ax1.plot(t, Uc, label=r"$U_C(t)$")
ax1.set_xlabel("time")
ax1.set_ylabel("voltage")
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(t, I, label=r"$I(t)$", linestyle="--")
ax2.set_ylabel("current")

# build a combined legend
lines = ax1.get_lines() + ax2.get_lines()
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc="best")
ax1.set_title("RC circuit: $U_C(t)$ and $I(t)$")

# --- Plot 2: phase-style (dUc/dt vs Uc) ---
dUc_dt = I / C
plt.figure()
plt.plot(Uc, dUc_dt, label=r"$\frac{d}{dt} U_C$ vs $U_C$")
plt.xlabel(r"$U_C$")
plt.ylabel(r"$\frac{d}{dt} U_C$")
plt.title("RC circuit phase plot")
plt.legend()
plt.grid(True)

plt.show()
