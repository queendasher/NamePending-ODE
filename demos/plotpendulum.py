# plot_pendulum_theta.py
import numpy as np
import matplotlib.pyplot as plt

# file produced by the C++ program: columns = [t, theta, theta_dot]
data = np.loadtxt("output_test_ode.txt", usecols=(0, 1, 2))
t = data[:, 0]
theta = data[:, 1]  # angular displacement (radians)

plt.figure()
plt.plot(t, theta, label=r"$\theta(t)$")
plt.xlabel("time [s]")
plt.ylabel("angle [rad]")
plt.title("Pendulum: Angular Displacement vs Time")
plt.grid(True)
plt.legend()
plt.show()
