import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("histogram_data.txt")

plt.hist(data, bins=100, density=True)

plt.xlabel("Position (x)")
plt.ylabel("Probability Density |psi(x)|^2")
plt.title("Feynman path integral Histogram")

plt.show()
