import matplotlib.pyplot as plt 
import numpy as np

A = 10
m = 5
s = 1.5

x = np.arange(m-10, m+10, 0.25)
y = A * np.exp(-(x-m)**2 / (2*(s**2)))

plt.figure()

plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Gaussian curve")

plt.plot(x, y, color="red")
plt.show()