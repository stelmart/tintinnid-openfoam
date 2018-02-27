from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

U_vals = np.genfromtxt('U_vals.csv', delimiter=',')

fig = plt.figure()
ax = fig.gca(projection='3d')

#X = U_vals.T[0]
#Y = U_vals.T[1]
#X, Y = np.meshgrid(X, Y)

#print(X, Y)

#Z = U_vals.T[2]

# Plot
for i in range(0,len(U_vals)):
    s = ax.scatter(U_vals[i][0], U_vals[i][1], -U_vals[i][2])

ax.set_xlabel('Oral Radius')
ax.set_ylabel('Cilia Length')
ax.set_zlabel('U')

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
