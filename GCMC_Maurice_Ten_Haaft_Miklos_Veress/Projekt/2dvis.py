#!@note courtesy of Yuding Ai
#author M. Klopotek
#2016-Jan
#Ein Python-Skript zur Visualisierung 2D Staebchen auf dem Gitter

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')


a = 0
with open("Senkrechte.dat", "r") as file:
    for line in file:
        a= a+1
xpos = np.zeros(a)
ypos = np.zeros(a)


i = 0
with open("Senkrechte.dat", "r") as file:
    for line in file:
        words = line.split()
        wx = words[0]
        wy = words[1]
        xpos[i] = wx
        ypos[i] = wy
        i = i+1

dx = np.ones(a)
dy = np.ones(a)

for y in range(0,a):
    dy[y] = 8 #Aspektverhaeltnis
    if a != 0:
        ax.add_patch(
            patches.Rectangle(
                (xpos[y], ypos[y]),
                dx[y],
                dy[y],
                facecolor="red",
                linewidth=0.3
            )
        )



a = 0
with open("Waagerechte.dat", "r") as file:
    for line in file:
        a= a+1
xpos = np.zeros(a)
ypos = np.zeros(a)

i = 0
with open("Waagerechte.dat", "r") as file:
    for line in file:
        words = line.split()
        wx = words[0]
        wy = words[1]
        xpos[i] = wx
        ypos[i] = wy
        i = i+1

dx = np.ones(a)
dy = np.ones(a)

for y in range(0,a):
    dx[y] = 8 #Aspektverhaeltnis
    if a != 0:
        ax.add_patch(
            patches.Rectangle(
                (xpos[y], ypos[y]),
                dx[y],
                dy[y],
                facecolor="blue",
                linewidth=0.3
            )

        )
        
plt.axis('equal')
plt.grid()
fig.savefig('Visual2DStaebchen.png', dpi=300, bbox_inches='tight')
# plt.show()
exit()