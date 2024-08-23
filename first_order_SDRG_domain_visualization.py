'''
visualize the 3D array outputted by the SDRG program in a 3D plot
the 3D array is a cube so all sides have the same length L 
L is equal to the number of numbers in the first line of the file
the file is a comma-separated list of numbers, with each row separated by a newline

Example:
1,2
3,4
5,6
7,8
this example should be read into a 2x2x2 3D array:
[[[1,2],[3,4]],[[5,6],[7,8]]]

first, read in the data from the file, placing L rows into each 2D matrix of the 3D array

then, plot the 3D array as a 3D grid using different colors for different values of the data
neighboring locations with the same value in the 3D array should have the same color
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
FILE_NAME = 'first_order_SDRG_largest_cluster_1.500000'

# Read data from file
with open(FILE_NAME + '.csv') as file:
    lines = file.readlines()

# Parse data into 3D array
data = []
square = []
for line in lines:
    row = [float(num) for num in line.strip().split(', ')]
    L = len(row)
    square.append(row)
    if len(square) == L:
        data.append(square)
        square = []
data = np.array(data)

# Create a 3D grid
x, y, z = np.meshgrid(range(L), range(L), range(L))

# Plot the 3D array
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
sc = ax.scatter(x, y, z, c=data.flatten(), marker='.',
                cmap=colors.ListedColormap([(0,0,0,0), "darkorange"]),
                norm=colors.CenteredNorm(vcenter=.5), linewidths=0.001)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar(sc)
plt.savefig(FILE_NAME + '.png')
plt.show()
