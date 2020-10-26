from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from DosData import *
import random

TiCoO3 = 'mp-19424'



#GAP: distance between the fermi level and the next nearest non-spaghetti band
#WIDTH: distance between fermi-level and first zero-density energy 

def plot_database():

	energies = np.arange(0, max_bound, precision)
	num_energies = len(energies)

	data = np.zeros((num_energies, num_energies))

	X, Y = np.meshgrid(energies, energies)

	with open('database.txt', 'r') as d:
		materials = d.readlines()
		random.shuffle(materials)
		for ID in materials:
			ID = ID.rstrip()
			dos = DosData(ID)
			width, gap = dos.get_parameters()
			data[round(width / precision), round(gap / precision)] += 1 

		fig = plt.figure()
		ax = Axes3D(fig)
		ax.set_xlabel('Gap')
		ax.set_ylabel('Spaghetti band width')

		ax.plot_surface(X, Y, data)
		plt.show()


if __name__ == "__main__":
		plot_database()