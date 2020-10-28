from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from DosData import *
import pandas as pd
import random
import itertools

TiCoO3 = 'mp-19424'



#GAP: distance between the fermi level and the next nearest non-spaghetti band
#WIDTH: distance between fermi-level and first zero-density energy 

def plot_database(lower = 0, upper = 5):

	energies = np.arange(0, max_bound, precision)
	num_energies = len(energies)

	data = np.zeros((num_energies, num_energies))

	X, Y = np.meshgrid(energies, energies)

	chunksize = 5000

	reader = pd.read_csv('dos.csv', chunksize=chunksize)

	for c, chunk in enumerate(itertools.islice(reader, 10)):
		for r, (index, row) in enumerate(chunk.iterrows()):
			if not r%1000:
				print(f'{c*chunksize + r} entries processed')
			dos = DosData(row)
			width, gap = dos.get_parameters()
			data[round(width / precision), round(gap / precision)] += 1

	fig = plt.figure()
	ax = Axes3D(fig)
	ax.set_xlabel('Gap')
	ax.set_ylabel('Spaghetti band width')
	start = round(lower/precision)
	data[:, :start] = np.zeros((num_energies, start))
	ax.axes.set_xlim(right = 5, left = lower)
	ax.axes.set_ylim(bottom = 5, top = 0)
	ax.plot_surface(X[:, start:], Y[:, start:], data[:, start:])
	plt.show()


if __name__ == "__main__":
		plot_database(lower = .5)