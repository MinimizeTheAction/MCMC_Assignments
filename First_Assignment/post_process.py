#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
import sys

def main(data_file, figure_file):
	
	data = np.genfromtxt(data_file, delimiter=" ")
	samples = data[:,0]
	accept = data[:,1]

	acceptance_rate = float(sum(accept))/len(accept)*100.

	plt.figure(1)
	plt.plot(samples)
	plt.savefig(figure_file)
	
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
