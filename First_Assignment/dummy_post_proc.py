import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss


data = np.genfromtxt('dummy.txt')
x = np.linspace(-5,5,100)
y = ss.t.pdf(x,3.,1.,1.)
plt.plot(x,y,'k--',lw=2)
plt.hist(data,bins=50,normed=True)
plt.show()