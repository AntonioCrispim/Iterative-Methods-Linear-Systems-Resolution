from numpy import*
import matplotlib.pyplot as plt

x = arange(0., 1., 0.01)
y1 = (3*x+1)/4*log2(3*x+1)+3*(1-x)/4*log2(1-x)
y2 = (x+1)/2*log2(x+1)+(1-x)/2*log2(1-x)
y3 = y1-y2
ll1 = plt.plot(x, y1)
ll2 = plt.plot(x, y2)
ll3 = plt.plot(x, y3)
plt.show()
