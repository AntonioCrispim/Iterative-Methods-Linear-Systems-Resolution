import matplotlib.pyplot as plt
from numpy import*

a = arange(0., 1., 0.001)
y3 = (a**2+0.3*((1-a**2)-2*a*(1-a**2)**(1/2)))/(1+0.3*(1-4*a*(1-a**2)**(1/2)))
y1 = (a**2+0.3*((1-a**2)+2*a*(1-a**2)**(1/2)))/(1+0.3*(1+4*a*(1-a**2)**(1/2)))
y2 = 1-y1
y4 = 1-y3
Ef = -y3*log2(y3)-y4*log2(y4)
Eb = -y1*log2(y1)-y2*log2(y2)
ll2 = plt.plot(a**2, Eb)
ll3 = plt.plot(a**2, Ef)
plt.show()
