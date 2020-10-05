import numpy as np
import matplotlib.pyplot as plt
a=np.arange(0.1,1.,0.001)
En=-(a**2)*np.log2(a**2)-(1-a**2)*np.log2(1-a**2)
plt.plot(a**2,En)
plt.show()
