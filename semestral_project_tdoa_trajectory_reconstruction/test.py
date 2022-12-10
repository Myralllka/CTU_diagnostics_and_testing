import numpy as np
import matplotlib.pylab as plt


def getlinear(x, y):

    def inner(x1):
        return m * x1 + b

    top = (len(x) * np.sum(x * y) - np.sum(x) * np.sum(y))
    btm = (len(x) * np.sum(x * x) - np.sum(x) * np.sum(x))
    m = top / btm
    b = (np.sum(y) - m * np.sum(x)) / len(x)
    # m = (len(x) * np.sum(x * y) - np.sum(x) * np.sum(y)) / (len(x) * np.sum(x * x) - np.sum(x) * np.sum(x))
    # b = (np.sum(y) - m * np.sum(x)) / len(x)
    return inner

x = np.array([3, 4])
y = np.array([110, 120])

predict = getlinear(x, y)
plt.scatter(x, y)
plt.scatter(x, y)
plt.scatter(np.array([0, 5.2]), predict(np.array([0, 5.2])))
plt.show()
