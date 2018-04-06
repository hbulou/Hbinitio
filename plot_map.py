#!/usr/bin/python
"""Simple matshow() example."""
import matplotlib.pyplot as plt
import numpy as np

matrix = np.loadtxt('wfc3.mat', usecols=range(40))
print(matrix)

## Display matrix
plt.matshow(matrix)

plt.show()
