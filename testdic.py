import numpy as np

d = {'a':[], 'b':[]}

names = np.array(['a','b'])
values = np.array([1,2])

def changedic(name, value, d):
    d[name] += [value]
    return d

change = map(changedic, names, values, d))

#for name, value in zip(names, values):
#   changedic(name, value, d)

print(d)

