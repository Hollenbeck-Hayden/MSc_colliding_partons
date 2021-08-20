import matplotlib.pyplot as plt
from sys import argv

experiment = argv[1]

data = [[], []]
with open(experiment + "/" + experiment + ".txt", 'r') as infile:
    i = 0
    for line in infile:
        if i < 5:
            i += 1
        else:
            items = [float(x) for x in line.strip().split()]
            data[0].append(items[1])
            data[1].append(items[2]*1e6)

plt.loglog(data[0], data[1], '.')
plt.show()
