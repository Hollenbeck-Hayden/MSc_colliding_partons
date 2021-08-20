import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [8,10]

from math import log


data = []
with open("output", 'r') as infile:
    for line in infile:
        items = line.strip().split()
        data.append((float(items[1]), float(items[2])))

sdata = sorted(data, key=lambda x: x[1])

xs = [i+1 for i in range(len(data))]
ws = [x[0] for x in sdata]
cs = [log(x[1]) for x in sdata]

plt.plot(xs, ws, '-', label="weight")
plt.plot(xs, cs, '-', label=r"$\log(\chi^2)$")
plt.xlabel('replica #')
plt.ylabel(r'weight and $\log(\chi^2)$')
plt.title(r"Sorted graph of (log) $\chi^2$ and corresponding weights")
plt.legend()
plt.savefig('weight_dist.png')
