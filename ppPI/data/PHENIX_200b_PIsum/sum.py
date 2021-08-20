
def readfile(filename):
    data = [[] for i in range(5)]
    with open(filename, 'r') as infile:
        i = 0
        for line in infile:
            if i < 5:
                i += 1
                continue
            else:
                items = [float(x) for x in line.strip().split()]
                for j in range(5):
                    data[j].append(items[j+1])
    return data

data_plus  = readfile('../PHENIX_200b_PIp/PHENIX_200b_PIp.txt')
data_minus = readfile('../PHENIX_200b_PIm/PHENIX_200b_PIm.txt')

data = [[data_plus[i][j] + data_minus[i][j] for j in range(len(data_plus[i]))] for i in range(len(data_plus))]
data[0] = data_plus[0]

with open('PHENIX_200b_PIsum.txt', 'w') as outfile:
    outfile.write('PHENIX_200b_PIsum\n')
    outfile.write('%d %d\n' % (len(data[0]), 2))
    outfile.write('200\n')
    outfile.write('-0.35 0.35\n')
    outfile.write('\tpT\tobs\tstat\tsys1\tsys2\n')
    for i in range(len(data[0])):
        outfile.write('%d\t%f\t%e\t%e\t%e\t%e\n' % (i, data[0][i], data[1][i], data[2][i], data[3][i], data[4][i]))
