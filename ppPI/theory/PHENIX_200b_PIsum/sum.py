def readfile(filename):
    data = [[], []]
    with open(filename, 'r') as infile:
        for line in infile:
            items = [float(x) for x in line.strip().split()]
            data[0].append(items[0])
            data[1].append(items[1])
    return data


for i in range(201):
    data_minus = readfile('../PHENIX_200b_PIm/results/1_3_3_1.0_1%03d.res' % i)
    data_plus  = readfile('../PHENIX_200b_PIp/results/1_3_3_1.0_1%03d.res' % i)

    data = [[data_minus[i][j] + data_plus[i][j] for j in range(len(data_plus[i]))] for i in range(len(data_minus))]
    data[0] = data_plus[0]

    with open('results/1_3_3_1.0_1%03d.res' % i, 'w') as outfile:
        for i in range(len(data[0])):
            outfile.write("%e\t%e\n" % (data[0][i], data[1][i]))

