import sys

dct = sys.argv[1]

T = open(dct + "finalTs.txt")
finalT = T.readline()
finalT = T.readline()[:-1]
U = open(dct + finalT + "/U")
    
for i in range(0,88141):
    q = U.readline()

qlist = U.readline()[1:-2].split()
qz = qlist[2]

params = dct[:-1].split("/")
params[0] = params[0][4:] + ','

with open('U_vals.csv', 'a') as U_vals:
    U_vals.write(params[0] + qz + '\n')

U.close()
T.close()
