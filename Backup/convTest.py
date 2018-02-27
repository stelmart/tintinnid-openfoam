#This is a program to determine whether a run has converged. If it has converged, it will copy the most recent U and p files to the "Update" folder in the main directory and output "true". If it has not converged, the program will update system/controlDict to prepare it for another run and output "false".
#We test for convergence by checking a number of random points to see if there is a difference to within a given tolerance.

import sys
from numpy import *

def updateControlDict(startTime, endTime):
    #This is to rewrite controlDict if simulation has failed to converge
    template = open("../controlDictTemp", "r")
    newFile = open("system/controlDict", "w")

    for i in range(0,21):
        newFile.write(template.readline())

    template.readline()
    newFile.write("startTime       " + startTime + ";\n")

    for i in range(0,3):
        newFile.write(template.readline())

    template.readline()
    newFile.write("endTime         " + endTime + ";\n")

    for line in template:
        newFile.write(line)
    template.close()
    newFile.close()


def copyToUpdate(ult):
    #This is to copy the bulk of the U and p files to the Update directory once each simulation is complete
    newU = open(ult + "/U", "r")
    newp = open(ult + "/p", "r")
    oldU = open("../0/U", "w")
    oldp = open("../0/p", "w")

    for i in range(0,12):
        oldU.write(newU.readline())
        oldp.write(newp.readline())

    newU.readline()
    newp.readline()

    test=1
    while test:
        line1 = newU.readline()
        line2 = newp.readline()
        if "boundaryField" in line1:
            test = 0
        else:
            oldU.write(line1)
            oldp.write(line2)

    newU.close()
    newp.close()
    oldU.close()
    oldp.close()


#Start by finding the two most recent updates
penult = sys.argv[1] 
ult = sys.argv[2]

UPen = open(penult + "/U", "r")
UUlt = open(ult + "/U", "r")

num_points = 100
tolerance = 0.0001 #This is a relative tolerance, not absolute



#Find correct position in files
fline = ""
while "internalField" not in fline:
    fline = UPen.readline()
fline = ""
while "internalField" not in fline:
    fline = UUlt.readline()

#Find total number of blocks
numBlocks = int(UPen.readline())
UUlt.readline()

#Generate random numbers from 0 to numBlocks
intList = sort(random.randint(num_points, size=numBlocks))

#Compare these lines to see if the elements differ to within the given tolerance
j = 0
diffs = 0
for i in range(0,num_points):
    while j <= intList[i]+1:
        penLine = UPen.readline()
        ultLine = UUlt.readline()
        j+=1
    penVec = penLine[1:-2].split()
    ultVec = ultLine[1:-2].split()
    diffx = abs((float(penVec[0]) - float(ultVec[0]))/float(ultVec[0]))
    diffy = abs((float(penVec[1]) - float(ultVec[1]))/float(ultVec[1]))
    diffz = abs((float(penVec[2]) - float(ultVec[2]))/float(ultVec[2]))
    if diffx > tolerance or diffy > tolerance or diffz > tolerance:
        diffs+=1

UPen.close()
UUlt.close()

#If diffs = 0, sim has converged. Otherwise, run it again.
if diffs == 0:
    #copyToUpdate(ult)
    print("true")
else:
    updateControlDict(ult,str(float(ult)+float(sys.argv[3])))
    print("false")
