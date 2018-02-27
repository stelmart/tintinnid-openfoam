#!/bin/bash 

saveTime="2e-06"

for i in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100; do
	NOTEXISTS="false"
        if [ ! -d "freq$i" ]; then
		NOTEXISTS="true"
        fi
        if $NOTEXISTS; then
        	cp -r Backup/ freq$i
                #cp -r 0/ freq$i
        fi
	cd freq$i/system
        if $NOTEXISTS; then
	        python meshGen.py $i
        fi
	cd ..
        if $NOTEXISTS; then
                blockMesh
                decomposePar
        fi
	notConverged="true"
	while $notConverged; do
	        mpirun -np 12 pimpleFoam -parallel
	        reconstructPar
	        ls -d 1* 2* 3* 4* 5* 6* 7* 8* 9* > timeList.txt
	        sort -g --output=timeList.txt timeList.txt
	        tail -n 2 timeList.txt > finalTs.txt
	        penult=`head -n 1 finalTs.txt`
	        ult=`tail -n 1 finalTs.txt`
	        if `python convTest.py $penult $ult $saveTime`; then
		        notConverged="false"
	        fi
        done
        mail -s Complete$i chico.martin1987@gmail.com </dev/null
        cd ..
done
