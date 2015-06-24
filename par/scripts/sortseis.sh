echo
echo =====================================================================
echo "Enter number of shots:"
read nshot1
echo

N=1

while test $N -le $nshot1
do
	f="./su/kugel/kugel.shot$N"
	
	echo $f

	susort < $f gy > $f.sort
	cp $f.sort $f
	rm $f.sort
	N=`expr $N + 1`
done
