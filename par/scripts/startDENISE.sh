#lamboot -v lamhosts  
lamboot
mpirun -np 8 nice -19 ../bin/denise in_and_out/MAMM_profilP2.json | tee in_and_out/MAMM_profilP2.out

#cd su

#for ((i=1; i < ($1+1); i++)) ; do

#	cat DENISE_y.su.shot$i.* > DENISE_y.su.shot$i
#        cat DENISE_x.su.shot$i.* > DENISE_x.su.shot$i

#done

#cd ..			
