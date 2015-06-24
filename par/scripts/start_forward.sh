rm -f su/DENISE* su/measured_data/DENISE*

lamboot
#lamboot -v lamhosts
mpirun -np 8 nice -19 ../bin/denise in_and_out/DENISE_HESSIAN.inp | tee in_and_out/DENISE_HESSIAN.out
cd model
cp model_Test_rho_it_0.bin model_Test_vp_it_0.bin model_Test_vs_it_0.bin ../model_true/.
cd ..

cd su

for ((i=1; i < ($1+1); i++)) ; do

	mv DENISE_y.su.shot$i.it1 measured_data/DENISE_y.su.shot$i
	mv DENISE_x.su.shot$i.it1 measured_data/DENISE_x.su.shot$i

done

cd ..
