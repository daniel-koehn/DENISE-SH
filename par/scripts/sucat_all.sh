#!/bin/csh

set x=1
while ( $x < 161)
echo "cat iteration ... $x"
cat su/$1_x.su.shot51_it$x.* > su/$1_x.su.shot51_it$x.tmp
susort < su/$1_x.su.shot51_it$x.tmp > su/$1_x.su.shot51_it$x gelev
rm su/$1_x.su.shot51_it$x.tmp
cp su/$1_x.su.shot51_it$x /fastfs/koehn/DENISE_results/vp_relation_stepl_450/seis
set x = `expr $x + 5` 
end
