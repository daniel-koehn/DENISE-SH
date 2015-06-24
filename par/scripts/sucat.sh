#!/bin/csh

set x=1
while ( $x < 51)

#cat su/$1_y.su.? > su/$1_y.su
cat su/full_wave_forward_y.su.shot$x.* > su/tmp.su
#cat su/$1_p.su.? > su/$1_p.su
#cat su/$1_div.su.? > su/$1_div.su
#cat su/$1_rot.su.? > su/$1_rot.su

susort < su/tmp.su gx > su/tmp/full_wave_forward_y.su.shot$x 

set x = `expr $x + 1`

end
rm su/tmp.su