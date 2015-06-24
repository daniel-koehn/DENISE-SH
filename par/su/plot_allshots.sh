#!/bin/csh

set nshots=43
set x=1
while ( $x < $nshots)

supswigp label1="Time [s]"  label2="Trace No." title="Shot no. $x" < Profil_2/P2_TD0p16_y.su.shot$x > test_$x.ps

set x = `expr $x + 1`

end

# merge all shots in one PS-file
cat test_1.ps > allshots.ps 

set x=2
while ( $x < $nshots)

cat test_$x.ps >> allshots.ps

set x = `expr $x + 1`

end

# convert PS->PDF
ps2pdf allshots.ps

# tidy up the place
rm *.ps
