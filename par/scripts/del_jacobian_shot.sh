#!/bin/csh
set NT=100

set x=1
while ( $x < 101 )
  
  mv jacobian/jacobian_p.shot$x jacobian/tmp
  #rm jacobian/jacobian_p.shot$x*
  
  echo "delete shot ... $x"
  set x = `expr $x + 1`
end
