#!/bin/csh
set NT=124

set x=1
while ( $x < 23 )
  cp jacobian_p_p.it$x jacobian_constant_eps
  echo "copy model it ... $x"
  set x = `expr $x + 2`
end
