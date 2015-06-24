#!/bin/csh
set NT=124

set x=1
while ( $x < 351 )
  rm jacobian/jacobian_p_p.it$x.*
  rm jacobian/jacobian_p_g.it$x.*
  rm jacobian/jacobian_p_c.it$x.*
  
  rm jacobian/jacobian_p_p_rho.it$x.*
  rm jacobian/jacobian_p_g_rho.it$x.*
  rm jacobian/jacobian_p_c_rho.it$x.*
  
  
  rm jacobian/jacobian_p_p_u.it$x.*
  rm jacobian/jacobian_p_g_u.it$x.*
  rm jacobian/jacobian_p_c_u.it$x.*

  
  rm model/waveform_test_model_rho_it_$x.bin.*
  rm model/waveform_test_model_vp_it_$x.bin.*
  rm model/waveform_test_model_vs_it_$x.bin.*
  
  rm model/waveform_test_model_rho_it_$x.bin
  rm model/waveform_test_model_vp_it_$x.bin
  rm model/waveform_test_model_vs_it_$x.bin
  
  echo "delete model it ... $x"
  set x = `expr $x + 1`
end
