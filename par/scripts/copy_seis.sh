#!/bin/csh
set NT=124

set x=1
while ( $x < 120 )
  cp su/full_wave_forward_x.su.shot10_it$x.* /fastfs/koehn/DENISE_results/crase_model_simple_lam_mu_rho_time_win/seis
  echo "copy seis it ... $x"
  set x = `expr $x + 5`
end
