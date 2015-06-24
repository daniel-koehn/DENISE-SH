#!/bin/csh
set NT=124

set x=4
while ( $x < 89)
  cp model/waveform_test_model_vs_it_$x.bin model/tmp/waveform_test_model_vs_it_$x.bin
  cp model/waveform_test_model_vp_it_$x.bin model/tmp/waveform_test_model_vp_it_$x.bin
  cp model/waveform_test_model_rho_it_$x.bin model/tmp/waveform_test_model_rho_it_$x.bin
  echo "copy model it ... $x"
  set x = `expr $x + 4`
end
