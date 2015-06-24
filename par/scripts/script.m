% Kind of script to get the start model, final model, comparision of
% measured and final seismograms and misfit history

clear all
close all

% number of iteration
nr_it=83;

nx=500;     % gridpoints x-axis
ny=100;     % gridpoints (vertical) y-axis
DH=0.2;     % gridspacing
bondary=5; % absorbing bondary in gridpoints
P1=20;
P2=80;

% true model v_p, v_s, rho
%model_2D('model_true/model_Test_vs_it_','V_s true model',nx,ny,DH,0,P1,P2)

%model_2D('model_true/model_Test_vp_it_','V_p true model',nx,ny,DH,0,bondary)

%model_2D('model_true/model_Test_rho_it_','Rho true model',nx,ny,DH,0,bondary)


% start model v_p, v_s, rho
model_2D('model/model_Test_vs_it_','V_s start model',nx,ny,DH,0,P1,P2)

%model_2D('model/model_Test_vp_it_','V_p start model',nx,ny,DH,0,bondary)

%model_2D('model/model_Test_rho_it_','Rho start model',nx,ny,DH,0,bondary)


% final model v_p, v_s, rho
model_2D('model/model_Test_vs_it_','V_s final model',nx,ny,DH,nr_it,P1,P2)

%model_2D('model/model_Test_STF_vs_it_','V_p final model',nx,ny,DH,31,P1,P2)

%model_2D('model/model_Test_rho_it_','Rho final model',nx,ny,DH,nr_it,bondary)


%vel_profile('model_true/model_Test_vs_it_0.bin','model/model_Test_vs_it_0.bin','model/model_Test_vs_it_',nx,ny,nr_it,P1,P2,DH)

% true model v_p, v_s, rho
%model_2D('model/model_Test_vs_smoothed_it_','V_s smoothed model',nx,ny,DH,nr_it,bondary)

%model_2D('model/model_Test_vp_smoothed_it_','V_p smoothed model',nx,ny,DH,nr_it,bondary)

%model_2D('model_true/model_Test_rho_smoothed_it_','Rho smoothed model',nx,ny,DH,nr_it,bondary)



% seismograms
%seismograms(['su/DENISE_y.su.shot1.it',num2str(nr_it)],'su/measured_data/DENISE_y.su.shot1')


% misfit history
%misfit('L2_LOG_STF.dat',73)
misfit('L2_LOG.dat',84)