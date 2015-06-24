% transform ASCII-Data to mat-file


! /opt/CWPsu/bin/sustrip <  ../su/mil1c_z.su | /opt/CWPsu/bin/b2a n1=1 > mil1c_z.asc

% file name:
load mil1c_z.asc;
A=mil1c_z;
clear mil1c_z;

% Number of traces in file:
ntr=144;                       
% Number of samples in file:
ns=400;                     

% samplingrate:
dt=0.0025;                     
% distance between traces:
dx=1.0;

% time window:
% tbeg=0.0;
% tend=0.6999;

% reduce number of samples:
ins=1;
% reduce number of traces:
idx=1;

%---------------------------------------------------------------------------


% sort data into traces
clear sek;
for i=1:ntr,
   sek(1:ns,i)=A(((i-1)*ns+1):i*ns);
end
clear A;

% cut of values not within time window [tbeg, tend]:
%sek=sek(tbeg/dt+1:tend/dt,:);
% ns=tend/dt;

% reduce number of samples:
sek=sek(1:ins:ns,:);
% reduce number of traces:
sek=sek(:,1:idx:ntr);

% modifying number of samples and samplingrate
ns=size(sek,1);
ntr=size(sek,2);

dt=dt*ins;
dx=dx*idx;

save mil1c_z.mat sek ntr ns dt dx
!rm mil1c_z.asc
