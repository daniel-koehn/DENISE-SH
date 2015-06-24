close all
clear all

% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
% Here the basic name of the binary snapshot files must be given:
% (The default extension is *.bin. The increasing number of the snapshot
% is added to the basic name, e.g. if ten snapshots were computed
% and the basic name is snap, then the filenames for the snapshots
% are snap1.bin, snap2.bin, ..., snap10.bin.
filerot='snap/gewoehnlich_1.bin.y';


% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=600;
NY1=1; NY2=200; 
IDX=1; IDY=1;
dh=0.2;

% time increment for snapshots:
TSNAPINC=0.005; TSNAP1=0.005;
FW=10.0;


% seismic traces
ns=6000;
ntr=101;
dt=0.0001;
seis_file='su/gewoehnlich_1_y_qagc_0.2.bin';
dx=1.0; XREC1=10.0; XREC2=110.0;
xcur=1.0;
clip=1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;


firstframe=1;
lastframe=120;
snapno=0;

caxis_value=3.0e-2;

load /home/tbohlen/util/seismic.map
colormap(seismic)



% load seismic traces

 disp([' loading seismic data from file ' seis_file]);
 fid=fopen(seis_file,'r','ieee-le');
 sek=fread(fid,[ns,ntr],'float');
 fclose(fid);
 sm=max(max(sek));
  tm=max(sek);
   for i=1:ntr;
     if tm(i), 
     	sek(:,i)=sek(:,i)/tm(i)/clip; 
	end
   end
   
   
  % clip 
  sek(find(abs(sek)>1.0))=1.0;
  
   
   % offset:
   xo=XREC1:dx:XREC2;
   % TWT:
   t=dt:dt:ns*dt;

% load model
 file='model/gewoehnlich_1.bin.mu';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 model=fread(fid,[ny,nx],'float');
 fclose(fid);
 %model=model(1:2:size(model,1),1:2:size(model,2));
   
 disp(['opening file ' filerot]);
 fid_rot=fopen(filerot,'r','ieee-le');
 

 for i=firstframe:1:lastframe,
 
 hold off
   disp(['loading snapshot no ',int2str(i)]);
   % loading data:
    tsnap=(i-1)*TSNAPINC+TSNAP1;
   offset=4*nx*ny*(i-1);
   fseek(fid_rot,offset,-1);
   velrot=fread(fid_rot,[ny,nx],'float');
   
%   vmp=max(max(abs(veldiv(30:ny,:))));
%   vms=max(max(abs(velrot(30:ny,:))));
   vms=max(max(abs(velrot)));
   disp([' Maximum amplitude of snapshot: ', num2str(vms)]);
  % veldiv=veldiv/vmp;
  % velrot=velrot/vmp;
	


% plotting seismic section
    subplot(211),
     axis([min(x)+FW max(x)-FW min(t) max(t)]);
     hold on
     for tr=1:size(sek,2),
        a=sek(1:round(tsnap/dt),tr)*xcur*dx;
        plot(a+xo(tr),t(1:round(tsnap/dt)),'k','Linewidth',1.0);
     end
    % xlabel('X [m]'); 
    ylabel('Time [s]');
    title(['T=',sprintf('%1.2f',tsnap),' s']);
    set(get(gca,'Title'),'FontSize',12);
    set(get(gca,'Title'),'FontWeight','bold');
   set(get(gca,'Ylabel'),'FontSize',12);
    set(get(gca,'Ylabel'),'FontWeight','bold');    
    set(get(gca,'Xlabel'),'FontSize',12);
    set(get(gca,'Xlabel'),'FontWeight','bold');
    set(gca,'FontSize',12);
    set(gca,'FontWeight','bold');
    set(gca,'Linewidth',2.0);
    set(gca,'Box','on');
    set(gca,'xtick',[]);
    pos=get(gca,'Position');
    pos(2)=pos(2)-0.12;
    if (i==1) set(gca,'Position',pos); end
    hold off
    
    
   % plotting in window:
   % set(gcf,'Position',[100 400 800 550]);
   subplot(212), 
      imagesc(x,y,velrot);  
		 hold on
		contour(x,y,model,2,'k-'); 
	 	%plot(3.0,0.014,'w*','linewidth',4);
 		hold off
    	if (i>0), caxis([-caxis_value caxis_value]); end
				
    
%		set(text(1000,11000,['Max Amp=',sprintf('%1.2g',vms)]),...
%		 'FontSize',14,'FontWeight','bold','color','k');    
%       title('S-waves')
% 		set(text(5000,8000,' x20') ,'FontSize',12,'FontWeight','bold');      
		%axis off
       xlabel('Distance [m]')
       ylabel('Depth [m]')
       set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Linewidth',1.0);
       set(gca,'Box','on');
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij
       %set(gca,'YTick',[400 800 1200])
		 
	

	pause(0.1)
	snapno=snapno+1;
   % Saving the snapshot:
     %epsfile=['eps/ktbsnap',int2str(i),'.eps'];
     %eval(['print -deps ' epsfile]);
    %jpgfile=['ppm/snap_',int2str(i),'.png'];
    %eval(['print -djpeg100 ' jpgfile]);
	pause(10.0)
        tiffile=['ppm/snap_',int2str(i),'.tif'];
        eval(['print -dpng ' tiffile]);
   	pause(10.0)
%ppmfile=['ppm/snap_',int2str(i),'.ppm'];
   % eval(['print -dppmraw ' ppmfile]);
 end
fclose(fid_rot);











