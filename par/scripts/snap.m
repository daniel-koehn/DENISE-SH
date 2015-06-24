%close all
%clear all

% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
% Here the basic name of the binary snapshot files must be given:
% (The default extension is *.bin. The increasing number of the snapshot
% is added to the basic name, e.g. if ten snapshots were computed
% and the basic name is snap, then the filenames for the snapshots
% are snap1.bin, snap2.bin, ..., snap10.bin.
filerot='snap/ktb10.bin.rot';
filediv='snap/ktb10.bin.div';


% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=400;
NY1=1; NY2=600; 
IDX=1; IDY=1;
dh=0.005;

% time increment for snapshots:
TSNAPINC=0.01; TSNAP1=0.01;
FW=0.100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;


firstframe=1;
lastframe=50;
snapno=0;

caxis_value=1.0e-4;

%load /home/tbohlen/util/seismic.map
colormap(gray)

% load model
 file='model/ktb10.pi';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 model=fread(fid,[ny,nx],'float');
 fclose(fid);
 %model=model(1:2:size(model,1),1:2:size(model,2));
   
 disp(['opening file ' filerot]);
 fid_rot=fopen(filerot,'r','ieee-le');
 
 disp(['opening file ' filediv]);
 fid_div=fopen(filediv,'r','ieee-le');

 for i=firstframe:1:lastframe,
 
 hold off
   disp(['loading snapshot no ',int2str(i)]);
   % loading data:
    tsnap=(i-1)*TSNAPINC+TSNAP1;
   offset=4*nx*ny*(i-1);
   fseek(fid_div,offset,-1);
   fseek(fid_rot,offset,-1);
   veldiv=fread(fid_div,[ny,nx],'float');
   velrot=fread(fid_rot,[ny,nx],'float');
   
%   vmp=max(max(abs(veldiv(30:ny,:))));
%   vms=max(max(abs(velrot(30:ny,:))));
   vmp=max(max(abs(veldiv)));
   vms=max(max(abs(velrot)));
   disp([' Maximum amplitude of P-snapshots: ', num2str(vmp)]);
   disp([' Maximum amplitude of S-snapshots: ', num2str(vms)]);
  % veldiv=veldiv/vmp;
  % velrot=velrot/vmp;
	

   % plotting in window:
   % set(gcf,'Position',[100 400 800 550]);
   subplot(121), 
      imagesc(x,y,velrot);  
		 hold on
		contour(x,y,model,2,'w-'); 
	 	%plot(3.0,0.014,'w*','linewidth',4);
		 xvsp=[1.3, 1.3]; yvsp=[0.002,2.9]; 
		 plot(xvsp,yvsp,'w:','linewidth',1);
 		hold off
    	if (i>0), caxis([-caxis_value caxis_value]); end
				
		set(text(1.4,-0.0,['T=',sprintf('%1.2f',tsnap),' s']),...
		 'FontSize',12,'FontWeight','bold','color','k');
    
%		set(text(1000,11000,['Max Amp=',sprintf('%1.2g',vms)]),...
%		 'FontSize',14,'FontWeight','bold','color','k');    
       title('S-waves')
% 		set(text(5000,8000,' x20') ,'FontSize',12,'FontWeight','bold');      
		%axis off
       xlabel('Distance [km]')
       ylabel('Depth [km]')
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
		 
	
   subplot(122), 
	set(gcf,'Color',[1 1 1])
       %surf(x,y,veldiv); shading interp; 
       imagesc(x,y,veldiv);  
		 hold on
		%contour(x,y,model,'w-'); 
		 contour(x,y,model,2,'w-'); 
	 	%plot(3.0,0.014,'w*','linewidth',4);
		 xvsp=[1.3, 1.3]; yvsp=[0.002,2.9]; 
		 plot(xvsp,yvsp,'w:','linewidth',1);
 		hold off
		 if (i>0), caxis([-caxis_value caxis_value]); end

%		set(text(1000,11000,['Max Amp=',sprintf('%1.2g',vmp)]),...
%		 'FontSize',12,'FontWeight','bold');    
       title('P-waves')
% 		set(text(5000,8000,' x20') ,'FontSize',12,'FontWeight','bold');      
       %axis off
       xlabel('Distance [km]')
       %ylabel('Depth [m]')
       set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij

        %set(gca,'YTick',[400 800 1200])
 

   	pause(0.1);
	%set(gcf,'Renderer','zbuffer')
	%M(i-firstframe+1)=getframe(gcf);

	%brighten(0.5)
	

	snapno=snapno+1;
   % Saving the snapshot:
     %epsfile=['eps/ktbsnap',int2str(i),'.eps'];
     %eval(['print -deps ' epsfile]);
    jpgfile=['eps/ktb11_',int2str(i),'.jpg'];
    eval(['print -djpeg100 ' jpgfile]);
%    ppmfile=['eps/ktbsnap',int2str(i),'.ppm'];
%    eval(['print -dppmraw ' ppmfile]);
 end
fclose(fid_div);
fclose(fid_rot);











