% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
file_basic='../14/snap_bin/snapaz';

% gridsize and grid spacing 
nx1=250; nx2=965; ny1=1; ny2=100; dh=0.2;
IDX=1; IDY=1;

REFR_START=50;
LINE=4;


% time increment for snapshots:
dtsnap=0.005;

% reduce numer of values:
ix=1; iy=1;

% amplitude over seismogram over xcur traces
% no seismogram visualization when xcur<0
xcur=-5.0;

load  milref_new_0.2.dat
refr=-milref_new_0; clear milref_new_0;

Z1=1.6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ((xcur>0)&(~exist('sek'))),
   su2mat;
   tm=max(sek);
   for i=1:ntr;
     if tm(i), sek(:,i)=sek(:,i)/tm(i); end
   end
   % offset:
   xo=dx:dx:(ntr)*dx;
   % TWT:
   t=dt:dt:ns*dt;
end


x=(nx1:IDX:nx2)*dh-REFR_START;
y=-(ny1-1:IDY:ny2-1)*dh;
rzline=refr(1:IDX:size(refr,1),LINE);
rzline(find(rzline>0.0))=0.0;
rzline(find(abs(rzline)<Z1))=-Z1;
rzline_new=x*0.0;
nrs=abs(x(1)/(IDX*dh));
if nrs==0, nrs=1; end
nre=size(rzline,1)+nrs;
rzline_new(1:nrs-1)=rzline(1);
rzline_new(nrs:nre-1)=rzline;
rzline_new(nre:size(x,2))=rzline(size(rzline,1));
rzline=rzline_new;
i=1;



% How many snapshot files are there ?
   % Creating first name of snapshot file:
   i=1;
   file=[file_basic,int2str(i),'.bin'];
   while (exist(file)==2),
      % Creating next filename:
      i=i+1;
      file=[file_basic,int2str(i),'.bin'];
   end
   lastframe=i-1;
  disp(['There are ' num2str(lastframe) ' snap files']);

colormap(jet);

firstframe=30;
lastframe=30;

% for i=1:lastframe,
 for i=firstframe:1:lastframe,
   % loading data:
     tsnap=i*dtsnap;
     file=[file_basic,int2str(i),'.bin'];
     disp([' loading file ' file]);
     fid=fopen(file,'r');
     vel=fread(fid,[nx2-nx1+1,ny2-ny1+1],'float');
     fclose(fid);
                                %load(file);
     % reduce number of values:
     xn=x(1:ix:size(x,2));
     yn=y(1:iy:size(y,2));
     veln=vel(1:iy:size(vel,1),1:ix:size(vel,2));
     a=max(max(abs(vel)));
     disp([' Maximum amplitude of snapshot: ', num2str(a)]);
      veln=veln/a;

   % plotting in window:
    if xcur>0, subplot(212), end
    surf(xn,yn,veln'); shading interp;
   %  caxis([-1e-9 1e-9]);
     caxis([-0.008 0.008]);
    view(2);
    set(text(110,-18,['T=',num2str(tsnap),' s']),'FontSize',18,'FontWeight','bold');
     xlabel('X[m]')
     ylabel('Depth [m]')
   set(gca,'DataAspectRatio',[2 1 1])
   set(get(gca,'Ylabel'),'FontSize',18);
   set(get(gca,'Ylabel'),'FontWeight','bold');
   set(get(gca,'Xlabel'),'FontSize',18);
   set(get(gca,'Xlabel'),'FontWeight','bold');
   axis([min(xn) max(xn) min(yn) max(yn)])
   hold on
   p_h=plot3(x,rzline,x*0.0+10,'w-');
   set(p_h(:),'LineWidth',2.0);
   set(gca,'FontSize',18);
   set(gca,'FontWeight','bold');
   set(gca,'Linewidth',1.0);
   hold off
   if xcur>0,
     subplot(211),
     sm=max(max(sek));
     axis([(min(xo)-dx) (max(xo)+dx) min(t) max(t)]);
     hold on
     for tr=1:size(sek,2),
        a=sek(1:round(tsnap/dt),tr)*xcur*dx/sm;
        plot(a+xo(tr),t(1:round(tsnap/dt)),'k','Linewidth',2.0);
     end
    % xlabel('X [m]'); 
    ylabel('TWT [s]');
    set(get(gca,'Ylabel'),'FontSize',18);
    set(get(gca,'Ylabel'),'FontWeight','bold');
    set(get(gca,'Xlabel'),'FontSize',18);
    set(get(gca,'Xlabel'),'FontWeight','bold');
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    set(gca,'Linewidth',2.0);
    set(gca,'Box','on');
    set(gca,'xtick',[]);
    pos=get(gca,'Position');
    pos(2)=pos(2)-0.15;
    if xcur>0, set(gca,'Position',pos); end
    hold off
   end
    
    % set(gcf,'Position', [400 383 640 220])
    % set(gca,'Position', [0.13 0.11 0.775 0.4])
    % set(gca,'Units','pixels')
    % set(gca,'Position', [40 1 560 220])
    % set(gcf,'PaperPositionMode','auto')
     pause(0.1)
    % eval(['print -zbuffer -dppmraw ',file_basic,int2str(i),'.ppm']);
     eval(['print -zbuffer -deps snapmil01_',int2str(i),'.ps']);

   % Saving the snapshot:
    %if i==1,  M=moviein(lastframe,gcf); end
    % if i==firstframe,  M=moviein(lastframe-firstframe+1,gcf); end
    % M(:,i)=getframe(gcf);
    % M(:,i-firstframe+1)=getframe(gcf);
 end

% Playing the movie:
%  movie(M,5,12)

% save movie to play later again with play_movie.m
% save mov_mil33a.mat M
