% shots and receiver
clear all
close all

nr_shots=20;
DG=1; % distance between receiver


shot_x(1)=5.0;
shot_y(1)=0.2;

shot_x(2)=16.2;
shot_y(2)=0.2;

shot_x(3)=26.0;
shot_y(3)=0.2;

shot_x(4)=34.0;
shot_y(4)=0.2;

shot_x(5)=45.2;
shot_y(5)=0.2;

shot_x(6)=55.0;
shot_y(6)=0.2;
 
shot_x(7)=64.6;
shot_y(7)=0.2;

shot_x(8)=75.6;
shot_y(8)=0.2;

shot_x(9)=85.0;
shot_y(9)=0.2;

shot_x(10)=94.0;
shot_y(10)=0.2;

shot_x(11)=38.6;
shot_y(11)=0.2;

shot_x(12)=19.8;
shot_y(12)=0.2;

shot_x(13)=30.6;
shot_y(13)=0.2;

shot_x(14)=40.2;
shot_y(14)=0.2;

shot_x(15)=50.2;
shot_y(15)=0.2;

shot_x(16)=60.8;
shot_y(16)=0.2;
 
shot_x(17)=70.6;
shot_y(17)=0.2;

shot_x(18)=80.6;
shot_y(18)=0.2;

shot_x(19)=48.0;
shot_y(19)=0.2;

shot_x(20)=88.6;
shot_y(20)=0.2;

receiver_x=6:DG:93;
receiver_y=length(receiver_x);
receiver_y(:)=0.2;


figure
hold on
for i=1:nr_shots
    plot(shot_x(i),shot_y(i),'*r','MarkerSize',11)
end
plot(receiver_x,receiver_y,'+k')

save shots_and_receiver