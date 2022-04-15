clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','docked')
T = 300;
q = 1.602e-19;
kb = 1.38064852*10^-23;
kt = 1.38064852*10^-23*T;
m = (9.10938356*10^-31)*0.26;
vth = sqrt(2*kt/m);
width = 200e-9;
height = 100e-9;
n = 1000;
xp = zeros(n,1);
x = width*rand(n,1);
y = height*rand(n,1);
yp = zeros(n,1);
Nt = 1000;
tmn = 0.2e-12;
dt = (1/100)*(width/vth) ;

%calculate vx and vy
vx = ((rand(n,1)* 2) -1) * vth;
vy = sqrt(vth^2 - vx.^2).*(round(rand(n,1))*2 -1) ;

hold on
axis([0 width 0 height])


for t =1:Nt  
    %position calculations 
    xp = x;
    yp = y;

    x = x + vx*dt;
    y = y + vy*dt;
    
    %particle hitting left and right
    xp(x>width) = 0;
    x(x>width) = x(x>width)-width;
    xp(x<0) = width;
    x(x<0) = width+x(x<0);

    %particle hitting top and bottom
    vy(y>height) = -vy(y>height);
    vy(y<0) = -vy(y<0);
    
      
    v_ave(t) = sum(sqrt((vx.^2)+(vy.^2)))/n;
    temperature(t) =  (((v_ave(t)*m)/2)/kb);
    
    %ploting 
    for i = 1:n
       plot([xp(i),x(i)],[yp(i),y(i)],'Seriesindex',i)
        hold on 
    end
    
     pause(0.1)

end
time = linspace(1,50,Nt);
figure(2) 
plot(time,temperature); 
title('Temperature')
