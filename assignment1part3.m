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
x = ((0.1e-7)-0.05e-7).*rand(n,1)+0.05e-7;
y = ((0.55e-7)-0.45e-7).*rand(n,1)+0.45e-7;
yp = zeros(n,1);
Nt = 1000;
tmn = 0.2e-12;
dt = (1/100)*(width/vth) ;

%maxwell-boltzmann distribution
vx = randn(n,1) *vth/sqrt(2);
vy = randn(n,1) *vth/sqrt(2);
vx_thermalized = randn(n,1) *vth/sqrt(2);
vy_thermalized = randn(n,1) *vth/sqrt(2);

boxes = {};
boxes{1}.X = [0.8 0.4]*1e-7;
boxes{1}.Y = [0.6 0.4]*1e-7;
boxes{2}.X = [0.8 0.4]*1e-7;
boxes{2}.Y = [0 0.4]*1e-7;
boxes{3}.X = [1.4 0.3]*1e-7;
boxes{3}.Y = [0.4 0.2]*1e-7;

%retangle(xmin ymin howwide howlong)
for z = 1:3
hold on
rectangle('position',[boxes{z}.X(1), boxes{z}.Y(1),boxes{z}.X(2),boxes{z}.Y(2)])
axis([0 width 0 height])
end

for t =1:Nt  
    %position calculations 
    xp = x;
    yp = y;

    x = x + vx*dt;
    y = y + vy*dt;

    %scattering probabilty
    Pscat = 1-exp(-dt/tmn);
    scat = Pscat > rand(n,1);
    vx(scat) = vx_thermalized(scat);
    vy(scat) = vy_thermalized(scat);
      
    %Boundary conditions 
    %inside top and bottom box and right box
    InTopBox = x < 1.2e-7 & x > 0.8e-7 & y > 0.6e-7;
    InBotBox = x < 1.2e-7 & x > 0.8e-7 & y < 0.4e-7;
    InrightBox = x < 1.7e-7 & x > 1.4e-7 & y < 0.6e-7 & y > 0.4e-7;
    
    %right and left of the boxes
    reflectxtopright = xp > 1.2e-7 & yp > 0.6e-7;
    reflextxtopleft = xp < 0.8e-7 & yp > 0.6e-7;
    reflectxbotright = xp > 1.2e-7 & yp < 0.4e-7;
    reflextxbotleft = xp < 0.8e-7 & yp < 0.4e-7;
    reflectmiddle = xp < 1.2e-7 & xp > 0.8e-7;

    %reverse xspeed for topbox 
    vx(reflectxtopright & InTopBox) = -vx(reflectxtopright & InTopBox);
    vx(reflextxtopleft & InTopBox) = -vx(reflextxtopleft & InTopBox);
   
    %reverse xspeed for bottom box
    vx(reflectxbotright & InBotBox) = -vx(reflectxbotright & InBotBox);
    vx(reflextxbotleft & InBotBox) = -vx(reflextxbotleft & InBotBox);
    
    %reverse x for bottom box
    x(reflectxbotright & InBotBox) = xp(reflectxbotright & InBotBox);
    x(reflextxbotleft & InBotBox) = xp(reflextxbotleft & InBotBox);
    
    %reverse x for top box
    x(reflectxtopright & InTopBox) = xp(reflectxtopright & InTopBox);
    x(reflextxtopleft & InTopBox) = xp(reflextxtopleft & InTopBox);
    
    %reverse y for top and bottom box
    y(reflectxbotright & InBotBox) = yp(reflectxbotright & InBotBox);
    y(reflextxbotleft & InBotBox) = yp(reflextxbotleft & InBotBox);
    y(reflectxtopright & InTopBox) = yp(reflectxtopright & InTopBox);
    y(reflextxtopleft & InTopBox) = yp(reflextxtopleft & InTopBox);
    
    %reverse yspeed and y for middle 
    vy(reflectmiddle & (InTopBox|InBotBox)) = -vy(reflectmiddle & (InTopBox|InBotBox));
    y(reflectmiddle & (InTopBox|InBotBox)) = yp(reflectmiddle & (InTopBox|InBotBox));

    %pparticle hitting left and right
    vx(x>width) = -vx(x>width);
    vx(x<0) = -vx(x<0);

    %particle hitting top and bottom
    vy(y>height) = -vy(y>height);
    vy(y<0) = -vy(y<0);
    
    
    %reflect top and bottom of right box
    toprightBox = xp < 1.7e-7 & xp > 1.4e-7 & yp > 0.6e-7;
    botrightBox = xp < 1.7e-7 & xp > 1.4e-7 & yp < 0.4e-7;
    
    %reflect right and left of right box
    rightrightBox =  xp < 1.4e-7 & yp < 0.6e-7 & yp > 0.4e-7;
    leftrightBox = xp > 1.7e-7  & yp < 0.6e-7 & yp > 0.4e-7;
    
    %reverse x and y speed for right box
    vy(toprightBox& InrightBox) = -vy(toprightBox & InrightBox);
    vy(botrightBox & InrightBox) = -vy(botrightBox & InrightBox);
    vx(rightrightBox& InrightBox) = -vx(rightrightBox & InrightBox);
    vx(leftrightBox & InrightBox) = -vx(leftrightBox & InrightBox);
  
    v_ave(t) = sum(sqrt((vx.^2)+(vy.^2)))/n;
 
    %ploting 
    for i = 1:n
       plot([xp(i),x(i)],[yp(i),y(i)],'Seriesindex',i)
        hold on 
    end
    
     pause(0.1)
end
 

Nbinx = 40;
Nbiny = 30;
binx = ceil(x*Nbinx/width);   % 100nm -> bin 5
% ->     [ 1, 4, 7, 7 .. bin_of_last_part ]
biny = ceil(y*Nbiny/height);

for  i=1:Nbinx
    for j = 1:Nbiny
    match = binx==i & biny==j;
    ebox(i,j) = sum(match);
    end
end

for  i=1:Nbinx
    for j = 1:Nbiny
    match = binx==i & biny==j;
    temperature(i,j) = sum(((sqrt(((vx(match)).^2)+((vy(match)).^2))*m)/2)/kb);
    end
end

z = 1:1:Nt;
figure(2)
surf(ebox)
title('Density')
figure(3)
surf(temperature)
title('Temperature')