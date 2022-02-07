clearvars; close all
set(0,'DefaultFigureWindowStyle','docked')
T = 300;
kb = 1.38064852*10^-23;
kt = 1.38064852*10^-23*T;
m = (9.10938356*10^-31)*0.26;
vth = sqrt(2*kt/m);
width = 200e-9;
height = 100e-9;
n = 100;
x = rand(n,1)*width;
xp = zeros(n,1);
y = rand(n,1)*height;
yp = zeros(n,1);
Nt = 100;
tmn = 0.2e-12;
dt = (1/20)*(width/vth) ;
sample = randi(n,10,1);
vx =  ((rand(n,1)* 2) -1) * vth;
vy = sqrt(vth^2 - vx.^2).*(round(rand(n,1))*2 -1) ;
%maxwell-boltzmann distribution
nx = ((m/(2*kb*pi*T))^(1/2))*exp(-(m*vx.^2)/(2*kb*T));
ny = (m/(2*kb*pi*T))^(1/2)*exp(-(m*vy.^2)/(2*kb*T));
randomspeedx =randn(100,1) *sqrt(kt/m)+vth;
randomspeedy = randn(100,1)*sqrt(kt/m)+vth;

boxes = {};
boxes{1}.X = [0.8 0.4]*1e-7;
boxes{1}.Y = [0.6 1]*1e-7;
boxes{2}.X = [0.8 0.4]*1e-7;
boxes{2}.Y = [0 0.4]*1e-7;
% boxes{3}.X = [1.4 1]*1e-7;
% boxes{3}.Y = [0.4 0.2]*1e-7;







% retangle(xmin ymin howwide howlong)
for z = 1:2
hold on
rectangle('position',[boxes{z}.X(1), boxes{z}.Y(1),boxes{z}.X(2),boxes{z}.Y(2)])
axis([0 width 0 height])
end



for t =1:Nt  
    
    %scattering probabilty
    Pscat = 1-exp(-dt/tmn);
    if Pscat >rand()
        vx = randomspeedx;
        vy = randomspeedy;
    end 
    
    %position calculations 
    xp = x;
    yp = y;
    x = x + (randomspeedx)*dt ;
    y = y + (randomspeedy)*dt ;
    
    % Boundary conditions 
    xp(x>width) = 0;
    x(x>width) = x(x>width)-width;
    xp(x<0) = width;
    x(x<0) = width+x(x<0);
    %vy(y>height) = -vy(y>height);
    %vy(y<0) = -vy(y<0);
    randomspeedy(y>height) = -randomspeedy(y>height);
    randomspeedy(y<0) = -randomspeedy(y<0);
    

    %Vmean = mean(vx.^2 + vy.^2);
    Velo = vx.^2 + vy.^2;
    temperature = (Velo.^2.*m)/(2*kb);
      
    %ploting 
    for i = 1:length(sample)
        plot([xp(sample(i)),x(sample(i))],[yp(sample(i)),y(sample(i))],'Seriesindex',i)
        title(['Temperature = ' num2str(temperature(i))])
        hold on 
    end
    
    pause(0.1)
end


