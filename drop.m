%Project Heavy rain
%Group 10
%Isabell Jansson, Ronja Grosz, Anna Georgellis, Caroline Gard

%Simulation of the falling motion for a drop approximated with Euler

n = 50;

Vt = 1.3;
Ttot = 10/Vt;
dt = Ttot/50;
s = 10;

%Add a static wind
windX = 0.1;
windY = 0.2;

[X,Y,Z] = sphere; 

%Three drops are falling separately
for i = 1:3
    r = 20.*rand(2,1)-10;
    
    %Set the initial height for the drop
    s = 10;
    for t = 0:dt:Ttot;
        
        %Set axis for the scene
        axis([-15 15 -15 15 0 10]);
        
        %Redraws the scene
        drawnow
   
        %Approximate the new position with Euler
        s = s - dt*Vt 
        %Add the new position to the sphere
        surf(X+r(1)+windX*s,Y+r(2)+windY*s,Z+s); 
        
        %Axis have to be reset
        axis([-15 15 -15 15 0 10]);
    end

end

close all;