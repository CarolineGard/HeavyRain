clc
% For recording
aviobj = VideoWriter('example.avi');
aviobj.FrameRate = 60;
open(aviobj);
tjohej= 0.1;

g = 9.81; %Gravity
dx = 1.0; %x movement  
dy = 1.0; %y movement
dt = 0.05; % time step

H = ones(41:41);  %Height of plane 
U = zeros(41:41); %Velocity in x-direction
V = zeros(41:41); %Velocity in y-direction (negative)

%intended for the calculations in first step, Lax Wendroff
Hx = zeros(40:40); 
Ux = zeros(40:40); 
Vx = zeros(40:40); 

Hy = zeros(40:40);
Uy = zeros(40:40);
Vy = zeros(40:40);

%Color of the plane, dependent on velocity
C = zeros(40:40);

%Create initial displacement
[x, y] = meshgrid( linspace(-4, 4, 10) );
R = sqrt(x.^2 + y.^2) + eps;
Z = (sin(R)./R);
Z = max(Z,0); 

%Create displacement to the height matrix
w = size(Z,1);
t = 10:w+9;
l = 20:w+19;
H(t,l) = H(t, l) + Z;

%Grids
clf 
shg
n = 40;
grid = surf(H);
axis([0 41 0 41 -1 3]);
hold all;

% To restore like it was before, set while 1 == 1
 while tjohej < 10   
     % Fixed timestep so that the animation will end
     tjohej = tjohej + 0.01;
     
     % Reflective boundary conditions 
     % for height H and velocities U and V
     H(:,1) = H(:,2);             
     H(:,n+2) = H(:,n+1);     
     H(1,:) = H(2,:);           
     H(n+2,:) = H(n+1,:);    
    
     U(:,1) = U(:,2);
     U(:,n+2) = U(:,n+1);
     U(1,:) = -U(2,:);
     U(n+2,:) = -U(n+1,:);
     
     V(:,1) = -V(:,2);
     V(:,n+2) = -V(:,n+1);
     V(1,:) = V(2,:);
     V(n+2,:) = V(n+1,:);
     
    %Lax- Wendroff Step method
    
    %FIRST STEP

    %X-direction
    i = 1:38+1;
    j = 1:38;

    Hx(i, j) = ((H(i+1, j+1) + H(i, j+1))/2) - (dt/(2*dx))*(U(i+1, j+1) - U(i, j+1));%.*H(i+1,j+1)-.*H(i, j+1));

    Ux(i, j) = ((0.5*(U(i+1,j+1) + U(i,j+1)) ...
                - (dt/(2*dx))*(U(i+1,j+1).^2  ... 
                + ((0.5*g*H(i+1,j+1) - (U(i,j+1).^2 ... 
                + 0.5*g*H(i,j+1)))))));
            
    Vx(i, j) = 0.5*(V(i,j+1) + V(i+1,j+1)) ...
                - (dt/2*dx)*((U(i+1,j+1).*V(i+1,j+1) ...
                - U(i,j+1).*V(i,j+1)));    

    %Y-direction
    i = 1:38;
    j = 1:38+1;
    
    Hy(i, j) = 0.5*(H(i+1, j+1) + H(i+1, j)) - (dt/(2*dy))*(V(i+1, j+1) - V(i+1, j));%.*H(i+1,j+1)- .*H(i+1, j));
            
    Uy(i, j) = 0.5*(U(i+1,j) + U(i+1,j+1)) ...
               - (dt/2*dy)*((V(i+1,j+1).*U(i+1,j+1) ...
               - V(i+1,j).*U(i+1,j)));
    
    Vy(i, j) = 0.5*(V(i+1,j+1) + V(i+1,j)) ...
               - (dt/(2*dy))*(V(i+1,j+1).^2 ...
               + 0.5*g*H(i+1,j+1) - (V(i+1,j).^2 ...
               + 0.5*g*H(i+1,j)));

       
    %SECOND STEP
    i = 2:38+1;
    j = 2:38;
    
    %Height of the plane 
    H(i, j) = H(i, j) - (dt/dx)*(Ux(i,j-1) - Ux(i-1, j-1)) - (dt/dy)*(Vy(i-1,j) - Vy(i-1,j-1));

    %Velocity in x-direction
    U(i, j) = U(i, j) - (dt/dx)*(Ux(i,j-1).^2 ... 
                      + (0.5*g*Hx(i,j-1) - (Ux(i-1,j-1).^2 + 0.5*g*Hx(i-1,j-1) ))) ...
                      - (dt/dy)*(Vy(i-1,j).*Uy(i-1,j) - Vy(i-1,j-1).*Uy(i-1,j-1)); 
    
    %Velocity in y-direction
    V(i, j) = V(i, j) - (dt/dy)*((Vy(i-1,j).^2 ... 
                      + 0.5*g*Hy(i-1,j)) - (Vy(i-1,j-1).^2 + 0.5*g*Hy(i-1,j-1))) ...
                      - (dt/dx)*(Ux(i,j-1).*Vx(i,j-1) - Ux(i-1, j-1).*Vx(i-1,j-1)); 
    
     %draw the plane and the colour in each point
     C = abs(U(i,j)) + abs(V(i,j)); %Calculate the colors for all points in the plane
     set(grid, 'zdata', H(i,j), 'cdata', C);
     axis off;
     
     %For recording
     frame = getframe;
     writeVideo(aviobj,frame);       
     drawnow             
                  
    %If H is not a number
    if all(all(isnan(H))), break, end
    
 end
 close(aviobj);
