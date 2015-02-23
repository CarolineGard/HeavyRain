clc

g = 9.81;
dx = 1.0;
dy = 1.0;
dt = 0.05; 

H = ones(41:41); 
U = zeros(41:41); %hastighet i x-led
V = zeros(41:41); %hastighet i y-led (obs negativ!)

Hx = zeros(40:40);
Ux = zeros(40:40);
Vx = zeros(40:40);

Hy = zeros(40:40);
Uy = zeros(40:40);
Vy = zeros(40:40);

C = zeros(40:40);

%Skapa initial displacement
[x, y] = meshgrid( linspace(-4, 4, 10) );
R = sqrt(x.^2 + y.^2) + eps;
Z = (sin(R)./R);
Z = max(Z,0);

%Skapa displacement to the height matrix
w = size(Z,1);
t = 10:w+9;
l = 20:w+19;
H(t,l) = H(t, l) + Z;

%Sätt upp grids
clf 
shg
n = 40;
grid = surf(H);
axis([0 41 0 41 -1 3]);
hold all;

 while 1 == 1  
    
    C = abs(U(i,j)) + abs(V(i,j));
    set(grid, 'zdata', H(i,j), 'cdata', C);
    axis off;
    drawnow   
    
     % Reflective boundary conditions Inte vår kod!!
     H(:,1) = H(:,2);      U(:,1) = U(:,2);       V(:,1) = -V(:,2);
     H(:,n+2) = H(:,n+1);  U(:,n+2) = U(:,n+1);   V(:,n+2) = -V(:,n+1);
     H(1,:) = H(2,:);      U(1,:) = -U(2,:);      V(1,:) = V(2,:);
     H(n+2,:) = H(n+1,:);  U(n+2,:) = -U(n+1,:);  V(n+2,:) = V(n+1,:);
    
    %Lax- Wendroff 2-stegsmetod
    
    %FIRST STEP
    %Dela allt med H i samma punkt

    %X-LED
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

    %Y-LED
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


                  
    %ANDRA STEGET
    i = 2:38+1;
    j = 2:38;

    H(i, j) = H(i, j) - (dt/dx)*(Ux(i,j-1) - Ux(i-1, j-1)) - (dt/dy)*(Vy(i-1,j) - Vy(i-1,j-1));


    U(i, j) = U(i, j) - (dt/dx)*(Ux(i,j-1).^2 ... 
                      + (0.5*g*Hx(i,j-1) - (Ux(i-1,j-1).^2 + 0.5*g*Hx(i-1,j-1) ))) ...
                      - (dt/dy)*(Vy(i-1,j).*Uy(i-1,j) - Vy(i-1,j-1).*Uy(i-1,j-1)); 
              
    V(i, j) = V(i, j) - (dt/dy)*((Vy(i-1,j).^2 ... 
                      + 0.5*g*Hy(i-1,j)) - (Vy(i-1,j-1).^2 + 0.5*g*Hy(i-1,j-1))) ...
                      - (dt/dx)*(Ux(i,j-1).*Vx(i,j-1) - Ux(i-1, j-1).*Vx(i-1,j-1)); 


    if all(all(isnan(H))), break, end
    
 end
