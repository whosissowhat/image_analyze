%function euler2D()

%-----------------------------
% initial data:
%----------------------------
%  - rho inside  circle = 1
%  - u inside   circle = 0
%  - v inside   circle = 0
%  - p inside   circle = 1
%  - rho outside circle = 0.125
%  - u outside circle = 0
%  - v outside circle = 0
%  - p outside circle = 0.1
% radius of circle r=0.5
% final time T=0.2
% CLF constant =0.85;
% Poisson adibatic constant gamma=1.4
%----------------------------
clear
close all;
clc;

xMin= -1; xMax= 1;
yMin= -1; yMax= 1;
%---------------------------

Nx=50; % number of control volumes in x-direction, i.e. number of columns in matrix w !!!
Ny=50; % number of control volumes in y-direction, i.e. number of rows in matrix w !!!

hx=(xMax-xMin)/Nx; % length of control volume (CV) in x-direction
hy=(yMax-yMin)/Ny; % length of control volume (CV) in y-direction

%x=...  % vector of x-coordinadets of midpoints CVs (fictitious too), i.e. x=(x0,..,x(Nx+2))
x = zeros(1,Nx+2);
x(1) = xMin-hx/2;
for i=1:Nx+1
    x(i+1) = x(i) + hx;
end

%y=...  % vector of x-coordinadets of midpoints CVs (fictitious too), i.e. y=(y0,..,y(Ny+2))
y = zeros(1,Ny+2);
y(1) = yMin-hy/2;
for i=1:Ny+1
    y(i+1) = y(i) + hy;
end
% counter time steps

dI = 1;
uI = 0;
vI = 0;
pI = 1;

dO = 0.125;
uO = 0;
vO = 0;
pO = 0.1;
r = 0.5;


T=0.2;                           % final time
gamma=1.4;                       % Poisson adibatic constant
g = gamma-1;
CFL=0.85;                        % CLF constant
tn=0;                            % initial time
t=0;
% writing of initial conditions
wInside =[dI; dI*uI; dI*vI; pI/(gamma-1) + 1/2*dI*(uI^2 + vI^2)];    % column vector of conservative variables, inside circle
wOutside=[dO; dO*uO; dO*vO; pO/(gamma-1) + 1/2*dO*(uO^2 + vO^2)]; % column vector of conservative variables, outside circle

w=zeros(Ny+2,Nx+2,4);  % declaration w

for i=1:Ny+2 % writing of initial state of conservative variables into matrix w
    for j=1:Nx+2
        if x(i)^2 + y(j)^2 <= r^2  % if I'm inside circle with radius r
            w(i,j,:)=wInside; % vector of conservative variables inside circle
        else
            w(i,j,:)=wOutside; % vector of conservative variables outside circle
        end
    end
end

% main loop
while abs(T-tn)>10^(-10)
    [Fx,Fy,LamMax] = Vijayasundaram2D(w,Nx,Ny,gamma);

    tau = min(hx,hy)*CFL/(sqrt(2)*LamMax); % new length of time step

    if (tn+tau) > T
        tau = T - tn;
    end
    tn = tn + tau;

    for i = 2:Ny+1
        for j = 2:Nx+1


            w(i,j,:) = w(i,j,:) - tau*((Fx(i,j,:)-Fx(i,j-1,:))/hx + (Fy(i-1,j,:)-Fy(i,j,:))/hy);

        end
    end

    %% extrapolate Boundary Conditions
    %w(1,2:end-1,:)= w(2,2:end-1,:)
    %w(end,2:end-1,:)= w(end-1,2:end-1,:)
    %w(1:end,1,:)= w(1:end,2,:)
    %w(1:end,end,:)= w(1:end,end-1,:)
    
    %North Side
    w(1,:,:)=w(2,:,:);
    % South Side
    w(end,:,:)=w(end-1,:,:);
    % West Side
    w(:,1,:)=w(:,2,:);
    % East Side
    w(:,end,:)=w(:,end-1,:);

    t=t+1; % counter time steps
end % end of main loop while

%---------------!!! do not change !!!! -------->beginning----------------------
% comparing with 1D Riemann solver:
fid=fopen('e1godf00.ini','rt');
tline=fgetl(fid); L=sscanf(tline,'%e',1);
tline=fgetl(fid); DIAPH=sscanf(tline,'%e',1);
tline=fgetl(fid); CELLS=sscanf(tline,'%d',1);
tline=fgetl(fid); GAMMA=sscanf(tline,'%e',1);
tline=fgetl(fid); T=sscanf(tline,'%e',1);
tline=fgetl(fid); DL=sscanf(tline,'%e',1);
tline=fgetl(fid); UL=sscanf(tline,'%e',1);
tline=fgetl(fid); PL=sscanf(tline,'%e',1);
tline=fgetl(fid); DR=sscanf(tline,'%e',1);
tline=fgetl(fid); UR=sscanf(tline,'%e',1);
tline=fgetl(fid); PR=sscanf(tline,'%e',1);
tline=fgetl(fid); CFL=sscanf(tline,'%e',1);
tline=fgetl(fid); PSCALE=sscanf(tline,'%e',1);
fclose(fid);
fprintf('pocita se presne reseni Eulerovych rovnic\n');
Xin={L,DIAPH,1000,GAMMA,T,DL,UL,PL,DR,UR,PR,PSCALE};
[XE,~]=e1r(Xin,1);
%---------------!!! do not change !!!! -------->end------------------------------

%=========================
% picture in 1D
for k=1:4
    for j=1:Ny+2
        XaxisData(j,k)=w((Nx+2)/2,j,k);
    end
end

% calculate primitive variables from XaxisData (for 1D plot)
rho1D = XaxisData(:,1); % density
u1D = XaxisData(:,2)./rho1D;  % velocity u
p1D = g*(XaxisData(:,4)-0.5*rho1D.*u1D.*u1D); % pressure

figure(1)
subplot(2,2,1);
plot(XE(1,:),XE(2,:),'b-',x,rho1D,'r')
title(sprintf('density, time =%1.4g',T));
subplot(2,2,2);
plot(XE(1,:),XE(3,:),'b-',x,u1D,'r')
title(sprintf('velocity u, time =%1.4g',T));
subplot(2,2,3);
plot(XE(1,:),XE(4,:),'b-',x ,p1D,'r')
title(sprintf('pressure p, time =%1.4g',T));

%=========================
% picture in 2D
%=========================

% calculate primitive variables from w (for 2D plot)
rho2D = w(:,:,1);        % density
u2D =   w(:,:,2)./rho2D; % velocity
v2D =   w(:,:,3)./rho2D; % velocity
p2D =   g*(w(:,:,4)-0.5*rho2D.*(u2D.*u2D+v2D.*v2D)); % pressure


figure(2)
title(sprintf('2D Euler Equations, time = %6.2f'));
[xx,yy]=meshgrid(x,y);
subplot(2,2,1)
surf(xx,yy,rho2D)
title(sprintf('density, time =%1.4g',T));
subplot(2,2,2)
surf(xx,yy,u2D)
title(sprintf('velocity u, time =%1.4g',T));
subplot(2,2,3)
surf(xx,yy,v2D)
title(sprintf('velocity v,, time =%1.4g',T));
subplot(2,2,4)
surf(xx,yy,p2D)
title(sprintf('pressure p, time =%1.4g',T));


%---------------------------------------------------------------------------
%   calculation of  fluxes
%---------------------------------------------------------------------------

function [fluxX,fluxY,LamMax]=Vijayasundaram2D(w,Nx,Ny,gamma)
% fluxX = zeros(Nx,Ny+1,4);
% fluxY = zeros(Nx+1,Ny,4);

fluxX = zeros(Ny+2,Nx+1,4);
fluxY = zeros(Ny+1,Nx+2,4);
LamMax=0;
g = gamma-1;

%% fluxes in x-direction
n1 = 1; n2 = 0;
for i = 1:Ny+2
    for j = 1:Nx+1

        wL(:,1) = w(i,j,:);% values of conserv. variables, left side of interface x_{i-1/2}

        dL = wL(1);
        uL= wL(2)/dL;
        vL = wL(3)/dL;
        EL = wL(4);
        pL = (gamma-1)*(EL - 0.5*dL*(uL^2 + vL^2));
        aL = sqrt(gamma*pL/dL);

        if ~isreal(aL)
            error('Program terminates, complex');
        end
        HL = (EL + pL)/dL;
        vnL = uL*n1 + vL*n2;
        v_magL = uL^2 + vL^2;

        TL = [1       1           0          1;
            uL-aL*n1  uL           n2         uL+aL*n1;
            vL-aL*n2  vL           -n1        vL+aL*n2;
            HL-aL*vnL  0.5*v_magL  n2*uL-n1*vL  HL+aL*vnL];

        TinvL = 1/(2*aL^2)*...
            [0.5*g*v_magL+aL*vnL  -aL*n1-g*uL  -aL*n2-g*vL  g;
            2*aL^2-g*v_magL     2*g*uL      2*g*vL      -2*g;
            2*aL^2*(vL*n1-uL*n2)  2*aL^2*n2   -2*aL^2*n1  0;
            0.5*g*v_magL-aL*vnL  aL*n1-g*uL   aL*n2-g*vL   g];

        DpL = diag([max(vnL-aL,0), max(vnL,0), max(vnL,0), max(vnL+aL,0)]);
        PpL = TL * DpL * TinvL;


        % For the Right Side of Interface
        wR(:,1) = w(i,j+1,:); % values of conserv. variables, right side of interface x_{i-1/2}
        dR = wR(1);
        uR= wR(2)/dR;
        vR = wR(3)/dR;
        ER = wR(4);
        pR = (gamma-1)*(ER - 0.5*dR*(uR^2 + vR^2));
        aR = sqrt(gamma*pR/dR);
        if ~isreal(aR)
            error('Program terminates, complex value for speed of sound');
        end
        HR = (ER + pR)/dR;
        vnR = uR*n1 + vR*n2;
        v_magR = uR^2 + vR^2;

        TR = [1       1           0          1;
            uR-aR*n1  uR           n2         uR+aR*n1;
            vR-aR*n2  vR           -n1        vR+aR*n2;
            HR-aR*vnR  0.5*v_magR  n2*uR-n1*vR  HR+aR*vnR];

        TinvR = 1/(2*aR^2)*...
            [0.5*g*v_magR+aR*vnR  -aR*n1-g*uR  -aR*n2-g*vR  g;
            2*aR^2-g*v_magR     2*g*uR      2*g*vR      -2*g;
            2*aR^2*(vR*n1-uR*n2)  2*aR^2*n2   -2*aR^2*n1  0;
            0.5*g*v_magR-aR*vnR  aR*n1-g*uR   aR*n2-g*vR   g];

        DmR = diag([min(vnR-aR,0), min(vnR,0), min(vnR,0), min(vnR+aR,0)]);
        PmR = TR * DmR * TinvR;

        % Calculate middle state variables
w_hat = (wL + wR)/2;
d_hat = w_hat(1);  
u_hat = w_hat(2)/d_hat;
v_hat = w_hat(3)/d_hat;
E_hat = w_hat(4);
p_hat = (gamma-1)*(E_hat - 0.5*d_hat*(u_hat^2 + v_hat^2));
a_hat = sqrt(gamma*p_hat/d_hat);

% Compute fluxes using the middle state
vn_hat = u_hat*n1 + v_hat*n2;
v_mag_hat = u_hat^2 + v_hat^2;

TL_hat = [1           1                0             1;
          u_hat-a_hat*n1  u_hat           n2            u_hat+a_hat*n1;
          v_hat-a_hat*n2  v_hat           -n1           v_hat+a_hat*n2;
          E_hat-a_hat*vn_hat  0.5*v_mag_hat  n2*u_hat-n1*v_hat  E_hat+a_hat*vn_hat];

Tinv_hat = 1/(2*a_hat^2)*...
    [0.5*g*v_mag_hat+a_hat*vn_hat  -a_hat*n1-g*u_hat  -a_hat*n2-g*v_hat  g;
     2*a_hat^2-g*v_mag_hat     2*g*u_hat      2*g*v_hat      -2*g;
     2*a_hat^2*(v_hat*n1-u_hat*n2)  2*a_hat^2*n2   -2*a_hat^2*n1  0;
     0.5*g*v_mag_hat-a_hat*vn_hat  a_hat*n1-g*u_hat   a_hat*n2-g*v_hat   g];

Dp_hat = diag([max(vn_hat-a_hat,0), max(vn_hat,0), max(vn_hat,0), max(vn_hat+a_hat,0)]);
Pp_hat = TL_hat * Dp_hat * Tinv_hat;

Dm_hat = diag([min(vn_hat-a_hat,0), min(vn_hat,0), min(vn_hat,0), min(vn_hat+a_hat,0)]);
Pm_hat = TL_hat * Dm_hat * Tinv_hat;

% Update fluxX using middle state
fluxX(i,j,:) = Pp_hat*wL + Pm_hat*wR;

        LamMax = max([LamMax,abs(vnL)+aL]);
    end
end

%% fluxes in y-direction
n1=0;n2=1;
for i = 1:Ny+1
    for j = 1:Nx+2
        wL(:,1) = w(i+1,j,:);% values of conserv. variables, left side of interface x_{i-1/2}

        dL = wL(1);
        uL= wL(2)/dL;
        vL = wL(3)/dL;
        EL = wL(4);
        pL = (gamma-1)*(EL - 0.5*dL*(uL^2 + vL^2));
        aL = sqrt(gamma*pL/dL);

        if ~isreal(aL)
            error('Program terminates, complex');
        end
        HL = (EL + pL)/dL;
        vnL = uL*n1 + vL*n2;
        v_magL = uL^2 + vL^2;

        TL = [1       1           0          1;
            uL-aL*n1  uL           n2         uL+aL*n1;
            vL-aL*n2  vL           -n1        vL+aL*n2;
            HL-aL*vnL  0.5*v_magL  n2*uL-n1*vL  HL+aL*vnL];

        TinvL = 1/(2*aL^2)*...
            [0.5*g*v_magL+aL*vnL  -aL*n1-g*uL  -aL*n2-g*vL  g;
            2*aL^2-g*v_magL     2*g*uL      2*g*vL      -2*g;
            2*aL^2*(vL*n1-uL*n2)  2*aL^2*n2   -2*aL^2*n1  0;
            0.5*g*v_magL-aL*vnL  aL*n1-g*uL   aL*n2-g*vL   g];

        DpL = diag([max(vnL-aL,0), max(vnL,0), max(vnL,0), max(vnL+aL,0)]);
        PpL = TL * DpL * TinvL;


        % For the West Side of Interface
        wR(:,1) = w(i,j,:); % values of conserv. variables, right side of interface x_{i-1/2}
        dR = wR(1);
        uR= wR(2)/dR;
        vR = wR(3)/dR;
        ER = wR(4);
        pR = (gamma-1)*(ER - 0.5*dR*(uR^2 + vR^2));
        aR = sqrt(gamma*pR/dR);
        if ~isreal(aR)
            error('Program terminates, complex value for speed of sound');
        end
        HR = (ER + pR)/dR;
        vnR = uR*n1 + vR*n2;
        v_magR = uR^2 + vR^2;

        TR = [1       1           0          1;
            uR-aR*n1  uR           n2         uR+aR*n1;
            vR-aR*n2  vR           -n1        vR+aR*n2;
            HR-aR*vnR  0.5*v_magR  n2*uR-n1*vR  HR+aR*vnR];

        TinvR = 1/(2*aR^2)*...
            [0.5*g*v_magR+aR*vnR  -aR*n1-g*uR  -aR*n2-g*vR  g;
            2*aR^2-g*v_magR     2*g*uR      2*g*vR      -2*g;
            2*aR^2*(vR*n1-uR*n2)  2*aR^2*n2   -2*aR^2*n1  0;
            0.5*g*v_magR-aR*vnR  aR*n1-g*uR   aR*n2-g*vR   g];

        DmR = diag([min(vnR-aR,0), min(vnR,0), min(vnR,0), min(vnR+aR,0)]);
        PmR = TR * DmR * TinvR;

        % Calculate middle state variables
w_hat = (wL + wR)/2;
d_hat = w_hat(1);  
u_hat = w_hat(2)/d_hat;
v_hat = w_hat(3)/d_hat;
E_hat = w_hat(4);
p_hat = (gamma-1)*(E_hat - 0.5*d_hat*(u_hat^2 + v_hat^2));
a_hat = sqrt(gamma*p_hat/d_hat);

% Compute fluxes using the middle state
vn_hat = u_hat*n1 + v_hat*n2;
v_mag_hat = u_hat^2 + v_hat^2;

TL_hat = [1           1                0             1;
          u_hat-a_hat*n1  u_hat           n2            u_hat+a_hat*n1;
          v_hat-a_hat*n2  v_hat           -n1           v_hat+a_hat*n2;
          E_hat-a_hat*vn_hat  0.5*v_mag_hat  n2*u_hat-n1*v_hat  E_hat+a_hat*vn_hat];

Tinv_hat = 1/(2*a_hat^2)*...
    [0.5*g*v_mag_hat+a_hat*vn_hat  -a_hat*n1-g*u_hat  -a_hat*n2-g*v_hat  g;
     2*a_hat^2-g*v_mag_hat     2*g*u_hat      2*g*v_hat      -2*g;
     2*a_hat^2*(v_hat*n1-u_hat*n2)  2*a_hat^2*n2   -2*a_hat^2*n1  0;
     0.5*g*v_mag_hat-a_hat*vn_hat  a_hat*n1-g*u_hat   a_hat*n2-g*v_hat   g];

Dp_hat = diag([max(vn_hat-a_hat,0), max(vn_hat,0), max(vn_hat,0), max(vn_hat+a_hat,0)]);
Pp_hat = TL_hat * Dp_hat * Tinv_hat;

Dm_hat = diag([min(vn_hat-a_hat,0), min(vn_hat,0), min(vn_hat,0), min(vn_hat+a_hat,0)]);
Pm_hat = TL_hat * Dm_hat * Tinv_hat;

% Update fluxY using middle state
fluxY(i,j,:) = Pp_hat*wL + Pm_hat*wR;

        LamMax = max([LamMax,abs(vnL)+aL]);
    end
end
end


