close all;
clear all;
clc;
%% input parameters
format long
lambda = 1;

% grid parameters
px = 111; % number of points of the source plane in x direction
py = 222; % number of points of the source plane in y direction

pzr = 101; % number of points in receiver horizontal plane in z direction
pxr = 101; % number of points in receiver horizontal plane in x direction

dx = 1*lambda; % spacing distance between source points (x direction) 
dy = 1*lambda; % spacing distance between source points (y direction)


dzr = 1*lambda; % spacing distance between receiver points (z direction)
dxr = 1*lambda; % spacing distance between receiver points (x direction)
%%
L = 1*(px-1)*dx; % longitudinal separation distance between the spaces
lambda = 1;
k = 2*pi/lambda; % wavenumber

Xds = (px-1)*dx; % dimensions of the spaces
Yds = (py-1)*dy;
Xdr = (pxr-1)*dxr;
Zdr = (pzr-1)*dzr;
Ydr = 0;

% maximum allowed value for the source spacing distances 
aux1 = 0.5*(Yds+Ydr);
aux2 = sqrt(0.25*(Xds+Xdr)^2+0.25*(Yds+Ydr)^2+(L)^2);
sin_th = aux1/aux2;
fprintf('Max source spacing dy: %f \n',(1/sqrt(2))*lambda/(sin_th))

% number of effective longitudinal Nz and transverse Nx modes 
Nz = (sqrt(L^2+(Yds/2).^2)-sqrt((L+Zdr)^2+(Yds/2).^2)+Zdr);
Nx = Xds*Xdr/(L+Zdr);
Nxp = Xds*Xdr/(L+Zdr*0.0);
fprintf('Nz: %f \n',Nz)
fprintf('Nx: %f \n',Nx)

% vertical positions of each longitudinal effective mode
nz = 1:ceil(Nz)+1;
Ynz = sqrt(nz).*sqrt(2*Zdr-nz).*sqrt(4*L*(L+Zdr)+2*Zdr.*nz-nz.^2)./(2*(Zdr-nz));
figure
stem(nz,Ynz,'k','LineWidth',1.5);grid
hold on
plot(nz,ones(length(nz))*(Yds/2),'r-.','LineWidth',1.5)

%%
% setting the source and receiver grids
xs = linspace(-(px-1)*dx/2,(px-1)*dx/2,px);
ys = linspace(-(py-1)*dy/2,(py-1)*dy/2,py);

zr = linspace(L,L+(pzr-1)*dzr,pzr);
xr = linspace(-(pxr-1)*dxr/2,(pxr-1)*dxr/2,pxr);

[Xs,Ys] = meshgrid(xs,ys);
[Zr,Xr] = meshgrid(zr,xr);
Zs = zeros(size(Xs));
Yr = zeros(size(Zr));

Ns = px*py;
Nr = pzr*pxr;

X1s = reshape(Xs,1,[]);
Y1s = reshape(Ys,1,[]);
Z1s = reshape(Zs,1,[]);

X1r = reshape(Xr,1,[]);
Y1r = reshape(Yr,1,[]);
Z1r = reshape(Zr,1,[]);

%% % eigenproblem computation
tic
dist = zeros(Nr,Ns);
for i=1:Nr
    for j=1:Ns
        dist(i,j) = sqrt((X1s(j)-X1r(i)).^2 + (Y1s(j)-Y1r(i)).^2 + (Z1s(j)-Z1r(i)).^2);
    end
end

g = -(1/(4*pi))*exp(1i*k*dist)./dist;
norm = -4*pi*L;
gn = g*norm;

gh = ctranspose(g);
ghg = gh*g;
toc

tic
Mpw = 1400; % number of communication modes
[psi,D,V] = svds(ghg,Mpw);
% psi - source eigenfunctions
s2 = diag(D); % squared amplitude of the eigenvalues
S =  sum(s2); % sum S 
toc

tic
gaux = g*psi;
toc
%%
% phi - receiver eigenfunctions
tic
phi = zeros(Nr,Mpw);
for i=1:Mpw
  phi(:,i) = (1/sqrt(s2(i)))*gaux(:,i);
end
toc
%%
% plotting of coupling strengths 
figure
plot(s2*100,'r-.','LineWidth',2);grid
xlabel('Communication mode index j')
ylabel('Absolute coupling strength (x 100)')
title('d_x = d_y = d_{x,r} = d_{z,r} = \lambda   p_x = 111   p_{x,r} = p_{z,r} = 101   L = 110\lambda');
xlim([0 Mpw])
%%
save('data/2D_svd_dx_dy_1lb_dxr_dzr_1lb_s_111_222_pt_r_101_101_pt.mat','psi','phi','s2','S','Mpw');
%%