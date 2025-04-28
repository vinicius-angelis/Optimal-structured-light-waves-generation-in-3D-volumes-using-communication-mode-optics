close all;
clear all;
clc
%% input parameters
format long

lambda = 1;

px = 101; % number of points in source line
py = 301; % number of points in source line

pzr = 51; % number of points in receiver line
pxr = 51; % number of points in receiver line
pyr = 10; % number of points in receiver line

dx = 0.5*lambda; % spacing between source points 
dy = 0.5*lambda; % spacing between source points 

dzr = 1*lambda; % spacing between receiver points
dxr = 1*lambda; % spacing between receiver points
dyr = 15*lambda; % spacing between receiver points

%%

k = 2*pi/lambda;

Xds = (px-1)*dx; % dimensions of the spaces
Yds = (py-1)*dy;
Xdr = (pxr-1)*dxr;
Ydr = (pyr-1)*dyr;
Zdr = (pzr-1)*dzr;

L = 1*Xds; % longitudinal separation distance between the spaces

%%

% maximum allowed value for the source spacing distances 
aux1 = 0.5*(Yds+Ydr);
aux2 = sqrt(0.25*(Xds+Xdr)^2+0.25*(Yds+Ydr)^2+(L)^2);
sin_th = aux1/aux2;
dymax = (1/sqrt(2))*lambda/(sin_th);
fprintf('Max source spacing dy (x lambda): %f \n', dymax)

%%
% setting the source and receiver grids
xs = linspace(-(px-1)*dx/2,(px-1)*dx/2,px);
ys = linspace(-(py-1)*dy/2,(py-1)*dy/2,py);


zr = linspace(L,L+(pzr-1)*dzr,pzr);
xr = linspace(-(pxr-1)*dxr/2,(pxr-1)*dxr/2,pxr);
yr = linspace(-(pyr-1)*dyr/2,(pyr-1)*dyr/2,pyr);

[Xs,Ys] = meshgrid(xs,ys);

[Xr,Yr,Zr] = meshgrid(xr,yr,zr);

Zs = zeros(size(Xs));

Ns = length(xs)*length(ys);
Nr = pzr*pxr*pyr;


X1s = reshape(Xs,1,[]);
Y1s = reshape(Ys,1,[]);
Z1s = reshape(Zs,1,[]);

X1r = transpose(reshape(Xr,1,[]));
Y1r = transpose(reshape(Yr,1,[]));
Z1r = transpose(reshape(Zr,1,[]));


%% %% % eigenproblem computation
tic
dist = zeros(Nr,Ns);
for i=1:Nr
    for j=1:Ns
        dist(i,j) = sqrt((X1s(j)-X1r(i)).^2 + (Y1s(j)-Y1r(i)).^2 + (Z1s(j)-Z1r(i)).^2);
    end
end

g = -(1/(4*pi))*exp(1i*k*dist)./dist;
ghg = ctranspose(g)*g;
toc

tic
Mpw = 3500; % number of communication modes
[psi,D,V] = svds(ghg,Mpw);
% psi - source eigenfunctions
s2 = diag(D); % squared amplitude eigenvalues
S =  sum(abs(s2)); % sum S 
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
Mpw = 3500;
% plotting of coupling strengths and a_j coefficients 
figure
plot(100*abs(s2),'o','LineWidth',2);grid
%hold on
%xline(64,'-.','LineWidth',1.5)
%xlim([0 size(phi,1)])
xlabel('Mode index')
ylabel('Absolute coupling strength (x 100)')
%ylabel('| s_j |^2 as % of sum rule S')
title(['Eigenvalues      S = ' num2str(S),'']);
xlim([0 Mpw])
%%
save('data/3D_svd_L_Xs_dx_dy_0_5lb_dxr_dzr_1lb_dyr_15lb_s_101_301_pt_r_51_51_10_pt.mat','psi','phi','s2','S','Mpw');
%%
