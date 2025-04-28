close all;
clear all;
clc
%% input parameters
format long

lambda = 1;

% grid parameters
px = 201; % number of points in source line
py = 201; % number of points in source line

pzr = 101; % number of points in receiver line
pphir = 101;

dx = 0.5*lambda; % spacing between source points 
dy = 0.5*lambda; % spacing between source points 

dzr = 0.5*lambda; % spacing between receiver points

%%
lambda = 1;
k = 2*pi/lambda;

Xds = (px-1)*dx; % dimensions of the spaces
Yds = (py-1)*dy;
Zdr = (pzr-1)*dzr;

L = 1*Xds; % longitudinal separation distance between the spaces
Rc = 0.5*Yds/2; % cylinder radius
%%

% maximum allowed value for the source spacing distances 
sin_th = (0.5*Yds+Rc)/sqrt((0.5*Xds)^2+(0.5*Yds+Rc)^2+L^2);
dymax = (1/sqrt(2))*lambda/(sin_th);
fprintf('Max source spacing dy (x lambda): %f \n',(1/sqrt(2))*lambda/(sin_th))
%%
% setting the source and receiver grids
xs = linspace(-(px-1)*dx/2,(px-1)*dx/2,px);
ys = linspace(-(py-1)*dy/2,(py-1)*dy/2,py);


zr = linspace(L,L+(pzr-1)*dzr,pzr);
delta_phi = 2*pi/pphir; 
phir = 0:delta_phi:2*pi-delta_phi;
rr = Rc*ones(1,length(phir));


[Xs,Ys] = meshgrid(xs,ys);
Zs = zeros(size(Xs));
[Phir,Zr] = meshgrid(phir,zr);
Rr = Rc*ones(size(Phir));
Xr = Rr.*cos(Phir);
Yr = Rr.*sin(Phir);

Ns = px*py;
Nr = pzr*pphir;

X1s = reshape(Xs,1,[]);
Y1s = reshape(Ys,1,[]);
Z1s = reshape(Zs,1,[]);

X1r = transpose(reshape(Xr,1,[]));
Y1r = transpose(reshape(Yr,1,[]));
Z1r = transpose(reshape(Zr,1,[]));

%% % eigenproblem computation
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
Mpw = 1000;
[psi,D,V] = svds(ghg,Mpw);
% psi - source eigenfunctions
s2 = diag(D); % squared amplitude of the eigenvalues
S =  sum(abs(s2)); % sum rule S 
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
% plotting of coupling strengths and a_j coefficients 
figure
plot(s2*100,'r-.','LineWidth',2);grid
xlabel('Communication mode index j')
ylabel('Absolute coupling strength (x 100)')
xlim([0 Mpw])
%%

save('data/cylindric_s_201_201_r_101_101_0_5lb_Rc_0_25Ys.mat','psi','phi','s2','S','Mpw');

%%

