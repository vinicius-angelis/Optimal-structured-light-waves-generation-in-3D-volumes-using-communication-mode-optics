close all;
clear all;
clc
%%
% this script compiles the results shown in Fig. S12 in the supplementary
% material for a transverse source plane and a receiving cylinder surface


%the code also computes the hologram phase profile and save it 
% in a png file (e.g., holog_test.png) to be encoded onto the SLM.

%Run the code: 'main_cylindric_eigenproblem.m' to solve the SVD 
% for a particular source and receiving configuration.

%Vinicius Soares de Angelis
%viniciusangelis@g.harvard.edu


%https://arxiv.org/abs/2411.10865


tic
% loading eigenfunctions and eigenvalues of the source-receiving
% configuration listed in Table 1
load('data/cylindric_s_201_201_r_101_101_0_5lb_Rc_0_25Ys.mat')
toc
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
k = 2*pi/lambda;

Xds = (px-1)*dx; % dimensions of the spaces
Yds = (py-1)*dy;
Zdr = (pzr-1)*dzr;

L = 1*Xds; % longitudinal separation distance between the spaces
Rc = 0.5*Yds/2; % cylinder radius

%%
% maximum allowed value for the source spacing distances 
aux1 = 0.5*(Yds+2*Rc);
aux2 = sqrt(0.25*(Xds+2*Rc)^2+0.25*(Yds+2*Rc)^2+(L)^2);
sin_th = aux1/aux2;
dymax = (1/sqrt(2))*lambda/(sin_th);
fprintf('Max source spacing dy (x lambda): %f \n',dymax)
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

% plotting of coupling strengths
figure
plot(s2*100,'r-.','LineWidth',2);grid
xlabel('Communication mode index j')
ylabel('Absolute coupling strength (x 100)')
%xlim([0 Mpw])
%%

Ao1 = imread('templates/template_numbers.png');
Ao1 = flip(Ao1,1);

Ao = imresize(Ao1, [pphir pzr]); % resize input image
Ao = rgb2gray(Ao);
I1 = single(Ao); % converts gray scale values to float
I1 = I1/max(max(I1)); % normalization of input image intensity 
Faux = sqrt(I1); % normalized field amplitude of input image  

% modulating target intensity profile by a complex exponential exp(iQz)
Q = 0.95*k;
F = Faux.*exp(1i*(Q).*Zr);

R = linspace(0,0,256);
G = linspace(0,1,256);
B = linspace(0,0,256);
map = [R',G',B'];

figure
surf(Zr,Xr,Yr,abs(F).^2)
shading interp
colormap(map)
colorbar
set(gca,'FontSize',10,'FontWeight','bold')
title('Desired receiver values');
%view([0 90])
xlabel('z_r / \lambda')
ylabel('x_r / \lambda')
xlim([zr(1) zr(end)])


F1 = reshape(F,1,[]);
% coefficients of desired profile in receiver space
a = zeros(size(phi,2),1);
for i=1:length(a)
a(i) = dot((phi(:,i)),(F1));
end
 
figure
plot(abs(a),'o','LineWidth',1.5);grid
hold on
plot(max(abs(a))*s2/max(s2),'r-.','LineWidth',1.5)
xlabel('Communication mode index j')
ylabel('| a_j |')
xlim([0 Mpw])     
lg = legend('| a_j |','|s_j|^2 (reference)');
set(lg,'FontSize',15) 

%%
% required source function computation

tic
Mp0 = 1;  % initial and final communication index values assigned to summation
Mp = 800;  

SM_psi = zeros(size(psi,1));
for i=Mp0:Mp
SM_psi(:,i) = (1/sqrt(s2(i)))*(a(i))*(psi(:,i));
end
S_psi = (sum(SM_psi,2));


Sm = reshape(S_psi,py,px);
Imax = max(max(abs(Sm).^2));
figure
surf(Xs,Ys,abs(Sm).^2)
shading interp
colormap(map)
colorbar 
xlabel('x / \lambda')
ylabel('y / \lambda')
set(gca,'FontSize',10,'FontWeight','bold')
view([0 90])
xlim([xs(1) xs(end)])
ylim([ys(1) ys(end)])
title(['Source Squared amplitude    z = ' num2str(0),' \lambda']);
toc
%%
% resulting wave in the receiving cylinder surface
tic
phiphi = interpn(phir,2);
phiphi(end) = 2*pi;
Qz = 201;
zz = linspace(L,L+Zdr,Qz);
[PHI,ZZ] = meshgrid(phiphi,zz);
RHO = Rc*ones(size(PHI));
XX = RHO.*cos(PHI);
YY = RHO.*sin(PHI);

aux = zeros(length(zz),length(phiphi));

        for q=1:Ns
             da = sqrt((XX-X1s(q)).^2+(YY-Y1s(q)).^2+(ZZ-Z1s(q)).^2);
             aux = aux + (exp(1i*k*da)./da)*S_psi(q);
        end
        Uz = -(1/(4*pi))*aux; 

figure
surf(ZZ-L,XX,YY,abs(Uz).^2)
shading interp
colormap(map)
colorbar 
xlabel('z / \lambda')
ylabel('x / \lambda')
zlabel('\Phi')
set(gca,'FontSize',10,'FontWeight','bold')
title('Resulting wave');
xlim([zz(1)-L zz(end)-L])


figure
surf(ZZ-L,PHI*180/pi,abs(Uz).^2);grid
shading flat
colormap(map)
colorbar 
h = gca;
h.Visible = 'On';
axis square
shading flat
xlabel('z (mm)')
ylabel('\phi (deg)')
set(gca,'FontSize',10,'FontWeight','bold')
view([0 90])
xlim([zz(1)-L zz(end)-L])
ylim([0 360])
%%
% resulting wave in the mid longitudinal distance
% z = L + Zr/2 of receiving cylinder surface
N = 201;
Zp = L+0.5*Zdr;
xx = linspace(xs(1),xs(end),N);
yy = linspace(ys(1),ys(end),N);

[XX,YY] = meshgrid(xx,yy);
ZZ = ones(size(XX))*Zp;

aux = zeros(N,N);
        for q=1:Ns
             da = sqrt((XX-X1s(q)).^2+(YY-Y1s(q)).^2+(ZZ-Z1s(q)).^2);
             aux = aux + (exp(1i*k*da)./da)*(S_psi(q));
        end
 U = -(1/(4*pi))*aux; 
               
figure
surf(XX,YY,abs(U).^2)
shading interp
colormap(map)
colorbar 
xlim([xx(1) xx(end)])
ylim([yy(1) yy(end)])
xlabel('x / \lambda')
ylabel('y / \lambda')
set(gca,'FontSize',10,'FontWeight','bold')
title(['z = ' num2str(Zp),' \lambda']);
%title(['Mode ' num2str(Mp),'        z = ' num2str(Zp),' \lambda']);
view([0 90])
%%
% computing the complex field at z = L 

%SLM screen size: 600 pts in ys direction

tic
Ny = 600;
Nx = ceil(Ny/(Yds/Xds));
if mod(Nx,2)==1
    Nx = Nx+1;
end

xx = linspace(xs(1),xs(end),Nx);
yy = linspace(ys(1),ys(end),Ny);

[XX,YY] = meshgrid(xx,yy);
ZZ = ones(size(XX))*L;

aux = zeros(Ny,Nx);
        for q=1:Ns
             da = sqrt((XX-X1s(q)).^2+(YY-Y1s(q)).^2+(ZZ-Z1s(q)).^2);
             aux = aux + (exp(1i*k*da)./da)*(S_psi(q));
        end
 U = -(1/(4*pi))*aux;    

figure 
surf(XX,YY,abs(U).^2)
shading interp
colormap(map)
colorbar 
xlim([xx(1) xx(end)])
ylim([yy(1) yy(end)])
xlabel('x / \lambda')
ylabel('y / \lambda')
title(['z = ' num2str(L),' \lambda']);
view([0 90])

toc
%%
lb = 532e-9; % operating wavelength

xmax = 0.5*(Nx-1)*8e-6/lb;
ymax = 0.5*(Ny-1)*8e-6/lb;

xxa = linspace(-xmax,xmax,Nx);
yya = linspace(-ymax,ymax,Ny);

fsc = ymax/max(yy);
Zpf = (Zdr)*fsc^2;
[XXa,YYa] = meshgrid(xxa,yya);

tic
[Det,Energy,Init_Field] =  GenerateHologram(U,Nx,Ny,XXa*lb,YYa*lb);

% Normalize the phase function  
SLM1=Det/max(Det(:))*255;
H1 = SLM1;
% The aperture size of the SLM is 1920 X 1200
S1 = zeros(1200,1920);
% Assign the hologram to the center of the SLM
S1(600-floor(Ny/2)+1:600+floor(Ny/2),960-floor(Nx/2):960+floor(Nx/2)-1) = H1;

% Digitize the hologram into 256 levels 
H = uint8(((round(S1*255./max(max(S1))))));
% Save the hologram into png file to be addressed onto SLM 

imwrite(H,('holog_cylinder.png'),'png') 
% load this png file in the script 'SLM.m' to create a .csv file
% to be uploaded into the SLM
toc
%%
% complex field at the back focal plane of the first lens
figure
subplot(2,2,1)
surf(XXa*lb*1e3,YYa*lb*1e3,abs(Init_Field).^2);grid
shading interp
colormap(map)
colorbar 
xlabel('x (mm)')
ylabel('y (mm)')
subtitle('Back focal plane of the first lens');
view([0 90])
axis square
xlim([min(xxa)*1e3*lb max(xxa)*1e3*lb])
ylim([min(yya)*1e3*lb max(yya)*1e3*lb])

% Fourier plane - iris plane
nx = Nx;
ny = Ny;
F1 = ifftshift(fft2(fftshift(exp(1i*Det))));

dfx = 1/(nx*(xxa(2)-xxa(1))*lb);
dfy = 1/(ny*(yya(2)-yya(1))*lb);

fx = (-nx/2:nx/2-1)*dfx*2*pi;
fy = (-ny/2:ny/2-1)*dfy*2*pi;

[Fx,Fy] = meshgrid(fx,fy);

subplot(2,2,2)
surf(Fx,Fy,abs(F1).^2);grid
shading interp
colormap(map)
colorbar 
xlabel('k_x')
ylabel('k_y')
subtitle('Fourier plane');
view([0 90])
set(gca, 'CLim', [0, 1e6]);
axis square
xlim([min(fx) max(fx)])
ylim([min(fy) max(fy)])

% selecting the first diffraction order in the Fourier plane
gy = 2*pi*(Nx/2)/(2*Nx*8e-6);
gx = 2*pi*(Ny/2)/(2*Ny*8e-6);

rdy = max(fy)-gy;
rdx = (max(fx)-gx)/2;
gxp = (gx/max(fx))*0.5*nx;
gyp = (gy/max(fy))*0.5*ny;
rdyp = (rdy/max(fy))*0.5*ny;
rdxp = (rdx/max(fx))*0.5*nx;

xc = round(0.5*nx+gxp);
yc = round(0.5*ny+gyp);

F2 = F1(floor(yc-rdyp):floor(yc+rdyp),floor(xc-rdxp):floor(xc+rdxp));
F2b = padarray(F2,[floor((ny-size(F2,1))*0.5)+1,floor((nx-size(F2,2))*0.5)+1],0,'both');
F2a = imresize(F2b,[ny nx]);

subplot(2,2,3)
surf(Fx,Fy,abs(F2a).^2)
shading interp
colormap(map)
colorbar 
view([0 90])
set(gca, 'CLim', [0, 1e6]);
subtitle('Selecting First diffraction order');
axis square
xlim([min(fx) max(fx)])
ylim([min(fy) max(fy)])
xlabel('k_x')
ylabel('k_y')

% complex field at the front focal plane of the second lens
F3 = ifftshift(ifft2(fftshift(F2a)));

subplot(2,2,4)
surf(XXa*lb*1e3,YYa*lb*1e3,abs(F3).^2);grid
shading interp
colormap(map)
colorbar 
xlabel('x (mm)')
ylabel('y (mm)')
subtitle('Front focal plane of the second lens');
view([0 90])
axis square
xlim([min(xxa)*1e3*lb max(xxa)*1e3*lb])
ylim([min(yya)*1e3*lb max(yya)*1e3*lb])
%%