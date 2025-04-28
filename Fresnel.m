function [U,u,v]=Fresnel(U0,lambda,z,dx,dy,t)


% U0 = padarray(U0,[0,0],0, 'both');

[ny, nx] = size(U0); 
c = 3*10^8;
k= 2*pi/(lambda);
omega = 2*pi*c/lambda;
Lx = dx * (nx);
Ly = dy * (ny);

dfx = 1./Lx;
dfy = 1./Ly;

u = ones(ny,1)*((1:nx)-nx/2-1)*dfx;    
v = ((1:ny)-ny/2-1)'*ones(1,nx)*dfy;   

% size(u)
% size(v)

O = fftshift(fft2(U0));

%  size(u)
%  size(v)
%  imagesc(u)
% figure
%  imagesc(v)
H = exp(1i*k*z).*exp(1i*pi*(lambda)*z*(u.^2+v.^2));  

U = ifft2(ifftshift(O.*H))*exp(1i*omega*t); 