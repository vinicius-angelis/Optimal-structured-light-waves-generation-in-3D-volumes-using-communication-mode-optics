function [SLM0,Energy,L] =  GenerateHologram(Psi1,n_x,n_y,X,Y)

% Obtain the phase of the complex field
ph=angle(Psi1);
v=abs(Psi1)/max(max(abs(Psi1)));%%Normalize field to unity 
%divide the amplitude values into 800 levels
aux=round(v*800+1);
F = zeros(n_y,n_x);
% Load lookup table to perform the inversion of Eq. 18 of Arrizon et al.
load fx2.mat;  
% 
  for mh=1:n_y
      for nh=1:n_x
          temp=aux(mh,nh);
          F(mh,nh)=fx(temp);                                
      end                                                                    
  end 

%parameters for the added blazed grating to spatially shift
% the desired spectrum in k-space off-axis and filter it  
nx = n_y/2;
ny = n_x/2;
gy=ny/(2*n_x*8e-6);
gx=nx/(2*n_y*8e-6);

%Implementing Eq. 16 of Arrizon et al.
Hol0=F.*sin(ph+2*pi*(X*gx+Y*gy));
%calibrate hologram phase from 0 to 1.17 pi
Hol=Hol0-min(Hol0(:));
SLM0=mod(Hol,1.17*pi);

%Additional parameters for monitoring the diffraction efficiency of the
%hologram (not used here)
L = F.*exp(1i*ph);
Energy = sum(sum(L.*conj(L)));
