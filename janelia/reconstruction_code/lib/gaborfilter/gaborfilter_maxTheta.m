function [g_1, g_2] = gaborfilter_maxTheta(I, Sx, Sy, f, n_theta)
% [g_1, g_2] = gaborfilter_maxTheta(I, Sx, Sy, f, n_theta)
% Gives maximum Gabor filter response over all theta.
% Code borrowed from gaborfilter2.m
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% The Gabor filter is basically a Gaussian (with variances sx and sy along x and y-axes respectively)
% modulated by a complex sinusoid (with centre frequencies U and V along x and y-axes respectively)
% described by the following equation
%
%               1                -1     x  ^    y  ^
% Gi(x,y) = ---------- * exp ([----{(----) 2+(----) 2}])*Mi(x,y,f);
%            2*pi*sx*sy           2    sx       sy
% i =1,2
% M1(x,y,f) = cos[2*pi*f*sqrt(x^2+y^2)];
% M2(x,y,f) = cos[2*pi*f*(x*cos(theta) + y*sin(theta)];
%
% Describtion :
%
% I : Input image
% Sx & Sy : Variances along x and y-axes respectively
% f : The frequency of the sinusoidal function
% theta : The orientation of Gabor filter
%
% G1 & G2 : The output filters as described above
% gabout1 & gabout2 : The output filtered images
%
%  Author : Ahmad poursaberi  e-mail : a.poursaberi@ece.ut.ac.ir
%          Faulty of Engineering, Electrical&Computer Department,Tehran
%          University,Iran,June 2004
%

if isa(I,'double')~=1
  I = double(I);
end

G1 = [];
for x = -fix(Sx):fix(Sx)
  for y = -fix(Sy):fix(Sy)
    M1 = cos(2*pi*f*sqrt(x^2+y^2));
    G1(fix(Sx)+x+1,fix(Sy)+y+1) = (1/(2*pi*Sx*Sy)) * exp(-.5*((x/Sx)^2+(y/Sy)^2))*M1;
  end
end
Imgabout1 = conv2(I,double(imag(G1)),'same');
Regabout1 = conv2(I,double(real(G1)),'same');
g_1 = sqrt(Imgabout1.*Imgabout1 + Regabout1.*Regabout1);

g_2 = zeros(size(I));
for i = 1:n_theta
  theta = pi*(i-1)/n_theta;
  
  G2 = [];
  for x = -fix(Sx):fix(Sx)
    for y = -fix(Sy):fix(Sy)
      M2 = cos(2*pi*f*(x*cos(theta)+y*sin(theta)));
      G2(fix(Sx)+x+1,fix(Sy)+y+1) = (1/(2*pi*Sx*Sy)) * exp(-.5*((x/Sx)^2+(y/Sy)^2))*M2;
    end
  end
  
  Imgabout2 = conv2(I,double(imag(G2)),'same');
  Regabout2 = conv2(I,double(real(G2)),'same');

  g_2 = max(g_2, sqrt(Imgabout2.*Imgabout2 + Regabout2.*Regabout2));
end

return
end