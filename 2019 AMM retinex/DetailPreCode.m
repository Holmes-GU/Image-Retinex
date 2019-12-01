%Code for 'A Detail Preserving Variational Moadel For Image Retinex'
%The model:
%   min_(R, L) 1/2 ||RL - I||_2^2 + lambda1||grad R||_1 ...
%              + lambda2/2||grad L||_2^2


clear all
close all
%S = double(imread('6.png'));
%S = double(imread('3.png'));
%S = double(imread('5.png'));
S = double(imread('girl.jpg'));

lambda1 = 0.01; lambda2 = 0.1; sigma1 = 5; sigma2 = 5;
sigma3 = 0.02; sigma4 = 5;

[rows, cols] = size(S(:,:,1));
H = rgb2hsv(S);
I = H(:,:,3);
I = max(I, 1e-6); 

L = I;
R = zeros(rows,cols);
t_R = 1; t_L = 1;
dx = zeros(rows, cols);
dy = zeros(rows, cols);
bx = zeros(rows, cols);
by = zeros(rows, cols);
b1 = bx; b2 = bx; px = bx; py = bx;
otfDx = psf2otf([1,-1],[rows, cols]);
otfDy = psf2otf([1;-1],[rows, cols]);

while t_R>(1e-3)&& t_L>(1e-3) 
   R_n = R;
    L_n = L;
    qx = sigma4*(Dx(L) + px)/(lambda2 + sigma4);
    qy = sigma4*(Dy(L) + py)/(lambda2 + sigma4);
    
    Nomin1 = sigma2*fft2(L + b2) + sigma4*(otfDx.*fft2(qx -px )+otfDy.*fft2(qy - py));
    Denom1 = sigma4*(abs(otfDx).^2 + abs(otfDy).^2)+sigma2;
    Nomin2 = sigma1*fft2(R + b1) + sigma3*(otfDx.*fft2(dx -bx )+otfDy.*fft2(dy - by));
    Denom2 = sigma3*(abs(otfDx).^2 + abs(otfDy ).^2)+sigma1;
    v = real(ifft2(Nomin1./Denom1));
    u = real(ifft2(Nomin2./Denom2));
    
    [dx, dy] = shrink2anisotropic(Dx(u)+bx, Dy(u)+by , lambda1/sigma3);
    
    R2 = [L.*I + sigma1.*(u - b1)]./(L.*L + sigma1);
    R = min(R2,1); R = max(1e-4,R2);
    L = [I.*R + sigma2.*(v - b2)]./(R.*R + sigma2);
    L = max(I,L);
   bx = bx + Dx(R) - dx;
   by = by + Dy(R) - dy;
   b1 = b1 + R - u;
   b2 = b2 + L - v;
   px = px + Dx(v) - qx;
   py = py + Dx(v) - qy;
   
   t_R = norm(R - R_n,'fro')/norm(R_n,'fro');
   t_L = norm(L - L_n,'fro')/norm(L_n,'fro');
    
end
figure, imshow(R)

load modelparameters.mat
 
 blocksizerow    = 96;
 blocksizecol    = 96;
 blockrowoverlap = 0;
 blockcoloverlap = 0;
 
%R1 = max(R, 1e-6);
R_gamma = (R).^(1/2.3);
L_gamma = 255*((L./255).^(1/2.8));
enhanced_V = (L_gamma.*R_gamma);
enhanced_V = 255*(enhanced_V - min(min(enhanced_V)))./(max(max(enhanced_V)) - min(min(enhanced_V)));
H(:,:,3) = enhanced_V;
enhanced_result = hsv2rgb(H);

quality = computequality(uint8(enhanced_result),blocksizerow,blocksizecol,blockrowoverlap,blockcoloverlap, ...
    mu_prisparam,cov_prisparam)






function d = Dx(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
end


function d = Dy(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
end


function [xs,ys] = shrink2anisotropic(x,y,lambda)
xs = max(0, abs(x)-lambda).*sign(x);
ys = max(0, abs(y)-lambda).*sign(y);
end




