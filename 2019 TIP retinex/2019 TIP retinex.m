%Code for model : min_{R, L} lambda_1|gra^alpha R|||_1 + ...
%                            lambda_2||gra^betaL||_1 + 1/2*||RL - I||_2^2

close all
clear all

lambda1 = 0.01; lambda2 = 0.1; mu1 = 0.02; mu2 = 10;
gamma1 = 0.02; gamma2 = 10;
alpha = 1.3; beta = 1.7;

%S = double(imread('cameraman.tif'));
%S = double(imread('3.png'));
%S = double(imread('girl.jpg'));
%S = double(imread('5.png'));
%S = double(imread('59078.jpg'));
%S = double(imread('darkness.jpg'));
%S = double(imread('N-014-0.bmp'));
%S = double(imread('N-041-0.bmp'));
%S = double(imread('shadow3.bmp'));
%S = double(imread('N-025-0.bmp'));
%S = double(imread('N-003-0.bmp'));
%S = double(imread('N-002-0.bmp'));
%S = double(imread('N-051-0.bmp'));
%S = double(imread('YouTube_0062.png'));
%S = double(imread('HK_0136.png'));
%S = double(imread('10.bmp'));
%S = double(imread('9.bmp'));
%S = double(imread('7.bmp'));
%S = double(imread('5.jpg'));
%S = double(imread('cloudy (16).jpg'));
%S = double(imread('59078.jpg'));
%S = double(imread('12003.jpg'));
%S = double(imread('cloudy (16).jpg'));
%S = double(imread('35010.jpg'));
%S = double(imread('277095.jpg'));
%S = double(imread('69020.jpg'));%daishu
%S = double(imread('12003.jpg')); 
%S = double(imread('H-084-0.bmp'));
%S = double(imread('forest_recovered.jpg'));

%space transfer
H = rgb2hsv(S);
I = H(:,:,3);
[rows, cols] = size(S(:,:,1));

L = I; R = zeros(rows,cols);
d1x = zeros(rows, cols); d1y = zeros(rows, cols);
b1x = zeros(rows, cols); b1y = zeros(rows, cols);
b2x = zeros(rows, cols); b2y = zeros(rows, cols);
p1 = zeros(rows, cols); p2 = zeros(rows, cols);
u = d1x; v = d1x; 
iter = 0;

for i = 1:25 %25%while max(err_R, err_L)>1e-4 % 
    iter = iter + 1;
    R_n = R; L_n = L;
    
    % for d1 d2
    N1 = nabla_a(rows, alpha);
    M1 = nabla_a(cols, alpha);
    [d1x, d1y] = shrink2anisotropic(N1*u + b1x, u*M1' + b1y,lambda1/mu1);
    
    N2 = nabla_a(rows, beta);
    M2 = nabla_a(cols, beta);
    [d2x, d2y] = shrink2anisotropic(N2*v + b2x, v*M2' + b2y,lambda2/mu2);
    
    % for R L
    R = [I.*L + gamma1*(u - p1)]./(L.*L + gamma1);
    %[ss, R] = BM3D(1,R,10);
    R = max(1e-4, R);
    R = min(R, 1);
    
    L = [I.*R + gamma2*(v - p2)]./(R.*R + gamma2);
    L = max(I, L);
    
    % for u v
    A1 = N1'*N1 + gamma1/mu1*eye(rows);
    B1 = M1'*M1;
    C1 = N1'*(d1x - b1x) + (d1y - b1y)*M1 + gamma1/mu1*(R + p1);
    u = sylvester(A1, B1, C1);
    
    A2 = N2'*N2 + gamma2/mu2*eye(rows);
    B2 = M2'*M2;
    C2 = N2'*(d2x - b2x) + (d2y - b2y)*M2 + gamma2/mu2*(L + p2);
    v = sylvester(A2, B2, C2);
    
    % for b1 b2 p1 p2
    b1x = b1x + N1*u - d1x; b1y = b1y + u*M1' - d1y;
    
    b2x = b2x + N2*v - d2x; b2y = b2y + v*M2' - d2y;
    
    p1 = p1 + R - u; p2 = p2 + L - v;
     
end

%gamma correction
gamma_1 = 4.6;gamma_2 = 2.6;
R_gamma = (R).^(1/gg);
L_gamma = 255*((L./255).^(1/hh));
enhanced_V = (L_gamma.*R_gamma);
H(:,:,3) = enhanced_V;
 enhanced_result = hsv2rgb(H);


function [xs,ys] = shrink2anisotropic(x,y,lambda)
xs = max(0, abs(x) - lambda).*sign(x);
ys = max(0, abs(y) - lambda).*sign(y);
end

function M = nabla_a(n, alpha)
M = zeros(n, n);
cak = zeros(1, n);
cak(1,1) = 1;
for ii = 2:n
    cak(1, ii) = (1 - (alpha + 1)/(ii-1))*cak(1, ii-1);
end

for jj = 1:n
    M(jj, 1:jj) = fliplr(cak(1,1:jj));
end

end
