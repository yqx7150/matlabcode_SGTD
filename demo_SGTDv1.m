%% The Code is created based on the methods described in the following papers: 
%   [1] Qiegen Liu, Jianbo Liu, Pei Dong, Dong Liang. SGTD: Structure gradient and texture decorrelating regularization 
%       for image decomposition, The IEEE International Conference on Computer Vision (ICCV), 2013, pp. 1081-1088. . 
%   Author: Qiegen Liu, Jianbo Liu, Pei Dong and Dong Liang
%   Date  : 06/25/2013
%   Version : 1.0 
%   The code and the algorithm are for non-comercial use only.
%   Copyright 2013, Department of Electronic Information Engineering, Nanchang University.
%   The current version is not optimized.

% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% Please cite the paper when you use th code

clear; 
close all; 
path(path,'ultilies/');

%% generate data -- see nested function below
I = double(imread('198023.jpg'))/255;
Bn = I;
mu = 0.31;
level = .6;

F = Bn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m n  d3] = size(F);
opts = [];
opts = getopts(opts);
opts.maxitr = 45; 

D = @(U) ForwardD3(U);
Dt = @(X,Y) Dive3(X,Y);

% initialization
X = F;
Lam1 = zeros(m,n,d3);
Lam2 = Lam1;
Lam3 = Lam1;
beta1 = opts.beta1;
beta2 = opts.beta2;
gamma = opts.gamma;
print = opts.print;
eigsDtD = abs(psf2otf([1,-1],[m,n])).^2 + abs(psf2otf([1;-1],[m,n])).^2;
Denom = eigsDtD + beta2/beta1;  

% finite diff
[D1X,D2X] = D(X);
KXF = X - F;

out.snr = [];
out.nrmdX = [];
out.relchg = [];

%% Main loop
for ii = 1:opts.maxitr
    V1 = D1X + Lam1/beta1;
    V2 = D2X + Lam2/beta1;
    V3 = KXF + Lam3/beta2;
    sumV12 = sqrt(sum(V1.^2 + V2.^2,d3));
    sumV12(sumV12==0) = 1;
    % ==================    %   Shrinkage Step    % ==================
    V = zeros(m,n,d3);
    for jj = 1:d3;
        V(:,:,jj) = sumV12;
    end
    V = max(V - 1/beta1, 0)./V;
    Y1 = V1.*V;
    Y2 = V2.*V;
    Z = max(abs(V3) - mu/beta2, 0).*sign(V3);
    % ==================    %     X-subprolem    % ==================
    Xp = X;
    Temp = (beta2*Z - Lam3)/beta1; 
    X = Dt(Y1 - Lam1/beta1,Y2 - Lam2/beta1) + Temp + beta2/beta1*F;
    X(:,:,1) = fft2(X(:,:,1))./Denom;
    X(:,:,2) = fft2(X(:,:,2))./Denom;
    X(:,:,3) = fft2(X(:,:,3))./Denom;
    X = real(ifft2(X));
    [D1X,D2X] = D(X);% finite diff.
    KXF = X - F;
    % ==================    %    Update Lam    % ==================
    Lam1 = Lam1 - gamma*beta1*(Y1 - D1X);
    Lam2 = Lam2 - gamma*beta1*(Y2 - D2X);
    Lam3 = Lam3 - gamma*beta2*(Z - KXF);
end
out.sol = X;
out.itr = ii;
out.exit = 'Exist Normally';
if ii == opts.maxitr
    out.exit = 'Maximum iteration reached!';
end
figure(1);subplot(131); imshow(Bn,[]);title(sprintf('TV,Input'),'fontsize',13); 
subplot(132); imshow(out.sol,[]);title(sprintf('TV,Cartoon'),'fontsize',13); 
subplot(133); imshow(Bn-out.sol,[]);title(sprintf('TV,Texture'),'fontsize',13); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 900;TVmu = 0.0021;
tempdiag = beta2/beta1/(mu+beta2);
Denom = eigsDtD + beta2*tempdiag;  
opts.maxitr = 25; 
%% Main loop
for ii = 1:opts.maxitr
    V1 = D1X + Lam1/beta1;
    V2 = D2X + Lam2/beta1;
    V3 = (beta2/(mu+beta2))*(KXF + Lam3/beta2);
    sumV12 = sqrt(sum(V1.^2 + V2.^2,d3));
    sumV12(sumV12==0) = 1;
    % ==================    %   Shrinkage Step    % ==================
    V = zeros(m,n,d3);
    for jj = 1:d3;
        V(:,:,jj) = sumV12;
    end
    V = max(V - (abs(Z)+TVmu)/beta1, 0)./V;
    Y1 = V1.*V;
    Y2 = V2.*V;
    Z = max(abs(V3) - V./(mu+beta2), 0).*sign(V3);
    % ==================    %     X-subprolem    % ==================
    Xp = X;
    Temp = ((mu+beta2)*Z - Lam3)*tempdiag; 
    X = Dt(Y1 - Lam1/beta1,Y2 - Lam2/beta1) + Temp + beta2*tempdiag*F;
    X(:,:,1) = fft2(X(:,:,1))./Denom;
    X(:,:,2) = fft2(X(:,:,2))./Denom;
    X(:,:,3) = fft2(X(:,:,3))./Denom;
    X = real(ifft2(X));
    [D1X,D2X] = D(X);% finite diff.
    KXF = X - F;   
    % ==================    %    Update Lam    % ==================
    Lam1 = Lam1 - gamma*beta1*(Y1 - D1X);
    Lam2 = Lam2 - gamma*beta1*(Y2 - D2X);
    Lam3 = Lam3 - gamma*beta2*(Z - KXF);
end
out.sol = X;
out.itr = ii;
out.exit = 'Exist Normally';
if ii == opts.maxitr
    out.exit = 'Maximum iteration reached!';
end
%% Plot result
figure(2);subplot(131); imshow(Bn,[]);title(sprintf('SGTD,Input'),'fontsize',13); 
subplot(132); imshow(out.sol,[]);title(sprintf('SGTD,Cartoon'),'fontsize',13); 
subplot(133); imshow(Bn-out.sol,[]);title(sprintf('SGTD,Texture'),'fontsize',13); 

