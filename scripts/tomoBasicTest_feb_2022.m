% tomography example of white noise reconstruction in linear solvers for
% post process denoising
clear; 
d = 256; % reconstruction dimension
SNR = 20; % SNR in data
% angles = -89:0.5:90; % tomo angles
angles = -88:1.5:90;
rng(2102);

% load data, generate sinogram, blur and downsample, add noise
% load('C:\Users\toby.sanders\Documents\miscellaneous\phantoms\circles-denoised');
% load('C:\Users\toby.sanders\Documents\miscellaneous\phantoms\bone');
%  load('/Users/tobysanders/Dropbox/archives/data/tomography/phantoms/circles-denoised');
% load('C:\Users\toby.sanders\Documents\miscellaneous\phantoms\trippy-circles');
load('C:\Users\toby.sanders\Dropbox\archives\data\tomography\phantoms\circles-denoised');
% X = Z;
% load('C:\Users\toby.sanders\Documents\miscellaneous\phantoms\gamecock');
% load('/Users/tobysanders/Dropbox/archives/data/tomography/phantoms/gamecock');
% load('/Users/tobysanders/Dropbox/archives/data/tomography/phantoms/circles-denoised');
P = imresize(X,[d,d]);
W = radonmatrix(angles,d,d);
W1 = radonmatrix(angles,512,d);
b = W1*X(:);
% b = W*P(:);
sino = reshape(b,d,numel(angles));
[sinoD,sigma] = add_Wnoise(sino,SNR); % add noise

gpuID = gpuDevice(2);
reset(gpuID);

hopts.mu = 150;
hopts.order = 1;
hopts.levels = 1;
hopts.nonneg = true;
[UTV,outTV] = HOTV3D(W,sinoD(:),[d,d],hopts);

popts.iter = 35;
popts.tol = 1e-4;
% Reg = getBM3D(C*sigma/d);
popts.sigma0 = 3.0/255;
popts.init = UTV;
Reg = getGBM3D;% getLTBM3D;
[p3rec,p3out] = PnP3_prox(W,Reg,sinoD(:),[d,d,1],popts);


%%

figure(890);colormap(gray)
tiledlayout(2,3,'tilespacing','compact');
t1 = nexttile;imagesc(P,[0 1]);title('original image');
t2 = nexttile;imagesc(UTV,[0 1]);
title(sprintf('TV solution, PSNR = %g',myPSNR(P,UTV,1)));
t3 = nexttile;imagesc(p3rec,[0 1]);
title(sprintf('P3 solution, PSNR = %g',myPSNR(P,p3rec,1)));
linkaxes([t1 t2 t3]);
set([t1 t2 t3],'fontsize',16);

t1 = nexttile;imagesc(P,[0 1]);
t2 = nexttile;imagesc(UTV,[0 1]);
t3 = nexttile;imagesc(p3rec,[0 1]);
linkaxes([t1 t2 t3]);
set([t1 t2 t3],'fontsize',16);


