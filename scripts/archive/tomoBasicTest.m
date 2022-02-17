% tomography example of white noise reconstruction in linear solvers for
% post process denoising
clear; 
d = 256; % reconstruction dimension
SNR = 20; % SNR in data
% angles = -89:0.5:90; % tomo angles
angles = -88:4:90;
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

popts.iter = 75;
popts.tol = 1e-4;
% Reg = getBM3D(C*sigma/d);
popts.sigma0 = 1/255;
Reg = getGBM3D;% getLTBM3D;
[p3rec,p3out] = PnP3_prox(W,Reg,b,[d,d,1],popts);

popts.mu = 4000/sigma^2;
popts.beta = 25;
popts.sigma0 = 10/255;
R = getLTBM3D;
[p2rec,p2out] = PnP3(W,R,b,[d,d,1],popts);

popts.beta = popts.beta/2;
popts.mu = popts.mu;
[p4rec,p4out] = PnP3(W,R,b,[d,d,1],popts);

%%
myrel(p2rec,p4rec)

figure(890);colormap(gray)
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(P);
t2 = nexttile;imagesc(p2rec);
t3 = nexttile;imagesc(p3rec);
t4 = nexttile;imagesc(p4rec);
linkaxes([t1 t2 t3 t4]);

figure(891);hold off;
semilogy(p2out.rel_chg_inn);hold on;
semilogy(p4out.rel_chg_inn);
semilogy(p3out.rel_chg);
hold off;

