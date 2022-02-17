clear;
d = 256;
SNR = 10;
omega = 1;
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
% I = im2double(imread('cameraman.tif'));
% I = im2double(rgb2gray(imread('peppers.png')));
I = im2double(rgb2gray(imread([path,'AK.jpg'])));
I = im2double(rgb2gray((imread([path,'lena.png']))));

[d1,d2] = size(I);
[h,hhat] = makeGausPSF([d1,d2],omega);
b = ifft2(fft2(I).*hhat);
[b,sigma] = add_Wnoise(b,SNR);


opts.nonneg = true;
opts.inner_iter = 2;
opts.iter =50;
opts.tol = 1e-4;
opts.mu = 2.4/sigma;% 500;
opts.sigma0 = 8*sigma;

% gpuID = gpuDevice(2);
% reset(gpuID);
    
% LTopts.Wiener = true;
% R = getLTBM3D(LTopts);
R = getGBM3D;
[rec1,out1] = PnP3_deblur(h,R,b,opts);

gpuID = gpuDevice(2);
reset(gpuID);

% LTopts.Wiener = false;
% R2 = getLTBM3D(LTopts);
% [rec2,out2] = PnP3_deblur(h,R2,b,opts);

%%
addpath('C:\Users\toby.sanders\Dropbox\TobySharedMATLAB\TobyShared\solvers\DenoisingEngines\BM3D\bm3d');
rec3 = BM3DDEB(b,sigma,fftshift(h));
%%
myrel(rec1,I)
myrel(rec3,I)

figure(77);tiledlayout(2,2,'tilespacing','compact');colormap(gray);
t1 = nexttile;imagesc(b,[0 1]);
t2 = nexttile;imagesc(rec1,[0 1]);
t3 = nexttile;imagesc(rec1 + out1.sigma/out1.beta);
nexttile;semilogy(out1.rel_chg);
linkaxes([t1 t2 t3]);