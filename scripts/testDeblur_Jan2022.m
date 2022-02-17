clear;
d = 256;
SNR = 20;
omega = 1.5;
K = 2;
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
% I = im2double(imread('cameraman.tif'));
I = im2double(rgb2gray(imread('peppers.png')));
% I = im2double(rgb2gray(imread([path,'AK.jpg'])));
% I = im2double(rgb2gray((imread([path,'lena.png']))));

[d1,d2] = size(I);

g = zeros(d1,d2);
g([1:K],[1:K]) = 1/K^2;
g = fraccircshift(g,[-K/2 + 1/2, -K/2 + 1/2]);
ghat = fft2(g);



[h,hhat] = makeGausPSF([d1,d2],omega);
b = ifft2(fft2(I).*hhat.*ghat);
b = b(1:K:end,1:K:end);
[b,sigma] = add_Wnoise(b,SNR);


opts.nonneg = true;
opts.iter = 20;
opts.tol = 1e-4;
opts.sigma0 = 1.9/255;

% gpuID = gpuDevice(2);
% reset(gpuID);
    
% blur kernel for downsampling


R = getGBM3D;% getLTBM3D(LTopts);
A = getSuperResDeblurOpers(d1,d2,K,hhat,ghat);
[rec1,out1] = PnP3_prox(A,R,b,[d1,d2],opts);

%%
hopts.mu = 150;
hopts.order = 1;
hopts.levels = 1;
hopts.nonneg = true;
[UTV,outTV] = HOTV3D(A,b(:),[d1,d2],hopts);

%%
mm1 = 0;
mm2 = 1;
gam = .75;
figure(77);tiledlayout(2,4,'tilespacing','compact');colormap(gray);
t0 = nexttile;imagesc(I.^gam,[mm1 mm2].^gam);title('original image');
t1 = nexttile;imagesc(max(b,0).^gam,[mm1 mm2].^gam);title(sprintf('blurry/noisy and low resolution'));
t2 = nexttile;imagesc(rec1.^gam,[mm1 mm2].^gam);title(sprintf('P3 solution, PSNR = %g',myPSNR(I,rec1,1)));
t3 = nexttile;imagesc(UTV.^gam,[mm1 mm2].^gam);title(sprintf('TV solution, PSNR = %g',myPSNR(I,UTV,1)));
linkaxes([t0 t2 t3]);
set([t0 t1 t2 t3],'fontsize',16);
set(t1,'Xtick',50:50:150)
t0 = nexttile;imagesc(I.^gam,[mm1 mm2].^gam);
t1 = nexttile;imagesc(max(b,0).^gam,[mm1 mm2].^gam);
t2 = nexttile;imagesc(rec1.^gam,[mm1 mm2].^gam);
t3 = nexttile;imagesc(UTV.^gam,[mm1 mm2].^gam);
linkaxes([t0 t2 t3]);
set([t0 t1 t2 t3],'fontsize',16);


figure(78);semilogy(out1.rel_chg);