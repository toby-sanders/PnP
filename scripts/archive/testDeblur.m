clear;
d = 256;
SNR = 100;
omega = 2;
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
% I = im2double(imread('cameraman.tif'));
% I = im2double(rgb2gray(imread('peppers.png')));
I = im2double(rgb2gray(imread([path,'AK.jpg'])));
I = im2double((imread([path,'surfer.png'])));

[d1,d2] = size(I);
[h,hhat] = makeGausPSF([d1,d2],omega);
b = ifft2(fft2(I).*hhat);
[b,sigma] = add_Wnoise(b,SNR);


opts.nonneg = true;
opts.inner_iter = 2;
opts.iter =10;
opts.tol = 1e-4;
opts.mu = 2.4/sigma;% 500;
opts.sigma0 = 8*sigma;
R = getLTBM3D;
R2 = getBM3D;
gpuID = gpuDevice(2);
reset(gpuID);
    
[rec1,out1] = PnP3_deblur(h,R,b,opts);
[rec2,out2] = PnP3_deblur2(h,R,b,opts);
%%
[rec3,out3] = PnP3_deblur(h,R2,b,opts);
LTopts.blockSize = 32;
rec4 = LTBM3D_deconv(b,hhat,sigma,LTopts);
rec5 = BM3DDEB(b,sigma,fftshift(h));
%%
myrel(rec1,rec2)
fprintf('LTBM3D P3 error: %g\n',myrel(rec1,I));
% myrel(rec2,I)
fprintf('BM3D P3 error: %g\n',myrel(rec3,I));
fprintf('LTBM3D colored: %g\n',myrel(rec4,I));
fprintf('BM3D colored: %g\n',myrel(rec5,I));

figure(344);colormap(gray);
tiledlayout(2,3);
t1 = nexttile;imagesc(I,[0 1]);title('original');
t2 = nexttile;imagesc(b,[0 1]);title('blurry');
t3 = nexttile;imagesc(rec1,[0 1]);title('LTBM3D P3');
t4 = nexttile;imagesc(rec3,[0 1]);title('BM3D P3');
t5 = nexttile;imagesc(rec4,[0 1]);title('LT BM3D colored');
t6 = nexttile;imagesc(rec5,[0 1]);title('BM3D colored');
linkaxes([t1 t2 t3 t4 t5 t6]);

figure(88);hold off;
semilogy(out1.rel_chg);hold on;
semilogy(out2.rel_chg_inn);hold off;