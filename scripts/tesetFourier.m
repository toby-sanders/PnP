clear;
SNR = 5;
iter = 20;
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
pp = .2;

I = im2double(rgb2gray(imread([path,'lena.png'])));

[d1,d2] = size(I);
S = rand(d1,d2);
[~,S] = sort(S(:));
S = S(1:round(d1*d2*pp));
S = sort(S);
if S(1)~=1
    S = [1;S];
end
Sp = 1:d1*d2;
Sp(S) = '';

b = fft2(I);
b = b(S);
b = add_Wnoise(b,SNR);

F = @(x,mode)fft_partial_transforms(x,mode,[d1,d2],S);


popts.iter = 25;
popts.tol = 1e-4;
% Reg = getBM3D(C*sigma/d);
popts.sigma0 = 5/255;
popts.mu = 20;
popts.nonneg = true;
Reg = getLTBM3D;

gpuID = gpuDevice(2);
reset(gpuID);
[p3rec,p2out] = PnP3_prox(F,Reg,b,[d1,d2,1],popts);
%%
popts.mu = 20;
popts.sigma0 = 5/255*sqrt(25);
[p2rec,p3out] = PnP3(F,Reg,b,[d1,d2,1],popts);
[p4rec,p4out] = PnP3_Fourier(S,Reg,b,[d1,d2,1],popts);
%%
bstr = zeros(d1,d2);
bstr(S) = b;
bprec = real(ifft2(bstr));
sZ = 32; % size of blocks
LTopts.blockSize = sZ;
LTopts.numMax = 16; % max number of blocks to match for each ref. block
LTopts.numMin = 16; % min number of blocks to match for each ref. block
LTopts.wname = 'bior';
LTopts.wnamez = 'db';
LTopts.order = 13;
LTopts.orderz = 1;
LTopts.levels = 3;
LTopts.cycleSpin = 2;
LTopts.matchSpins =  [0 4 sZ/2-1];
LTopts.wiener = false;
LTopts.tauMode = 2;
LTopts.filtType = 'ht';
mySig = 5/255;
for i = 1:iter
    y = real(ifft2(bstr));
    y = LTBM3D(y,mySig,LTopts);
    tmp = fft2(y);
    tmp = tmp(Sp) + randn(1,numel(Sp)) + 1i*randn(1,numel(Sp));
    bstr(Sp) = tmp;
    figure(77);colormap(gray);imagesc(y);title(i);
end





%%

figure(127);colormap(gray);tiledlayout(3,3,'tilespacing','none');
t1 = nexttile;imagesc(p3rec,[0 1]);title('p3 - prox');
t2 = nexttile;imagesc(bprec,[0 1]);title('ifft rec');
t3 = nexttile;imagesc(I,[0 1]);title('original');
t4 = nexttile;imagesc(p2rec,[0 1]);title('P3 - admm');
t5 = nexttile;imagesc(p4rec,[0 1]);title('P3 - fourier');
t6 = nexttile;imagesc(y,[0 1]);
t7 = nexttile;hold off;
semilogy(p3out.rel_chg_inn);hold on;
semilogy(p4out.rel_chg);hold off;
linkaxes([t1 t2 t3 t4 t5 t6]);