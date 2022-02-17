% a demo for the SAR Gotcha data

% The GOTCHA data needed for this code is provided from the publication:
% Casteel, Curtis H., et al. "A challenge problem for 2D/3D imaging of 
% targets from a volumetric data set in an urban environment." Algorithms 
% for Synthetic Aperture Radar Imagery XIV. Vol. 6568. International
% Society for Optics and Photonics, 2007.

% One also needs to download and set up Fessler's code for the NUFFT: 
% https://web.eecs.umich.edu/~fessler/code/index.html

% Finally, one needs to download and set up my imaging software:
% https://www.toby-sanders.com/software

% Some of the methods used here are given in:
% Sanders, Toby, Anne Gelb, and Rodrigo B. Platte.
% "Composite SAR imaging using sequential joint sparsity." 
% Journal of Computational Physics 338 (2017): 357-370.

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 12/4/2018

% add sar setup directory
addpath([pwd,'/SAR_setup']);

clear;
% base path where the data is located
opt.basePath = 'C:\Users\toby.sanders\Dropbox\archives\data\SAR';

% '/home/toby/Documents/MATLAB/heavy/radar/data/Disk4/DATA';

d = 500;
% Define input data parameters here
opt.min_azimuth  = 1;       % Minimum azimuth angle (degrees)  
opt.max_azimuth = 10;       % Maximum azimuth angle (degrees)
opt.min_bandwidth = 0;      % Minimum usable bandwidth
opt.max_bandwidth = 10^11;  % Maximum usable bandwidth
opt.taper_flag = 0;         % Add a hamming taper for sidelobe control
opt.target = 'pass1';       % What target to image
opt.polarization = 'HH';    % What polarization to image (HH,HV,VV)

% Define image parameters here
opt.range_x = 100;           % Scene extent x (m)
opt.range_y = 100;           % Scene extent y (m)
opt.Nfftx = 2*d;
opt.Nffty = 2*d;
opt.Nx = d;          % Number of pixels in x direction
opt.Ny = d;          % Number of pixels in y direction
opt.x0 = 0;            % Center of image scene in x direction (m)
opt.y0 = 0;            % Center of image scene in y direction (m)

%Imaging method (NUFFT or projection matrix)
%set to either 'nufft' or 
%opt.im_method = 'matrix';
opt.im_method = 'nufft';
%opt.im_method = 'matrix';

% joint sparsity options
opt.angle_ranges = 8;
opt.overlap = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS END HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get operator A and data vector bb from the specified opt parameters
[A,bb,Sout] = get_SAR_GOTCHA(opt);

%% basic nufft reconstructions
Ny = opt.Ny;Nx = opt.Nx;
rb = reshape(A(bb,2),Ny,Nx);
% rBP = bpBasicFarField_gotcha(opt);
%% regularized solutions
pat.mu = 80;
pat.mu0 = pat.mu;
pat.order = 1;
pat.iter = 20;
pat.inner_iter = 2;
pat.isreal = false;
pat.disp = false;
pat.scale_A = false;
pat.wrap_shrink = false;
pat.smooth_phase = false;
pat.phase_angles = angle(rb);
[rec,out] = HOTV3D(A,bb,[Ny,Nx,1],pat);




R = getLTBM3D;

gpuID = gpuDevice(2);
reset(gpuID);


opts.iter = 50;
opts.mu = 50;
opts.inner_iter = 2;
opts.disp = false;
opts.beta = 25;
opts.sigma0 = 45/255;
[p3rec,out] = PnP3_SAR(A,R,bb,[Ny,Nx,1],opts);
%%
R = getLTBM3D;
opts.tol = 1e-50;
opts.iter = 25;
opts.init = rb;
opts.sigma0 = 1.5/255;
[recprox,outprox] = Prox_median_SAR(A,bb,[Ny,Nx,1],opts);


%%
figure(101);tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;
displaySAR(log(abs(rec)));
t2 = nexttile;
displaySAR(max(log(abs(p3rec)),-10));
t3 = nexttile;
displaySAR(max(log(abs(rb)),-10));
t4 = nexttile;
displaySAR(max(log(abs(recprox)),-10));
linkaxes([t1 t2 t3 t4]);
figure(102);hold off;
semilogy(out.rel_chg_inn);hold on;
semilogy(outprox.rel_chg); hold off;
title('conv of P3')



