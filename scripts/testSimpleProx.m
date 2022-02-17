% testing out ideas for "multistep" proximal gradient methods
% the hypothesis is that if I take K gradient steps on the one part of the
% object function, then perform the proximal step on the other part of the
% objective function, is this equivalent to increasing the weight on the
% first part of the object function by a factor of K?

% similarly, if I take K gradient steps of length tau, then take the
% proximal step with strength tau*K, does this solve the same problem as an
% ordinary proximal gradient method?

clear;
d = 512;
m = d/2;
mu = .0001;
SNR = 1;
k = 1;
K = 20;
iter = 5000;

x = sin(4*pi*linspace(0,1,d)');
A = randn(m,d);
b = A*x;
b = add_Wnoise(b,SNR);

V = my_Fourier_filters(k,1,d,1,1);

opts.scale_A = false;
opts.mu = mu;
opts.order = k;
opts.nonneg = false;
opts.tol = 1e-7;
opts.iter = 200;
[U,out] = Tikhonov(A,b,[d,1,1],opts);

tau = 1;
U1 = zeros(d,1);U2 = U1;U3 = U1;U4 = U1;
Atb = A'*b;

% standard PGM
for i = 1:iter
    g = A*U1;
    g = A'*g - Atb;
    U1 = U1 - tau*mu*g;
    U1 = real(ifft(fft(U1)./(1 + tau*V)));
end

% modified PGM with 2 gradient steps
tau2 = tau/2;
for i = 1:iter
    g = A*U2;
    g = A'*g - Atb;
    U2 = U2 - tau2*mu*g;
    
    g = A*U2;
    g = A'*g - Atb;
    U2 = U2 - tau2*mu*g;
    
    U2 = real(ifft(fft(U2)./(1 + tau*V)));
end

% modified PGM with K gradient steps
tauK = tau/K;
for i = 1:iter
    for j = 1:K
        g = A*U3;
        g = A'*g - Atb;
        U3 = U3 - tauK*mu*g;
    end
    U3 = real(ifft(fft(U3)./(1+tau*V)));
end
  
% modified PGM with 2 gradient steps
for i = 1:iter
    g = A*U4;
    g = A'*g - Atb;
    U4 = U4 - tau*mu*g;
    
    g = A*U4;
    g = A'*g - Atb;
    U4 = U4 - tau*mu*g;
    
    U4 = real(ifft(fft(U4)./(1 + tau*V)));
end
opts.mu = mu*2; % compare U4 with this solution (doubling mu same as taking two grad steps)
[U0,out] = Tikhonov(A,b,[d,1,1],opts);
%%
myrel(U,U1)
myrel(U,U2)
myrel(U,U3)
myrel(U0,U4)
myrel(U,U0)

figure(132);hold off;
plot(x);hold on;
plot(U);
plot(U1);
plot(U2);
plot(U3);
hold off;

figure(134);hold off;
plot(U0);hold on;
plot(U4);hold off;