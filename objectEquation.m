%Solving the discrete system (for inverse scattering problem)
clear;
close all;

%% Defining inputs
%Using si = [cb/c(rho)]^2 - 1 as the contrast function
%Solving the inverse scattering problem

syms lamB cb
%lamB = 0.03;
%cb = 2.8e8;
% Number of points
n = 256;

%Position of the delta source
a = -lamB;
%z = -a;
alpha = 1;
c = (1/8).*cb;
%c = (0.95).*cb;

%ndelZ = d
%d = alpha.*lamB
d = alpha.*lamB;
delZ = d./n;
delZ = simplify(delZ);

%Grid nodes
k = 1:1:n;
zk = (k - 1/2).*delZ;
zk = simplify(zk);

%for loop to reduce the computation
%sub = delZ./lamB
substitute = alpha./n;
Zsub = (k - 1/2);

%Contrast funtion
xiZ = (1 - (cb./c).^2).*ones(size(zk));
xiZ = double(simplify(xiZ));

%Constrast matrix
C = diag(xiZ);

%uinc matrix
%fb = cb./abs(lamB);
kb = 2*pi./abs(lamB);
gamma = 1j.*kb;
uinc = exp(-gamma.*abs(zk+a));
uinc = double(simplify(uinc));

%G matrix
const = -1j.*pi.*substitute;
%const = double(simplify(const));
G = zeros(n);

%%
for ind = 1:n
    %adding abs to denom to make it easy to simplify
    vec1 = exp(-1j.*2.*pi.*abs(Zsub(ind) - Zsub).*substitute);
    %vec1 = double(simplify(vec1));
    G(ind, :) = vec1.*const;
end

%What is required?
%Plot of the solution and of the incident field

%Identity matrix
I = eye(n);

%L.H.S of system of equation
A = I + G*C;

%R.H.S of system of equation
B = uinc';

%Solving the system of equation
U = linsolve(A,B);

%% Plotting
plotVec = double(simplify(zk./lamB));

figure();
plot(plotVec, abs(U), 'LineWidth', 1.5);
xlabel('Normalized range 0<z<d');
ylabel('Magnitude of the solution');
title('Magnitute of the solution w.r.t. z');
grid on;

figure();
plot(plotVec, abs(uinc), 'LineWidth', 1.5);
xlabel('Normalized range 0<z<d');
ylabel('Magnitude of the incident field');
title('Magnitute of the incident field w.r.t. z');
grid on;