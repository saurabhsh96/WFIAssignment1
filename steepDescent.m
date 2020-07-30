%Steepest descent algorithm
clear;
close all;

%% Steepest descent
syms x1 x2 tau

a = 16;
f = @(x1, x2) 1./2.*x1.^2 + a./2.*x2.^2;
g = gradient(f, [x1, x2]);

%initial guess
[X, Y] = meshgrid(-2:.01:2,-2:.01:2);
G1 = subs(g(1), [x1 x2], {X,Y});
G2 = subs(g(2), [x1 x2], {X,Y});
%quiver(X, Y, G1, G2)

%Initial guess
x0 = [1.8 1.4];

%Stepwise vectors, first value in x0
x(1) = x0(1); 
y(1) = x0(2);

%Initial value of f
nextVal = f(x(1), y(1));

%Defining xd and fd
xd = @(x, d, s) x - d.*s;
fd = @(x, y, d, s) f(xd(x, d, s(1)), xd(y, d, s(2))); 

%for loop to calculate the vectors
k = 2; 

%figure();
while nextVal > 10^-6
    %Gradient calculation
    s = double(subs(g, [x1 x2], {x(k-1),y(k-1)}));
    
    if(s==0)
        break;
    end
    
    fdG = fd(x(k-1), y(k-1), tau, s);
    %tauStar = fminsearch(fdG, 0);
    fdG = symfun(fdG, tau);
    fdG(tau) = simplify(fdG);
    
    %differentialtion of fdG
    fdG_Dash = diff(fdG);
    
    %Solving for minimum
    tauStar = solve(fdG_Dash, tau);

    %Getting the tauStar
    tauStar = double(tauStar);
    im = imag(tauStar);
    re = real(tauStar);
    ri = im./re;
    
    tauStar(ri>1e-3) = [];
    tauStar = real(tauStar);
   
    temp = double(fdG(tauStar));
    [~, te] = min(abs(temp));
    tauStar = tauStar(te);
    
    %prevStep = [xi;
    
    x(k) = xd(x(k-1), tauStar, s(1));
    y(k) = xd(y(k-1), tauStar, s(2));
    
    nextVal = f(x(k), y(k));
    %plot(x(k), y(k), 'o'); hold on;
    
    k = k + 1;
end

%% Plotting
figure();
z = f(X, Y); %Level Curves
contour(X,Y,z,15); hold on;

%Plotting level curves and steepese descent points
plot(x, y, '-o', 'LineWidth', 1.5);
grid on;
xlabel('x1', 'FontWeight', 'bold');
ylabel('x2', 'FontWeight', 'bold');
title('Level curves and steepest descent points');