clear all
close all
clc

% Parameters
nbf = 101;                  % number of frequency samples
f = linspace(15e9,20e9,nbf); % frequency vector
c = 3e8;                     % speed of light
k = 2*pi*f/c;                % wavenumber
df = f(2)-f(1);              % frequency resolution
t = linspace(0,1/df,nbf);    % time vector

% Transmit antenna
xT = 0;
yT = -0.05;

% Receive array
dxR = c/f(end);            % receive antenna array sampling
nbR = 31;                   % number of antennas
xR = (1:nbR)*dxR;            % receive array coordinates
xR = xR - mean(xR);          % centered on x = 0;
yR = xR .*0;                 % receive array on the plane y = 0;

% Target
xC = 0;
yC = 0.5;
%SigC = ones(numel(x)); % Reflectivity

% Line target

xa = 0:1;
ya = 0.8*xa + 0.8;
za = zeros(numel(xa));
SigC = ones(numel(xa));
% for loops are used to facilitate understanding. It's a very bad idea if you try to optimize computation times

S = zeros(nbR, numel(f));
for mC = 1:numel(xa)
    rT = sqrt((xT-xa(mC)).^2 + (yT-ya(mC)).^2); % distance from transmit antennas to target(s)
    rR = sqrt((xR-xa(mC)).^2 + (yR-ya(mC)).^2); % distance from targets to the receive antenna
    for mf = 1:numel(f)
        S(:,mf) = S(:,mf) + (exp(-1j*k(mf)*rT)./sqrt(rT) .* SigC(mC) .* exp(-1j*k(mf)*rR)./sqrt(rR)).'; % interaction of the transmit waves with the target(s) and measurement
    end
end

s = (ifft(ifftshift(S,2),[],2));

% Reconstruction grid (should be defined according to the resolution)
x = linspace(-0.4,0.4,41);
y = linspace(0.1,1.1,41);

% Computation of the transfer (Green) matrix
[XR,F,X,Y] = ndgrid(single(xR),single(f),single(x),single(y));
G = exp(-1j*(2*pi*F/c).* ((sqrt((xT-X).^2 + (yT-Y).^2)+(sqrt((XR-X).^2 + (yR(1)-Y).^2)))));
G = reshape(G,[numel(xR)*numel(f),numel(x)*numel(y)]);

% Reconstruction by matched filtering (phase compensation)
Im = G'*S(:);
Im = reshape(Im,[numel(x),numel(y)]); %from vector to matrix format

%fig1 = figure(1); clf()
%fig1.Position = [680 542 871 420];
%subplot(2,4,[1 2 5 6])

figure;
plot3(xT,yT,0,'kv')
hold on
plot(xR,yR,'rv')
scatter3(xC ,yC, 0,'bo')
%plot3([xa xb], [ya yb], [za zb])
plot(xa, ya)
hold off
legend('Transmit antenna','Receive antennas','Target(s)')
grid on
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
daspect([1 1 1])
ylim([-0.1 1])
xlim([-0.4 0.4])
zlim([-0.4 0.4])

figure;
subplot(2,2,1)
pcolor(xR,f/1e9,abs(S).')
shading flat
xlabel('x (m)')
ylabel('f (GHz)')
title('Magnitude of S(f)')

subplot(2,2,2)
pcolor(xR,f/1e9,angle(S).')
shading flat
xlabel('x (m)')
ylabel('f (GHz)')
title('Phase of S(f)')

subplot(2,2,3)
pcolor(xR,t*1e9,abs(s).')
shading flat
xlabel('x (m)')
ylabel('t (ns)')
title('Magnitude of s(t)')

subplot(2,2,4)
pcolor(x,y,abs(Im).')
shading flat
xlabel('x (m)')
ylabel('y (m)')
title('Reconstructed image')
