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
zT = 0;

xT2 = 0.1;
yT2 = -0.05;
zT2 = 0;

xF = [xT xT2];
yF = [yT yT2];
zF = [zT zT2];

% Receive array horizontal
dxR = c/f(end);            % receive antenna array sampling
nbR = 51;                   % number of antennas
xR = (1:nbR)*dxR;            % receive array coordinates
xR = xR - mean(xR);          % centered on x = 0;
yR = xR .*0;                 % receive array on the plane y = 0;

% Receive array vertical
dzR = c/f(end);            % receive antenna array sampling
nbR = 51;                   % number of antennas
zR = (1:nbR)*dzR;            % receive array coordinates
zR = zR - mean(zR);          % centered on x = 0;
yR = zR .*0;                 % receive array on the plane y = 0;

zees = zeros(numel(xR));

% Circle target

xCenter = 0;
yCenter = 0.9;
zCenter = 0;
theta = 0 : 0.01 : 2*pi;
radius = 0.2;
xa = radius * cos(theta) + xCenter;
za = radius * sin(theta) + zCenter;
ya = zeros(numel(xa)) + yCenter;
SigC = ones(numel(xa));

S = zeros(nbR, numel(f));
for mC = 1:numel(xa)
    rT = sqrt((xF(1)-xa(mC)).^2 + (yF(1)-ya(mC)).^2 + (zF(1)-za(mC)).^2); % distance from transmit antennas to target(s)
    rR = sqrt((xR-xa(mC)).^2 + (yR-ya(mC)).^2 + (zR-za(mC)).^2); % distance from targets to the receive antenna
    for mf = 1:numel(f)
        S(:,mf) = S(:,mf) + (exp(-1j*k(mf)*rT)./sqrt(rT) .* SigC(mC) .* exp(-1j*k(mf)*rR)./sqrt(rR)).'; % interaction of the transmit waves with the target(s) and measurement
    end
end

s = (ifft(ifftshift(S,2),[],2));

% Reconstruction grid (should be defined according to the resolution)
x = linspace(-0.4,0.4,41);
y = linspace(0.1,1.1,41);
z = linspace(-0.4,0.4,41);

% Computation of the transfer (Green) matrix (horizontal)
[XR,F,X,Y] = ndgrid(single(xR),single(f),single(x),single(y));
G = exp(-1j*(2*pi*F/c).* ((sqrt((xF(1)-X).^2 + (yF(1)-Y).^2)+(sqrt((XR-X).^2 + (yR(1)-Y).^2)))));
G = reshape(G,[numel(xR)*numel(f),numel(x)*numel(y)]);

% Computation of the transfer (Green) matrix (vertical)
[ZR,F,Y,Z] = ndgrid(single(zR),single(f),single(y),single(z));
G_d = exp(-1j*(2*pi*F/c).* ((sqrt((zF(1)-Z).^2 + (yF(1)-Y).^2)+(sqrt((ZR-Z).^2 + (yR(1)-Y).^2)))));
G_d = reshape(G_d,[numel(zR)*numel(f),numel(z)*numel(y)]);

% Reconstruction by matched filtering (phase compensation) (horizontal)
Im = G'*S(:);
Im = reshape(Im,[numel(x),numel(y)]); %from vector to matrix format

% Reconstruction by matched filtering (phase compensation) (vertical)
Im_d = G_d'*S(:);
Im_d = reshape(Im_d,[numel(z),numel(y)]); %from vector to matrix format

figure;
plot3(xF(1),yF(1),0,'kv')
hold on
plot3(xF(2),yF(2),0,'kv')
plot3(xR,yR,zees,'rv')
plot3(zees,yR,zR,'rv')
plot3(xa, ya, za)
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
subplot(2,3,1)
pcolor(xR,f/1e9,abs(S).')
shading flat
xlabel('x (m)')
ylabel('f (GHz)')
title('Magnitude of S(f)')

subplot(2,3,2)
pcolor(xR,f/1e9,angle(S).')
shading flat
xlabel('x (m)')
ylabel('f (GHz)')
title('Phase of S(f)')

subplot(2,3,3)
pcolor(xR,t*1e9,abs(s).')
shading flat
xlabel('x (m)')
ylabel('t (ns)')
title('Magnitude of s(t)')

subplot(2,3,4)
pcolor(x,y,abs(Im).')
shading flat
xlabel('x (m)')
ylabel('y (m)')
title('Reconstructed image')

subplot(2,3,5)
pcolor(y,z,abs(Im_d).')
shading flat
xlabel('z (m)')
ylabel('y (m)')
title('Reconstructed image')

figure;
plot3(xa, ya, za)
grid
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
