% Composite Plate
clear all, close all, clc

% Input Geometry (m)
Lx = 0.2;
Ly = 1;
h = 0.0002;
orientation = (pi/180)*[0 90 0 90];

% Input Loading Conditions Nx, Ny, Nxy, Mx, My, Mxy
loading = [0; 0; 0; 0; 0; 0]; %(MPa*mm)
p = 10000; %Pa
DT = [-100 -100]; %Celcius
MC = 0; % Percent Weight of H20

% Input Material Properties (Pa)
E1 = 181000e6;
E2 = 10300e6;
G12 = 7170e6;
v12 = 0.28;
v23 = 0.59;
a1 = -0.7e-6;
a2 = 25e-6;
a3 = 0;
B1 = 0;
B2 = 0;
B3 = 0;

% Ply Strength (MPa)
s1plus = 1500;
s1minus = 1500;
s2plus = 40;
s2minus = 246;
s12 = 68;

% Quadratic Strength Parameters
F1 = (1/s1plus) - (1/s1minus);
F2 = (1/s2plus) - (1/s2minus);
F11 = 1/(s1plus*s1minus);
F22 = 1/(s2plus*s2minus);
F66 = 1/(s12^2);
F12 = -0.5*sqrt(F11*F22);

% Calculate Q (Stiffness) Matrix
v21 = v12*E2/E1;
D = 1-v12*v21;
Q11 = E1/D;
Q22 = E2/D;
Q12 = v12*E2/D;
Q66 = G12;
Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];

% Calculate Q_bar Matrix
for var = 1:length(orientation)
    theta = orientation(var);
    c = cos(theta);
    s = sin(theta);
    Ts(:,:,var) = [c^2 s^2 2*c*s;
    s^2 c^2 -2*c*s;
    -c*s c*s (c.^2-s.^2)];
    Te(:,:,var) = [c^2 s^2 c*s;
    s^2 c^2 -c*s;
    -2*c*s 2*c*s (c.^2-s.^2)];
end

for num = 1:length(orientation)
    Q_bar(:,:,num) = inv(Ts(:,:,num))*Q*Te(:,:,num);
end

% Calculate [A], [B], and [D] Matrices
N = length(orientation)/2;
z = -N*h:h:h*N;
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
for k1 = 1:3
    for k2 = 1:3
        for k3 = 1:length(z)-1
            z0 = z(k3);
            z1 = z(k3+1);
            A(k1,k2) = A(k1,k2) + Q_bar(k1,k2,k3)*h;
            B(k1,k2) = B(k1,k2) + (1/2)*Q_bar(k1,k2,k3)*(z1^2-z0^2);
            D(k1,k2) = D(k1,k2) + (1/3)*Q_bar(k1,k2,k3)*(z1^3-z0^3);
        end
    end
end
ABD = [A B; B D];

% Get D Matrix Components
D11 = D(1,1);
D12 = D(1,2);
D16 = D(1,3);
D22 = D(2,2);
D26 = D(2,3);
D66 = D(3,3);

% Coefficients
a12 = [a1; a2; a3];
for l = 1:length(orientation)
    axy(:,l) = Te(:,:,l)\a12;
end
B12 = [B1; B2; B3];
for ll = 1:length(orientation)
    Bxy(:,ll) = Te(:,:,ll)\B12;
end

% Thermal/Moisture Loading
T0 = [0; 0; 0; 0; 0; 0];
T1 = [0; 0; 0; 0; 0; 0];
DT0 = (DT(1)+DT(2))/2;
DT1 = (DT(2)-DT(1))/(h*length(orientation));
M0 = [0; 0; 0; 0; 0; 0];

for q = 1:length(z)-1
    z0 = z(q);
    z1 = z(q+1);
    T0(1:3) = T0(1:3)+Q_bar(:,:,q)*axy(:,q)*DT0*(z1-z0);
    T0(4:6) = T0(4:6)+0.5*Q_bar(:,:,q)*axy(:,q)*DT0*(z1^2-z0^2);
    T1(1:3) = T1(1:3)+0.5*Q_bar(:,:,q)*axy(:,q)*DT1*(z1^2-z0^2);
    T1(4:6) = T1(4:6)+Q_bar(:,:,q)*axy(:,q)*DT1*(z1^3-z0^3)/3;
    M0(1:3) = M0(1:3)+Q_bar(:,:,q)*Bxy(:,q)*MC*(z1-z0);
    M0(4:6) = M0(4:6)+0.5*Q_bar(:,:,q)*Bxy(:,q)*MC*(z1^2-z0^2);
end
otherload = T0 + T1 + M0;
htepsilon0 = ABD\otherload;
midstrain = ABD\(loading+otherload);

% Stress & Strain through thickness
for u = 1:length(z)-1
    zbottom = z(u);
    ztop = z(u+1);
    zk(2*u-1) = zbottom;
    zk(2*u) = ztop;
    DTbottom = DT0 + DT1*(zbottom);
    DTtop = DT0 + DT1*(ztop);
    strain(:,2*u-1) = midstrain(1:3) + zbottom*midstrain(4:end);
    strain(:,2*u) = midstrain(1:3) + ztop*midstrain(4:end);
    stress(:,2*u-1) = Q_bar(:,:,u)*(strain(:,2*u-1)-DTbottom*axy(:,u));
    stress(:,2*u) = Q_bar(:,:,u)*(strain(:,2*u)-DTtop*axy(:,u));
end

% Set I & J to Arbitrary Value
I = 10;
J = 10;
% Set X Force to Arbitrary Value
Ny0 = 0;
Nx0 = 10;
Nxy0 = 0;
% Calculate Summations
for i = 1:I
    for j = 1:J
        for m = 1:I
            for n = 1:J
                L = (m-1)*J+n;
                k = (i-1)*J+j;
                if k==L
                    delta=1;
                else
                    delta=0;
                end
                % r values
                r_jn = 0;
                if rem((j-n),2)~=0
                    r_jn = 2*j/(j^2-n^2)/pi; % odd
                end
                r_nj = 0;
                if rem((n-j),2)~=0
                    r_nj = 2*n/(n^2-j^2)/pi; % odd
                end
                r_im = 0;
                if rem((i-m),2)~=0
                    r_im = 2*i/(i^2-m^2)/pi; % odd
                end
                r_mi = 0;
                if rem((m-i),2)~=0
                    r_mi = 2*m/(m^2-i^2)/pi; % odd
                end
                % G, P, & b Matrices
                G(L,k) = 0.25*Lx*Ly*pi^4*(D11*(i/Lx)^4+2*(D12+2*D66)*(i/Lx)^2*(j/Ly)^2+D22*(j/Ly)^4)*delta-2*Lx*Ly*pi^4*D16*((i/Lx)^2*(m/Lx)*(n/Ly)*r_im*r_jn+(m/Lx)^2*(i/Lx)*(j/Ly)*r_mi*r_nj)-2*Lx*Ly*pi^4*D26*((j/Ly)^2*(m/Lx)*(n/Ly)*r_im*r_jn+(n/Ly)^2*(i/Lx)*(j/Ly)*r_mi*r_nj);
                P(L) = 4*p*Lx*Ly/(pi^2*m*n);
                if rem(m,2)==0 || rem(n,2)==0
                    P(L) = 0;
                end
                bb(L,k) = 1/4*Lx*Ly*pi^2*(Nx0*(i/Lx)^2+Ny0*(j/Ly)^2)*delta+Lx*Ly*pi^2*Nxy0;
            end
        end
    end
end

% Omega-k Matrix
wk = G\P.';
% Max Deflection
X = 0:0.002:Lx;
Y = 0:0.01:Ly;
[x,y] = meshgrid(X,Y);
w0 = zeros(size(x));
kx = w0;
ky = w0;
kz = w0;

for a=1:I
    for b=1:J
        c = (a-1)*J+b;
        wij(a,b)=wk(c); % for orthotropic plate
        w0 = w0+wij(a,b)*sin(a*pi*x/Lx).*sin(b*pi*y/Ly);
        kx = kx+(a*pi/Lx)^2*wij(a,b)*sin(a*pi*x/Lx).*sin(b*pi*y/Ly);
        ky = ky+(b*pi/Ly)^2*wij(a,b)*sin(a*pi*x/Lx).*sin(b*pi*y/Ly);
        kz = kz+(a*pi/Lx)*(b*pi/Ly)*2*wij(a,b)*cos(a*pi*x/Lx).*cos(b*pi*y/Ly);
    end
end

% Required Force and Moment for 0 Curvature
R = ABD*[htepsilon0(1:3); 0; 0; 0] - otherload

% Plot Stresses & Strains
figure(1)
plot(strain,zk)
xlabel('\epsilon')
ylabel('z [m]')
legend('\epsilon_x', '\epsilon_y', '\gamma_{xy}', 'Location', 'Best')
grid on
figure(2)
plot(stress,zk)
xlabel('\sigma [Pa]')
ylabel('z [m]')
legend('\sigma_x', '\sigma_y', '\tau_{xy}', 'Location', 'Best')
grid on

% Plot Deflection
figure(3)
surf(x,y,w0)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Deflection [m]')