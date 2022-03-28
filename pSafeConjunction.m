%% Initial conditions with Random Keplerian elements
clear all;
format long;

%% Initial conditions for KOSMOS and IRIDIUM
clear all;
format long;

global m1 m2 e1 e2
% physical parameters
Re = 6378.135; % Earth radius in [m]
MU = 3.986*10^14; % [m^3/s^2]
% keep-out sphere parameters
sigma_a = 1500; % along-track standard deviation [m]
sigma_c = 500; % cross-track standard deviation [m]
n = 8; % represents 8-std

% KOSMOS 2251
e1 = 0.0016015; % eccentricity
i1 = 74.0357*pi/180; % inclination
w1 = 95.9865*pi/180; % argument of periapsis
o1 = 17.1729*pi/180; % RAAN
M1 = 264.3113*pi/180; % initial mean anomaly
E1 = -0.177528682723533; % initial eccentric anomaly
m1 = 14.31135598*2*pi/(24*3600); % mean angular velocity
a1 = nthroot((MU/m1^2),3); % altitude 787.9368km
R1 = n*sigma_c+17; % keep-out sphere radius
sigma1_Xo = (sigma_a/(a1*(1-e1*cos(E1))));
% IRIDIUM 33
e2 = 0.0002253; % eccentricity
i2 = 86.3989*pi/180; % inclination
w2 = 89.6115*pi/180; % argument of periapsis
o2 = 121.2960*pi/180; % RAAN
M2 = 270.5342*pi/180; % initial mean anomaly
E2 = 0.289902214841812; % initial eccentric anomaly
m2 = 14.34220263*2*pi/(24*3600); % mean angular velocity
a2 = nthroot((MU/m2^2),3); % altitude 777.6582km
R2 = n*sigma_c+25; % keep-out sphere radius
sigma2_Xo = (sigma_a/(a2*(1-e2*cos(E2)))); % sigma_X2

%% Compute angle of orbit intersection
% Vectors poiting towards orbital angular momentum
n1 = [sin(i1)*sin(o1);sin(i1)*cos(o1);cos(i1)];
n2 = [sin(i2)*sin(o2);sin(i2)*cos(o2);cos(i2)];
% Vector poiting towards apogee
h1 = [cos(o1)*cos(w1)-sin(o1)*sin(w1)*cos(i1); -cos(w1)*sin(o1)-sin(w1)*cos(o1)*cos(i1); sin(w1)*sin(i1)];
h2 = [cos(o2)*cos(w2)-sin(o2)*sin(w2)*cos(i2); -cos(w2)*sin(o2)-sin(w2)*cos(o2)*cos(i2); sin(w2)*sin(i2)];
% Vector pointing towards cross/intersection of both orbits
v = cross(n1,n2);
v = v/norm(v);
% Angles between apogee and cross/intersection of orbits
theta1_c = -acos(dot(h1,-v));
theta2_c = acos(dot(h2,-v));
E1_c = 2*atan2(sqrt(1-e1)*sin(theta1_c/2), sqrt(1+e1)*cos(theta1_c/2));
E2_c = 2*atan2(sqrt(1-e2)*sin(theta2_c/2), sqrt(1+e2)*cos(theta2_c/2));

E0 = [E1-E1_c; E2-E2_c]; % change of coordinates to enable small angle approx

miss_distance = abs(a2*(1-e2^2)/(1+e2*cos(theta2_c)) - a1*(1-e1^2)/(1+e1*cos(theta1_c)))

%% Yalmip initialization
yalmip('clear');

% yalmip variables
x = sdpvar(2,1);
p = sdpvar;

%% Initial set A
syms a b;
E = [a;b];
%fimplicit((E-E0)'*Ma*(E-E0)-1,'LineWidth',5,'Color','blue');

% Matrix of set A
Ma = diag([1/(n*sigma1_Xo)^2,1/(n*sigma2_Xo)^2]);
% Yalmip set
A = 1-(x-E0)'*Ma*(x-E0);

%% Unsafe set U
gamma = 0.8; % Parameter to enlarge the set found in this section

% Quadratic approximation of true unsafe set U by means of small angle approximation
c1 = cos(theta1_c);
c2 = cos(theta2_c);
s1 = sin(theta1_c);
s2 = sin(theta2_c);

Rot1 = Rotz(-o1)*Rotx(-i1)*Rotz(-w1);
Rot2 = Rotz(-o2)*Rotx(-i2)*Rotz(-w2);

a11 = Rot1(1,1);
a12 = Rot1(1,2);
a21 = Rot1(2,1);
a22 = Rot1(2,2);
a31 = Rot1(3,1);
a32 = Rot1(3,2);

b11 = Rot2(1,1);
b12 = Rot2(1,2);
b21 = Rot2(2,1);
b22 = Rot2(2,2);
b31 = Rot2(3,1);
b32 = Rot2(3,2);

% Matrix A
A11 = (a1^2*(1-e1^2)^2*((a12*c1-a11*s1)^2 + (a22*c1-a21*s1)^2 + (a32*c1-a31*s2)^2)*(1+e2*c2)^2 + ...
      a2^2*(1-e2^2)^2*((b11*c2+b12*s2)^2 + (b21*c2+b22*s2)^2 + (b31*c2+b32*s2)^2)*e1^2*s1^2 - ...
      2*a1*a2*(1-e1^2)*(1-e2^2)*((a12*c1-a11*s1)*(b11*c2+b12*s2) + (a22*c1-a21*s1)*(b21*c2+b22*s2) + (a32*c1-a31*s1)*(b31*c2+b32*s2))*e1*s1*(1+e2*c2) - ...
      (R1+R2)^2*(e1^2*s1^2*(1+e2*c2)^2));

A22 = (a2^2*(1-e2^2)^2*((b12*c2-b11*s2)^2 + (b22*c2-b21*s2)^2 + (b32*c2-b31*s2)^2)*(1+e1*c1)^2 + ...
      a1^2*(1-e1^2)^2*((a11*c1+a12*s2)^2 + (a21*c1+a22*s1)^2 + (a31*c1+a32*s1)^2)*e2^2*s2^2 - ...
      2*a1*a2*(1-e1^2)*(1-e2^2)*((b12*c2-b11*s2)*(a11*c1+a12*s1) + (b22*c2-b21*s2)*(a21*c1+a22*s1) + (b32*c2-b31*s2)*(a31*c1+a32*s1))*e2*s2*(1+e1*c1) - ...
      (R1+R2)^2*(e2^2*s2^2*(1+e1*c1)^2));

A12 = (a1^2*(1-e1^2)^2*4*((a12*c1-a11*s1)*(a11*c1+a12*s1) + (a22*c1-a21*s1)*(a21*c1+a22*s1) + (a32*c1-a31*s1)*(a31*c1+a32*s1))*e2*s2*(1+e2*c2) + ...
      a2^2*(1-e2^2)^2*4*((b12*c2-b11*s2)*(b11*c2+b12*s2) + (b22*c2-b21*s2)*(b21*c2+b22*s2) + (b32*c2-b31*s2)*(b31*c2+b32*s2))*e1*s1*(1+e1*c1) - ...
      2*a1*a2*(1-e1^2)*(1-e2^2)*(((b12*c2-b11*s2)*(a11*c1+a12*s1) + (b22*c2-b21*s2)*(a21*c1+a22*s1) + (b32*c2-b31*s2)*(a31*c1+a32*s1))*e1*s1*(1+e2*c2) + ...
                                 ((a12*c1-a11*s1)*(b11*c2+b12*s2) + (a22*c1-a21*s1)*(b21*c2+b22*s2) + (a32*c1-a31*s1)*(b31*c2+b32*s2))*e2*s2*(1+e1*c1) + ...
                                 ((a11*c1+a12*s1)*(b11*c2+b12*s2) + (a21*c1+a22*s1)*(b21*c2+b22*s2) + (a31*c1+a32*s1)*(b31*c2+b32*s2))*e1*e2*s1*s2 + ...
                                 ((a12*c1-a11*s1)*(b12*c2-b11*s2) + (a22*c1-a21*s1)*(b22*c2-b21*s2) + (a32*c1-a31*s1)*(b32*c2-b31*s2))*(1+e1*c1)*(1+e2*c2)) - ...
      (R1+R2)^2*(4*e1*e2*s1*s2*(1+e1*c1)*(1+e2*c2)));

% Vector B
B1 = (a1^2*(1-e1^2)^2*2*((a12*c1-a11*s1)*(a11*c1+a12*s1) + (a22*c1-a21*s1)*(a21*c1+a22*s1) + (a32*c1-a31*s1)*(a31*c1+a32*s1))*(1+e2*c2)^2 + ...
     a2^2*(1-e2^2)^2*2*((b11*c2+b12*s2)^2 + (b21*c2+b22*s2)^2 + (b31*c2+b32*s2)^2)*e1*s1*(1+e1*c1) - ...
     2*a1*a2*(1-e1^2)*(1-e2^2)*(((a12*c1-a11*s1)*(b11*c2+b12*s2) + (a22*c1-a21*s1)*(b21*c2+b22*s2) + (a32*c1-a31*s1)*(b31*c2+b32*s2))*(1+e1*c1)*(1+e2*c2) + ...
                                ((a11*c1+a12*s1)*(b11*c2+b12*s2) + (a21*c1+a22*s1)*(b21*c2+b22*s2) + (a31*c1+a32*s1)*(b31*c2+b32*s2))*e1*s1*(1+e2*c2)) - ...
     (R1+R2)^2*(2*e1*s1*(1+e1*c1)*(1+e2*c2)^2));

B2 = (a1^2*(1-e1^2)^2*2*((a11*c1+a12*c1)^2 + (a21*c1+a22*s1)^2 + (a31*c1+a32*s1)^2)*e2*s2*(1+e2*c2) + ...
     a2^2*(1-e2^2)^2*2*((b12*c2-b11*s2)*(b11*c2+b12*s2) + (b22*c2-b21*s2)*(b21*c2+b22*s2) + (b32*c2-b31*s2)*(b31*c2+b32*s2))*(1+e1*c1)^2 - ...
     2*a1*a2*(1-e1^2)*(1-e2^2)*(((b12*c2-b11*s2)*(a11*c1+a12*s1) + (b22*c2-b21*s2)*(a21*c1+a22*s1) + (b32*c2-b31*s2)*(a31*c1+a32*s1))*(1+e1*c1)*(1+e2*c2) + ...
                                ((a11*c1+a12*s1)*(b11*c2+b12*s2) + (a21*c1+a22*s1)*(b21*c2+b22*s2) + (a31*c1+a32*s1)*(b31*c2+b32*s2))*e2*s2*(1+e1*c1)) - ...
     (R1+R2)^2*(2*e2*s2*(1+e2*c2)*(1+e1*c1)^2));

% Constant c
c = a1^2*(1-e1^2)^2*((a11*c1+a12*s1)^2 + (a21*c1+a22*s1)^2 + (a31*c1+a32*s1)^2)*(1+e2*c2)^2 + ...
    a2^2*(1-e2^2)^2*((b11*c2+b12*s2)^2 + (b21*c2+b22*s2)^2 + (b31*c2+b32*s2)^2)*(1+e1*c1)^2 - ...
    2*a1*a2*(1-e1^2)*(1-e2^2)*((a11*c1+a12*s1)*(b11*c2+b12*s2) + (a21*c1+a22*s1)*(b21*c2+b22*s2) + (a31*c1+a32*s1)*(b31*c2+b32*s2))*(1+e1*c1)*(1+e2*c2) - ...
    (R1+R2)^2*(1+e1*c1)^2*(1+e2*c2)^2;

% Aq, Bq and c terms of ellipse in quadratic form:
% [dtheta1,dtheta2]*Aq*[dtheta1,dtheta2]' + Bq*[dtheta1,dtheta2]' + c = 0
Aq = [A11,A12;A12,A22];
Bq = [B1;B2];

% M matrix of ellipse in standard form:
%[theta1,theta2]*M*[theta1,theta2]' = 1
M = Aq/(Bq'*inv(Aq)*Bq/4-c);

% Transformation matrix to express ellipse in terms of E
diag1 = (1+tan(E1_c/2)*tan(theta1_c/2)*sqrt((1-e1)/(1+e1)))/(sqrt((1-e1)/(1+e1))+tan(E1_c/2)*tan(theta1_c/2));
diag2 = (1+tan(E2_c/2)*tan(theta2_c/2)*sqrt((1-e2)/(1+e2)))/(sqrt((1-e2)/(1+e2))+tan(E2_c/2)*tan(theta2_c/2));
T = diag([diag1,diag2]);

% Apply transformation and multiply by gamma to encapsulate original set
Mu = gamma*(T'*M*T);
% Yalmip set
U = 1-x'*Mu*x;

%% State-Space set S
Uev = min(eig(Mu));
Aev = min(eig(Ma));

% Radius of set S
alpha = sqrt(2);
Rs = alpha*(norm(E0) + 1/min(Uev,Aev));
% Yalmip set S
S = Rs^2 - (x)'*(x);

%% Bernstein Approximation
d_b = 5; % Bernstein degree polynomial
B1 = 0; % Bernstein polynomial B1 initialization
B2 = 0; % Bernstein polynomial B2 initialization

P1 = -Rs; % initial point on domain D
P2 = Rs; % final point on domain D

% baricentric coordinates
lam11 = (x(1)-P1)/(P2-P1);
lam12 = 1-(x(1)-P1)/(P2-P1);
lam21 = (x(2)-P1)/(P2-P1);
lam22 = 1-(x(2)-P1)/(P2-P1);
for i = 0:d_b
    k1 = P1 + (i/d_b)*(P2-P1);
    k2 = P1 + (i/d_b)*(P2-P1);
    bk1 = m1/(1 - e1*(cos(E1_c)*cos(k1)-sin(E1_c)*sin(k1)));
    bk2 = m2/(1 - e2*(cos(E2_c)*cos(k2)-sin(E2_c)*sin(k2)));
    B1 = B1 + bk1*nchoosek(d_b,i)*(lam11^i)*(lam12^(d_b-i));
    B2 = B2 + bk2*nchoosek(d_b,i)*(lam21^i)*(lam22^(d_b-i));
end

%% p-safety
% Barrier polynomial degree
d_h = 11
% SOS polynomials
[s1,coef1] = polynomial(x,d_h);
[s2,coef2] = polynomial(x,d_h);
[s3,coef3] = polynomial(x,d_h);
[s4,coef4] = polynomial(x,d_h);
[s5,coef5] = polynomial(x,d_h);

[h,coefh] = polynomial(x,d_h);

% Drift term via Bernstein approximation of f1 and f2
F = [B1;B2];
% Nonlinear dynamics
% f1 = m1./(1-e1.*(cos(E1_c)*cos(x(1))-sin(E1_c)*sin(x(1))));
% f2 = m2./(1-e2.*(cos(E2_c)*cos(x(2))-sin(E2_c)*sin(x(2))));

% Difussion term
G = diag([0.0001,0.0001]);

% Extended generator
Lh = jacobian(h,x)*F + 0.5*trace(G*G'*hessian(h,x));

% Constraints
constraints=[sos(s1);sos(s2);sos(s3);sos(s4);sos(s5)];
constraints=[constraints; sos(h-s1*S); sos(-Lh-s2*S); sos(p-h-s4*A); sos(h-1-s5*U)];
constraints=[constraints; p>=0];
constraints=[constraints; 1-p>=0];


%% SOLUTION
params = [coef1;coef2;coef3;coef4;coef5;coefh];
opts = sdpsettings('solver','sdpt3','sdpt3.maxit',2000);
[sol,u,Qmat,res] = solvesos(constraints,p,opts,params);

%% CHECK
%For checking which constraints are satisfied
checkset(constraints)
Q_all_eig = cellfun(@eig,Qmat,'UniformOutput',0);
Q_min_eig = cellfun(@min,Q_all_eig,'UniformOutput',1);
min(Q_min_eig)

%It shows what is p
realparam_p = double(p)

%% PREPARE ELEMENTS FOR PLOTING
%For plotting Stochastic barrier function
f = sdisplay(clean(replace(h,coefh,double(coefh)),1e-8));
f = f{1};
f = replace(f,'*','.*');
f = replace(f,'^','.^');
f = replace(f,'internal(1)','x1');
f = replace(f,'internal(2)','x2');
f = replace(f,'x(1)','x1');
f = replace(f,'x(2)','x2');

h_p = str2sym(f); % Super-martingale barrier function in symbolic

% Bernstein polynomials in symbolic
syms x1 x2;
X = [x1,x2];
B1 = 0;
B2 = 0;
lam11 = (X(1)-P1)/(P2-P1);
lam12 = 1-(X(1)-P1)/(P2-P1);
lam21 = (X(2)-P1)/(P2-P1);
lam22 = 1-(X(2)-P1)/(P2-P1);
for i = 0:d_b
    k1 = P1 + (i/d_b)*(P2-P1);
    k2 = P1 + (i/d_b)*(P2-P1);
    bk1 = m1/(1 - e1*(cos(E1_c)*cos(k1)-sin(E1_c)*sin(k1)));
    bk2 = m2/(1 - e2*(cos(E2_c)*cos(k2)-sin(E2_c)*sin(k2)));
    B1 = B1 + bk1*nchoosek(d_b,i)*(lam11^i)*(lam12^(d_b-i));
    B2 = B2 + bk2*nchoosek(d_b,i)*(lam21^i)*(lam22^(d_b-i));
end

% Drift in symbolic
Fsym = [B1,B2]';

% Construct Extended generator in symbolic
Lh_p=jacobian(h_p,[x1,x2])*Fsym + 0.5*trace(G*G'*hessian(h_p,[x1,x2]));
Lhfunc = symfun(Lh_p,[x1,x2]);

f = eval(['@(x1,x2)' f]);

%% PLOTS
figure(1);
hold on;
% SET A
Ai=fimplicit((E-E0)'*Ma*(E-E0)-1,[E0(1)-0.01 E0(1)+0.01 E0(2)-0.01 E0(2)+0.01],'LineWidth',5,'Color','blue');
xa = Ai.XData;
ya = Ai.YData;
% SET S
Si=fimplicit((E)'*(E) - 2*(norm(E0) + sqrt(1/min(Uev,Aev)))^2, [-0.25 0.25 -0.25 0.25], 'LineWidth',5,'Color','g');
xs = Si.XData;
ys = Si.YData;
% SET U
Ui = fimplicit(E'*Mu*E-1, [-0.01 0.01 -0.01 0.01],'LineWidth',5,'Color','red');
xu = Ui.XData;
yu = Ui.YData;

% SET U PROJECTED ONTO h(X)
zu = zeros(1,length(xu));
for i=1:length(xu)
    zu(i) = double(f(xu(i),yu(i)));
    if mod(i,2)==0
        plot3([xu(i),xu(i)],[yu(i),yu(i)],[0,zu(i)],'Color','r');
    end
end
plot3(xu,yu,zu,'LineWidth',5,'Color','r');

% BARRIER FUNCTION h
%Note: The Barrier function is positive on the semialgebraic set S, bigger
%than 1 in the semialgebraic set U and smaller than p in the semialgebraic
%set A
[X1,X2] = meshgrid(-0.22:0.01:0.22,-0.22:0.01:0.22);
Zh=double(f(X1,X2)); % evaluate barrier h
hf = mesh(X1,X2,Zh,'FaceColor','None','EdgeColor','k');
zlim([0,1]);

% EXTENDED GENERATOR Lh
figure(2)
hold on;
ZLh=double(Lhfunc(X1,X2)); % evaluate Lh
hLf = mesh(X1,X2,ZLh,'FaceColor','None','EdgeColor','k');
fimplicit((E)'*(E) - 2*(norm(E0) + sqrt(1/min(Uev,Aev)))^2, [-0.25 0.25 -0.25 0.25], 'LineWidth',5,'Color','g');
zlim([-0.1,0]);
view(45,45);
grid on
xlabel('X_1')
ylabel('X_2')
zlabel('Lh(X)')
title('Lh(X) function')

% CONTOUR PLOTS OF LEVEL SETS 1 and p
figure(1);
[X1c,X2c] = meshgrid(-0.22:0.001:0.22,-0.22:0.001:0.22);
Zc=double(f(X1c,X2c));

colormap jet;
contcol_b = [0,0,0.6];
contcol_r = [0.6,0,0];
Cp = contour3(X1c,X2c,Zc,[realparam_p,realparam_p],'ShowText','off','LineWidth',3,'Color','m');
C1 = contour3(X1,X2,Zh,[1,1],'ShowText','off','LineWidth',3,'Color','k');hold on;
C05 = contour3(X1,X2,Zh,[0.5,0.5],'ShowText','off','LineWidth',3,'Color','k');
C06 = contour3(X1,X2,Zh,[0.6,0.6],'ShowText','off','LineWidth',3,'Color','k');
C07 = contour3(X1,X2,Zh,[0.7,0.7],'ShowText','off','LineWidth',3,'Color','k');
C02 = contour3(X1,X2,Zh,[0.2,0.2],'ShowText','off','LineWidth',3,'Color','k');
C01 = contour3(X1,X2,Zh,[0.1,0.1],'ShowText','off','LineWidth',3,'Color','k');
C03 = contour3(X1,X2,Zh,[0.3,0.3],'ShowText','off','LineWidth',3,'Color','k');
C04 = contour3(X1,X2,Zh,[0.4,0.4],'ShowText','off','LineWidth',3,'Color','k');

contours = contour3(X1,X2,Zh,[0.2,0.4,0.6,0.8],'ShowText','off','LineWidth',3,'Color','k');

view(45,45);
clabel(Cp,'FontSize',15)
grid on
xlabel('X_1')
ylabel('X_2')
zlabel('h(X)')
title('h(X) function')

%% SAVE DATA
x_len = length(-0.22:0.01:0.22);
y_len = length(-0.22:0.01:0.22);
h_data = []; %zeros(floor(data_length*data_length/4),3);
Lh_data = [];
for i=1:x_len
    for j=1:y_len
        h_data = [h_data;[X1(j,i),X2(j,i),Zh(j,i)]];
        Lh_data = [Lh_data;[X1(j,i),X2(j,i),ZLh(j,i)]];
    end
end

%% Functions
% Rotz
function Rz = Rotz(t)
    Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
end
% Rotz
function Rx = Rotx(t)
    Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
end
%Newton-Raphson
function E = NewRap(M, e)
    E = M;
    TOL = 0.0001;
    while abs(E - e*sin(E) - M) > TOL
        E = E - 0.01*(E-e*sin(E)-M)/(1-e*cos(E));
    end
end
