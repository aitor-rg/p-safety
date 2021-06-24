%%Initilization of the algorithm
clear 
x = sdpvar(2,1); %the states
%sdpvar p; %the probability to be found
p = 0.1;
D = 6; % the degree of all involved SOS polynomials
%%We use 8 unknown (sum of squares) polynomials
[s1,coefs1] = polynomial(x,D);
[s2,coefs2] = polynomial(x,D);
[s3,coefs3] = polynomial(x,D);
[s4,coefs4] = polynomial(x,D);
[s5,coefs5] = polynomial(x,D);
[s6,coefs6] = polynomial(x,D);
[s7,coefs7] = polynomial(x,D);
[s8,coefs8] = polynomial(x,D);

%%The switched diffusion system is considered with two subsystems and
%exponencial switched defined by  two transition rates
%The subsystems:
%   dx = f_1 dt + sigma_1 dW
%   dx = f_2 dt + sigma_2 dW
%The transition rates lambda_{12} and lambda_{21}.


sigma1 = .5*eye(2); % the diffusion matrix of System #1
sigma2 = .5*eye(2); % the diffusion matrix of System #2

f1 = [1; 1.5]; %drift of System #1
f2 = [1.5; 1]; %drift of System #1

lambda12 = 10;
lambda21 = 10;

%The difinition of three involved sets
%The state space  Y = [Y(i,x) \geq 0, i = 1,2] 
Y1 = 10^2 - x(1)^2 - x(2)^2;
Y2 = 10^2 - x(1)^2 - x(2)^2;
%The initial set A = [A(i,x) \geq 0, i = 1,2]
A1 = 1 - x(1)^2 - x(2)^2;
A2 = 1 - x(1)^2 - x(2)^2;
%The forbidden state  U = [U(i,x) \geq 0, i = 1,2]
U1 = 1 - (x(1)-5)^2 - (x(2)-5)^2;
U2 = 1 - (x(1)-5)^2 - (x(2)-5)^2;


%The two test functions for Subsystem#1 and Subsystem#2 h = (h_1,h_2)
[h1,coefh1] = polynomial(x,D);
[h2,coefh2] = polynomial(x,D);

%The infinitesmall generator consists of two components
Lh1 = 0.5 * trace(sigma1*sigma1'*hessian(h1,x)) + jacobian(h1,x)*f1 + lambda12*(h2-h1);
Lh2 = 0.5 * trace(sigma2*sigma2'*hessian(h2,x)) + jacobian(h2,x)*f2 + lambda21*(h1-h2);
%Pure diffusion processes for the test
%Lh1 = 0.5 * trace(sigma1*sigma1'*hessian(h1,x)) + jacobian(h1,x)*f1;
%Lh2 = 0.5 * trace(sigma2*sigma2'*hessian(h2,x)) + jacobian(h2,x)*f2;
%Pure dynamical (deterministic) systems for the test
%Lh1 = jacobian(h1,x)*f1;
%Lh2 = jacobian(h2,x)*f2;



    %Compute trace by hand
    %DifPart1 = sigma1*sigma1'*hessian(h1,x);
    %DifPart2 = sigma2*sigma2'*hessian(h2,x);
    %Tr1 = DifPart1(1,1) + DifPart1(2,2);
    %Tr2 = DifPart2(1,1) + DifPart2(2,2);
    %Lh1 = 0.5 * Tr1 + jacobian(h1,x)*f1 + lambda12*(h2-h1);
    %Lh2 = 0.5 * Tr2 + jacobian(h2,x)*f2 + lambda21*(h1-h2);

%Constraints of the optimisation
%We use eight sos: s1, s2, s3, s4, s5, s6, s7, s8
constr = [sos(s1);sos(s2);sos(s3);sos(s4);sos(s5);sos(s6);sos(s7);sos(s8)];
%h >= 0 on the state space Y 
constr = [constr; sos(h1 - s1*Y1)];
constr = [constr; sos(h2 - s2*Y2)];

%The infnitesmall generator -Lh >= 0 on the state space Y 
constr = [constr; sos(-Lh1 - s3*Y1)];
constr = [constr; sos(-Lh2 - s4*Y2)];

% p - h >= 0 on A
constr = [constr; sos(p-h1 - s5*A1)];
constr = [constr; sos(p-h2 - s6*A2)];

% h - 1 >= 0 on U
constr = [constr; sos(h1-1 - s7*U1)];
constr = [constr; sos(h2-1 - s8*U2)];


%SOLUTION
params = [coefs1;coefs2;coefs3;coefs4;coefs5;coefs6;coefs7;coefs8;coefh1;coefh2];
%options = sdpsettings('solver','mosek');
options = sdpsettings('solver','sdpt3');
solvesos(constr,p,options,params);
%[sol,u,Qmat,res] = solvesos(constr,p,options,params);

%For checking which constraints are satisfied
%checkset(constr)

%Q_all_eig = cellfun(@eig,Qmat,'UniformOutput',0);

%Q_min_eig = cellfun(@min,Q_all_eig,'UniformOutput',1);
%min(Q_min_eig)

%It shows what is p, h_1, and h_2
realparam_p = double(p)
%realparam_h1 = double(coefh1);
%realparam_h2 = double(coefh2);
%x = sdpvar(1,1); y = sdpvar(1,1);
%vv = monolist([x; y],D);
%poly_h1 = vectorize(sdisplay(realparam_h1'*vv))
%poly_h2 = vectorize(sdisplay(realparam_h2'*vv))         
 
 
 
