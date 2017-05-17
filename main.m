clear all
clc
close all;
global x0 A B N G WW S xx
% s=tf('s');
% systf=(1-s)/(s^2+s+1);
% [A,B,C,D]=tf2ss([0 -1 1],[1 1 1]);
A=[.8,.1;-.2,1];
B=[0.01;.1];
 C=[1 0;0 1];
D=[0];
x0=[1;-1.5];

x(1:2,1)=x0;

%-----------------------------------------------------
% Construct model
%------------------------------------------------
%Specify system here in state space form with A B C D
Ts = 1.0; %Sample interval
% A=[0.9 0;0.9 0.9]; B=[0.1;0.2]; C=[1 1]; D=0;
%Specify system dimensions
n_in = size(B,2); %# inputs
n_out = size(C,1); %# outputs
n_states = size(A,1); %# states
%------------------------------------------------

% Specify horizons
%------------------------------------------------
M =3; %Control Horizon
N = M; %Prediction Horizon
%------------------------------------------------------
% Initial Matrix Filling
%------------------------------------------------------

%------------------------------------------------
% Quadratic cost weighting Q (error) and R (control) generation.
%------------------------------------------------
r = diag([10*ones(1,n_in)]); %Weights on input moves (i.e. ||u(t)-u(t-1)||)
q = diag([1,0.1]); %Weights on output deviation from setpoint
Q = sparse((eye(N*n_out)));
R = diag([.1*ones(1,N)]);
x1=x;
%-----------------------------------------------
% Create an observer
%-----------------------------------------------
[Lambda Phi]=largematrices(N,M,n_in,n_out,n_states,A,B,C,D);
global H f

H = full(Phi'*Q*Phi) + R;
%-----------------------------------------------
% Create plant constraints
%-----------------------------------------------
lb = -1*ones(n_in,1); ub = 1*ones(n_in,1);

Lb = kron(ones((M),1),lb); Ub = kron(ones((M),1),ub);
%------------------------------------------------
% Set initial conditions
%------------------------------------------------
sim_length = 500; %Simulation length in sec
T = sim_length/Ts; %number of samples in simulation
t = 0:Ts:sim_length-Ts; %time steps
%-----------------------------------------------
% Set system initial conditions......
%-----------------------------------------------
opts=optimset('MaxFunEvals',10000000000,'MaxIter',100000000);%,'FunValCheck','on');


u = zeros(n_in,T); %Input record
u_old = zeros(n_in,1);
y = zeros(n_out,T); %Output record

opt=optimset('display','off');
%------------------------------------------------

G=[-1 0;1 0;-1 0;0 1;0 -1];
W=ones(5,1);
S=zeros(5,2);
S(1,2)=1;
%------------------------------------------------
% Generate noise signal
%------------------------------------------------
noise = 1e-2*cumsum(randn(T,n_out)')' + 1e-2*randn(T,n_out);
[bfilt,afilt] = butter(4,0.9);
noise = filter(bfilt,afilt,noise)';
%------------------------------------------------
% MPC Evolution
%------------------------------------------------
for i = 1:T,
    %SIMULATE REAL PLANT-------------------------------
    y(:,i) = C*x(:,i) + D*u(:,i);
    %Add disturbance
    % if i>300, y(:,i)=y(:,i)+0.5*ones(n_out,1); end
    %CONTROLLER STARTS HERE----------------------------
    %Observe state
    fd = full([Phi'*Q*Lambda]);
    
    %Solve dynamic optimisation problem
    f = fd*[x(:,i)];
        f1 = fd*[x1(:,i)];

    % H=[4.2 2;2 2.2];f=[2 0;6 2]'*x(:,i);
    % U = quadprog(H,f,[],[],[],[],Lb,Ub,[],opt);
     U1 = quadprog(H,f1,[],[],[],[],Lb,Ub,[],opt);
    WW=W+S*x(:,i);
    xx=x(:,i);
       [U,j,F]=fmincon(@fun,zeros(1*N,1));
       l(i)=U'*(H)*U+f'*U;
i
%         ly(i)=(U(1));
%         lx(i)=(-f(1));
%        
%         ll(:,i)=-pinv(H)*[i;i];
%         ly1(i)=(U1(1));
%         lx1(i)=(-f1(1));
% if U'*(H)*U+f'*U>=1e-6
%     2222222222222222
% end
        u(:,i+1) = U(1:n_in); u_old = u(:,i);
        u1(:,i+1) = U1(1:n_in); u_old1 = u1(:,i);

        dd=[0;0];
    if  i<320
        dd=[0;0];
    elseif i<325
        dd=[0;+0.5];
    end
    x(:,i+1) = A*x(:,i) + B*u(:,i+1)+dd;% + 1.2*noise(:,i);
    x1(:,i+1) = A*x1(:,i) + B*u1(:,i+1)+dd;% + 1.2*noise(:,i);

end;
%------------------------------------------------
% Plots
%------------------------------------------------
figure

subplot(2,1,1);
% plot(ref','r--','Linewidth',2); hold on
plot(x(1,:)','b-','Linewidth',0.5); hold on;
plot(x(2,:)','k-','Linewidth',0.5); hold on;

ylabel('Output and reference');
plot([1 T],[-2 -2],'r--','Linewidth',0.5); hold on;
subplot(2,1,2);
% plot([1 T],[lb lb],'r--','Linewidth',2);
stairs(u(:,2:T+1)','b-','Linewidth',0.5); hold on;
xlabel('Time (sample number)');
ylabel('Input and constraints');
axis([1 T 2*lb 2*ub]);
plot([1 T],[-1 -1],'r--','Linewidth',0.5); hold on;
plot([1 T],[1 1],'r--','Linewidth',0.5); hold on;






















% opts=optimset('MaxFunEvals',1000000,'MaxIter',100000,'FunValCheck','on');
% for k=1:50
%     x0=x(:,k);
%     [U,j]=fminsearch(@fun,zeros(1*N,1),opts);
%     u=U(1:1);
%     j
%     uu(:,k)=u;
%         x(:,k+1)=A*x(:,k)+B*u(:);
% end
% y=C*x;