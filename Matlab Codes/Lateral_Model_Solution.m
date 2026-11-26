%% LATERAL MODEL OF AN AIRCRAFT
%% Temizlik
clear all;
clc;
close all;
%%
global A B K x_n u_n
A = [-0.3220 0.0640 0.0364 -0.9917; 0 0 1 0.0393; -30.6490 0 -3.6784 0.6646; 8.5395 0 -0.0254 -0.4764];
A_nom = [-0.3220 0.0640 0.0364 -0.9917; 0 0 1 0.0393; -22.7599 -6.9980 -5.1138 -0.3242; 8.8560 -0.3183 -0.1161 -2.3489];
B = [0 0; 0 0; -0.7331 0.1315; -0.0319 -0.0620];
K = [10.6901 -9.5824 -2.0328 -6.1944; -0.3982 -0.2043 -0.4170 -27.0142];
%K = [0 0 0 0;0 0 0 0]
% x = [pitch angle, ..., angle of attack, pitch rate]^T
C = [0 0 0 57.2958; 0 0 57.2958 0; 57.2958 0 0 0; 0 57.2958 0 0]
D = 0;
x_initial=[0; 5; 0; 0];
[x_n,u_n]=size(B);

x_save=[];
t_save=[];

% Initiliazing the cost function parameters Q and R
Q=eye(x_n);
R=eye(u_n);

P_old=zeros(x_n);

%Data is obtained from the system every 0.05 seconds
T=0.05;  
Dxx=[];
XX=[];
XU=[];

X=[x_initial;kron(x_initial',x_initial')';kron(x_initial,zeros(u_n,1))]';

P=eye(x_n)*10; 
it=0;           
p_save=[];      
k_save=[];      

[P0,L0,K0]=care(A_nom,B,Q); % Ideal solution 
k_save=[norm(K-K0)];


while norm(P-P_old)>10^-3   
    for i=1:20
    % Reinforcement on the system
    [t,X]=ode45(@mysys, [it+(i-1)*T,it+i*T],X(end,:));

    
    Dxx=[Dxx;kron(X(end,1:x_n),X(end,1:x_n))-kron(X(1,1:x_n),X(1,1:x_n))];
    XX=[XX;X(end,x_n+1:x_n+x_n^2)-X(1,x_n+1:x_n+x_n^2)];
    XU=[XU;X(end,x_n+x_n^2+1:end)-X(1,x_n+x_n^2+1:end)];

  
    x_save=[x_save;X];
    t_save=[t_save;t];
    end
P_old=P;                        
    
     
    QK=Q+K'*R*K;                    
    X2=XX*kron(eye(x_n),K');         
    X1=[Dxx,-X2-XU];                
    Y=-XX*QK(:);                    
    
    pp = pinv(X1)*Y;
	
	P = reshape(pp(1:x_n*x_n), [x_n, x_n]);
	P = (P +P')/2;
      
     figure(4)
    
    scatter(it,P_old(1,1),'o','filled','b')
    hold on
    scatter(it,P_old(2,3),'+','g','Linewidth',2)
    hold on
    scatter(it,P_old(2,4),'o','r','Linewidth',2)
    hold on
    scatter(it,P_old(3,3),'s','c','Linewidth',2)
    hold on
    legend('P(1,1)','P(2,3)','P(2,4)','P(4,4)')
 
    p_save=[p_save,norm(P-P0)];    
    
    K=inv(R)*transpose(B)*P;          
    k_save=[k_save,norm(K-K0)];     
  
    

    it=it+1     
end
fprintf('\nNumber of iterations to obtain optimal control solution: %d\n',it)



figure(4)

scatter(it,P0(1,1),'*','b')
hold on
scatter(it,P0(2,3),'*','g')
hold on
scatter(it,P0(2,4),'*','r')
hold on
scatter(it,P0(4,4),'*','c')
hold on
axis([0 it -300 500])
xlabel('Time (s)')
title('Cost Matrix Parameters (P)')
legend('P(1,1)','P(2,3)','P(2,4)','P(4,4)','P(1,1)-Optimal','P(2,3)-Optimal','P(2,4)-Optimal','P(4,4)-Optimal')
    
figure(1)
plot([0:length(p_save)-1],p_save,'o',[0:length(p_save)-1],p_save)
legend('||P_i-P^*||')
xlabel('Number of iterations')

figure(2)
plot([0:length(k_save)-1],k_save,'^',[0:length(k_save)-1],k_save)
legend('||K_i-K^*||')
xlabel('Number of iterations')


% Post-learning simulation
[tt,xx]=ode23(@mysys,[t(end) 20],X(end,:)');

% Keep track of the post-learning trajectories
t_final=[t_save;tt];
x_final=[x_save;xx];


figure(3)

plot(t_final,x_final(:,1:8),'Linewidth',2)
axis([0,7,-0.2,6])   
legend('x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8')
xlabel('Time (s)')
title('System States')

figure(5)

plot(t_final,x_final(:,7:8),'Linewidth',2)
axis([0,20,-0.2,6])   
legend('x_7','x_8')
xlabel('Time (s)')
title('System States (Close to Optimal Solution)')


   function dX=mysys(t,X)
       global A B K x_n
        x=X(1:x_n);

        u=-K*x;
        dx=A*x+B*u;
        dxx=kron(x',x')';
        dux=kron(x',u')';
        dX=[dx(:,1);dxx;dux];
	end