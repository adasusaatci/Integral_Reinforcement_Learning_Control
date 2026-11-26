%% LONGITUDINAL MODEL OF A-7A AIRCRAFT
%% Temizlik
clear all;
clc;
close all;
%%
global A B K x_n u_n

A = [-1.0190 1; 0.8223 -1.0774];
B = [0; -0.1756]
K = [16.8696 -10.5911]; 
x_initial=[0; 10];
% K = [1.0222 2.2693 -38.9085 -14.7734]

[x_n,u_n]=size(B);
x_save=[];
t_save=[];

% Initiliazing the cost function parameters Q and R
Q=eye(x_n);
R=eye(u_n);

P_old=zeros(x_n);

%Data is obtained from the system every 0.05 seconds
T=0.05;  
Dxx=[];XX=[];XU=[];
X=[x_initial;kron(x_initial',x_initial')';kron(x_initial,zeros(u_n,1))]';

P=eye(x_n)*10; 
it=0;           
p_save=[];      
k_save=[];      

[P0,K0]=care(A,B,Q); % Ideal solution 
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
    scatter(it,P_old(2,2),'+','g','Linewidth',2)
    hold on
  
 
    p_save=[p_save,norm(P-P0)];    
    
    K=inv(R)*transpose(B)*P;          
    k_save=[k_save,norm(K-K0)];     
  
    

    it=it+1     
end
fprintf('\nNumber of iterations to obtain optimal control solution: %d\n',it)



figure(4)

scatter(it,P0(1,1),'*','b')
hold on
scatter(it,P0(2,2),'*','g')
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
axis([0,20,-0.2,6])   
legend('x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8','x_9','x_1_0','x_1_1','x_1_2')
xlabel('Time (s)')
title('System States')

figure(5)

plot(t_final,x_final(:,7:8),'Linewidth',2)
axis([0,20,-0.2,6])   
legend('x_7','x_8','x_9','x_1_0','x_1_1','x_1_2')
xlabel('Time (s)')
title('System States (Close to Optimal Solution)')


   function dX=mysys(t,X)
       global A B K x_n
        x=X(1:x_n);

        u=-K*x;
        dx=A*x+B*u;
        dxx=kron(x',x')';
        dux=kron(x',u')';
        dX=[dx;dxx;dux];
	end