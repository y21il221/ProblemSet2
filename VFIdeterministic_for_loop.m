close all
clear all
%%%% Set up parameters
t=cputime
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
T_mat = [0.977, 1-0.977; 1 - 0.926, 0.926];
A = [1.1, 0.678];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);


value_mat_h=zeros(1,1)
value_mat_l=zeros(1,1)
value_mat_h11=zeros(1,num_k);
value_mat_l11=zeros(1,num_k);
k_h3=zeros(1,num_k)
k_l3=zeros(1,num_k)
value_mat_h1=zeros(1,num_k);
value_mat_l1=zeros(1,num_k);

dis = 1; tol = 1e-06;
p=1
while dis > tol
i=1;
K=k_max/(num_k-1)
for k_1=0:K:k_max
    j=1;
    value_mat_h3=-Inf;
    value_mat_l3=-Inf;    
    for k_2=0:K:k_max
        cons1 =  A(1)*k_1 ^ alpha + (1 - delta) * k_1 - k_2;
        cons2 =  A(2)*k_1 ^ alpha + (1 - delta) * k_1 - k_2;
        ret1 = cons1^ (1 - sigma) / (1 - sigma);
        ret2 = cons2^ (1 - sigma) / (1 - sigma);        
        ret1(cons1<0)=-Inf;
        ret2(cons2<0)=-Inf;
        value_mat_h = ret1 + beta *(T_mat(1,1)*value_mat_h1(1,j)+T_mat(1,2)*value_mat_l1(1,j));
        value_mat_l = ret2 + beta *(T_mat(2,1)*value_mat_h1(1,j)+T_mat(2,2)*value_mat_l1(1,j));
        
        
        if value_mat_h>value_mat_h3
            value_mat_h3=value_mat_h;
            k_h3(1,i)=k_2;
        else
            value_mat_h3=value_mat_h3;
        end
        
        if value_mat_l>value_mat_l3
            value_mat_l3=value_mat_l;
            k_l3(1,i)=k_2;
        else
            value_mat_l3=value_mat_l3;
        end
        
        
        j=j+1;
    end
    value_mat_h11(1,i)=value_mat_h3;
    value_mat_l11(1,i)=value_mat_l3;
        i=i+1;
end

value_mat11=[value_mat_h11;value_mat_l11];
value_mat1=[value_mat_h1;value_mat_l1];
dis=max(max(abs(value_mat11-value_mat1)));
value_mat_h1=value_mat_h11;
value_mat_l1=value_mat_l11;
p=p+1
end

value_mat_11=[value_mat_h11;value_mat_l11];
plot(k,value_mat_11);
xlabel('k')
ylabel('vfn')
legend('Ah','Al')
figure
k_3=[k_h3;k_l3];
plot(k,k_3);
xlabel('k_t')
ylabel('k_t+1')
legend('Ah','Al')
figure
plot(k,k_3-(1-delta)*repmat(k, [2,1]));
xlabel('k')
ylabel('saving')
legend('Ah','Al')

e=cputime-t




