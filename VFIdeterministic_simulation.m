close all
clear all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
T_mat = [0.977, 1-0.977; 1 - 0.926, 0.926];
A = [1.0082, 0.9736];

A_simulation_1=zeros(1,1000);

A_simulation_1(1,1)=1;

for i=2:1000
    ran=unifrnd(0,1);
    if A_simulation_1(1,i-1)==1
        if ran<=T_mat(1,1)
            A_simulation_1(1,i)=1;
        else
            A_simulation_1(1,i)=0;
        end
    else
        if ran<=T_mat(2,1)
            A_simulation_1(1,i)=1;
        else
            A_simulation_1(1,i)=0;
        end
    end
end

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons(:,:,1) = A(1) * k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
cons(:,:,2) = A(2) * k_mat .^ alpha + (1 - delta) * k_mat - k_mat';

ret = cons .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret(cons < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(2, num_k);
i=1
while dis > tol
    % an alternative, more direct way to compute the value array:
    value_mat_alt = ret + beta * ...
        repmat(permute((T_mat * v_guess), [3 2 1]), [num_k 1 1]);
    
    % compute the utility value for all possible combinations of k and k':
    value_mat(:,:,1) = ret(:,:,1) + beta * ( ...
        T_mat(1,1) * repmat(v_guess(1,:), [num_k 1]) + ...
        T_mat(1,2) * repmat(v_guess(2,:), [num_k 1]));
    
    value_mat(:,:,2) = ret(:,:,2) + beta * ( ...
        T_mat(2,1) * repmat(v_guess(1,:), [num_k 1]) + ...
        T_mat(2,2) * repmat(v_guess(2,:), [num_k 1]));
    
    
    
    % find the optimal k' for every k:
    [vfn, pol_indx] = max(value_mat, [], 2);
    vfn =permute(vfn, [3,1,2]);
    
    
    % what is the distance between current guess and value function
    dis = max(max(abs(vfn - v_guess)));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
    
    i=i+1
end

g = k(pol_indx); % policy function
g = permute(g, [3, 1, 2]);

plot(k,vfn);
xlabel('k')
ylabel('vfn')
legend('Ah','Al')
figure
plot(k,g);
xlabel('k_t')
ylabel('k_t+1')
legend('Ah','Al')
figure
plot(k,g-(1-delta)*repmat(k, [2,1]));
xlabel('k')
ylabel('saving')
legend('Ah','Al')

A_simulation=zeros(1,1000);
A_simulation(A_simulation_1==1)=A(1);
A_simulation(A_simulation_1==0)=A(2);



k_simulation=zeros(1,1000);
k_simulation(1,1)=35;
ik_simulation=zeros(1,1);
for i=2:1000
    if A_simulation(1,i-1)==A(1);
        ik_simulation=find(k==k_simulation(1,i-1));
        k_simulation(1,i)=g(1,ik_simulation);
    else
        ik_simulation=find(k==k_simulation(1,i-1));
        k_simulation(1,i)=g(2,ik_simulation);
    end
end

y=A_simulation.*(k_simulation.^alpha);

std(y)/mean(y)



        
        


