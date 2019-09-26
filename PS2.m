%%%%%% PS2 VFI
clc
close all

%%% Set up parameters
alpha = .35;
beta = .99;
delta = .025;
sigma = 2;
a_h = 1.1;
a_l = .678;

%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k
k = linspace(k_min, k_max, num_k);
k_mat = repmat(k', [1 num_k]);

%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim(cols): k' chosen for tomorrow
pi_h = .977;
pi_l = .926;
pi = [pi_h,1-pi_h; 1-pi_l,pi_l];
cons_h = a_h*k_mat .^alpha+(1-delta)*k_mat-k_mat';
cons_l = a_l*k_mat .^alpha+(1-delta)*k_mat-k_mat';
ret_h = cons_h .^(1-sigma)/(1-sigma);
ret_l = cons_l .^(1-sigma)/(1-sigma);

% negtive consumption is not possible, make it irrelevant by assigning
% it vary large negative utility
ret_h(cons_h < 0) = -Inf;
ret_l(cons_l < 0) = -Inf;

%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(1, num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    % value_mat = ret + beta * repmat(v_guess, [num_k 1]);
    value_mat_h = ret_h + beta * repmat(pi(1,:)*v_guess, [num_k,1]);
    value_mat_l = ret_l + beta * repmat(pi(2,:)*v_guess, [num_k,1]);
    % find the optimal k' for every k:
    % [vfn, pol_indx] = max(value_mat, [], 2);
    [vfn_h, pol_indx_h] = max(value_mat_h, [], 2);
    [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
    vfn = [vfn_h';vfn_l'];
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
end
% policy function
g_h = k(pol_indx_h); 
g_l = k(pol_indx_l);
g = [g_h';g_l'];

plot(k,vfn)
figure
plot(k,g)
figure
plot(k,g-(1-delta)*repmat(k,[2,1]));

%%% Simulation
T_sim = 5000;
rng(1);
rand_num = rand(T_sim,1);
A_sim = zeros(T_sim,1);
A_sim(1) = 1;
for t = 1:T_sim
    if A_sim(t) == 1
        if rand_num(t) < pi_h
            A_sim(t+1) = 1;
        else A_sim(t+1) = 2;
        end
    elseif rand_num(t) < pi_l
        A_sim(t+1) = 2;
    else A_sim(t+1) = 1;
    end
end
k_sim_index = zeros(T_sim,1);
k_sim_index(1) = 10;
for t = 1:T_sim
    k_sim_index(t+1) = pol_index(A_sim(t), k_sim_index(t));
end

%Then do for output Y = AK^alpha
y_sim = A_sim(t)*k_sim_index^alpha
std=std(y);








