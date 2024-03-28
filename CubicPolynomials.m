clear all;
clear variables;

addpath('Data\')
% addpath('functions')

Dat=load('dataSimHeights.mat');
t=Dat.t;
l=Dat.l;

len_t=length(t);

%% Generate equidistant nodes
del_k = 0.5;
nodes = 0:del_k:max(t)+0.5;
len_nodes = length(nodes);

% Create a Generate equidistant nodes
t_dense = linspace(0, max(t), 1000)';
A = design_matrix(nodes, t);
A2 = design_matrix(nodes, t_dense);



%% C1 approximation
N = A' * A;
n = A' * l;
x_cap_C1 = N \ n;
l_cap_C1 = A * x_cap_C1;
v_cap_C1 = l_cap_C1 - l;
l_den_C1 = A2 * x_cap_C1;

%% C2 approximation
len_node_c2 = length(nodes) * 2;
C2_B = zeros(len_nodes - 2, len_node_c2);
j = 1;
k = 1;
for i = 2:len_nodes - 1
    nn = nodes(i) - nodes(i-1);
    nm = nodes(i+1) - nodes(i);
    C2_B(j, k) = 6 / nn^2;
    C2_B(j, k+1) = 2 / nn;
    C2_B(j, k+2) = -6 / nn^2 + 6 / nm^2;
    C2_B(j, k+3) = 4 / nn + 4 / nm;
    C2_B(j, k+4) = -6 / nm^2;
    C2_B(j, k+5) = 2 / nm;
    j = j + 1;
    k = k + 2;
end

A_C2 = [A; C2_B];
n_C2 = [l; zeros(size(C2_B, 1), 1)]; 
x_cap_C2 = (A_C2' * A_C2) \ (A_C2' * n_C2);
l_cap_C2 = A * x_cap_C2;
v_cap_C2 = l_cap_C2 - l;
l_den_C2 = A2 * x_cap_C2;

%% Manual constraints
twice_diff_node_index = find(nodes == 5);  % Index of κ = 5 node

AA = zeros(2, len_node_c2);
AA(1, :) = 0;
AA(1, 1) = 1;
AA(2, :) = C2_B(twice_diff_node_index, :);

bb = [0; 0];  % Constraints: f(0) = 0, f"(5) = 0

combined_mat_lsa = [A' * A AA'; AA zeros(2, 2)] \ [A' * l; bb];
x_cap_lsa = combined_mat_lsa(1:len_node_c2); 
l_cap_lsa = A * x_cap_lsa;
v_cap_lsa = l_cap_lsa - l;

l_den_lsa = A2 * x_cap_lsa;

%% Visualization

figure(1)
xline(nodes, '-', 'HandleVisibility', 'off');
hold on;
plot(t, l, 'r*');
hold on;
plot(t_dense, l_den_C1, 'g-')
xlabel('Time in years');
ylabel('Height in meters');
title('Approximation of C1');
legend('Observations', 'C1 Approximation','Location','northwest');
hold off;


figure(2)
xline(nodes, '-', 'HandleVisibility', 'off');
hold on;
plot(t, l, 'g*');
hold on;
plot(t_dense, l_den_C2, 'r-')
xlabel('Time in years');
ylabel('Height in meters');
title('Approximation of C2');
legend('Observations', 'C2 Approximation','Location','northwest');
hold off;

figure(3)
xline(nodes, '-', 'HandleVisibility', 'off');
hold on;
plot(t, l, 'r*');
hold on;
plot(t_dense, l_den_C1, 'b-');
plot(t_dense, l_den_C2, 'g-')
plot(t_dense, l_den_lsa, 'black')
hold off;
xlabel('Time in years');
ylabel('Height in meters');
title('C1, C2 and Manual Constraints');
legend('Observations', 'C1 Approximation', 'C2 Approximation','Manual constraints','Location','northwest');
hold off;