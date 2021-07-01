%% this script is aim to simulate the random trend and attractive trend
%%property of the base particles and solvent particles
r_d = 0.5; %min unit
c_total =0.7;c_p = 0.20;c_base = 0.20;
R_eff_vol = 2*r_d; R_eff_attract = 2*r_d;
r_e = 0.01 * r_d;  % lattice or random
p = 0.75; %  desorption
L = 50; 
r_p = r_d;r_b = r_d;
dt = 0.02; %update step
N_t = floor(L^2 / (pi * r_d^2)*c_total);
%%create the init state base and particles
N_b = floor(N_t*c_total*c_base);
N_c =  floor(N_t*c_total*c_p);
pts_base = poissonDisc([L,L],2*r_d,N_b);
pts_pt = poissonDisc([L,L],2*r_d,N_c);
%%find_neighbor

%disp(base)
scatter(pts_base(:,1),pts_base(:,2),'filled')
hold on
scatter(pts_pt(:,1),pts_pt(:,2),'filled')
%%show the position