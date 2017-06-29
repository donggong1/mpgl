function [x, supp_map, time] = mpgenlasso_proj_mex(y, D, edge_num, ...
    edge_head, edge_tail, sI, lambda, opts)
%% using mex for subproblem
%% specifically optimized for A=I, solve flsa (projection problem).
% opts.type_A
% opts.type_D
%% prepare parameters
type_D = opts.type_D;
rho = opts.rho;
kappa = opts.kappa;
wa_lambda = opts.wa_lambda;
ite_max_mp = opts.ite_max_mp;
ite_max_inner = opts.ite_max_inner;
ite_max_inner_final = opts.ite_max_inner_final;
x_diff_tol_inner = opts.x_diff_tol_inner;
x_diff_tol_mp = opts.x_diff_tol_mp;
cg_max_ite_inner = opts.cg_max_ite_inner;
cg_max_ite_inner_final = opts.cg_max_ite_inner_final;
tic;
%% initialization
s = length(y);
% x init
x = zeros(s,1);
% support map init
xi = y-x;
g = vol_deg_cal_proj(xi, D, wa_lambda, type_D); % beta
gabs = abs(g);
if(kappa == -1)
    kappa_rate = opts.kappa_rate;
    kappa = length(find(gabs>kappa_rate*max(gabs)));
end
[supp_map] = suppmap_cal(gabs, kappa, [], true);
ATy = y;
%% to find a good starting point
%     x = y;
x_prev = x;
%% prepare for mex admm
% edge_num = D.edge_num;
% edge_head= int32(D.edge_head-1); 
% edge_tail= int32(D.edge_tail-1);
% sI=D.sI;
supp_map_idx = int32(find(supp_map)-1);
non_zero_num = length(supp_map_idx);
x_diff_tol_inner2=x_diff_tol_inner.^2;
for ite_mp = 1:ite_max_mp
    [x] = mex_subadmm_proj(x, ATy, supp_map_idx, non_zero_num, ...
        edge_head, edge_tail, edge_num, s, sI, lambda, rho, ...
        ite_max_inner, x_diff_tol_inner2, cg_max_ite_inner, 1e-10);      
    %% update status
    x_diff = x - x_prev;
    x_diff_rel = norm(x_diff)/norm(x_prev);
    %% checking stopping condition
    if(x_diff_rel<x_diff_tol_mp)
        break;
    end
    x_prev = x; 
    %% update supp map
    xi = y-x;
    g = vol_deg_cal_proj(xi, D, wa_lambda, type_D); % beta
    gabs = abs(g);
% % %     kappa = min(kappa/2,1);
    [supp_map] = suppmap_cal(gabs, kappa, supp_map, false);
    supp_map_idx = int32(find(supp_map)-1);
    non_zero_num = length(supp_map_idx);
end
%% final tuning
% fprintf('final tuning\n');
[x] = mex_subadmm_proj(x, ATy, supp_map_idx, non_zero_num, ...
        edge_head, edge_tail, edge_num, s, sI, lambda, rho, ...
        ite_max_inner_final, x_diff_tol_inner2, cg_max_ite_inner_final, 1e-10); 

time=toc;
return





