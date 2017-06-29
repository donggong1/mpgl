function [map_out_1c] = suppmap_cal(g, kappa, map_prev_1c, is_init)
%(g, kappa, map_prev_1c, maprefine_thre, mapext_thre, rmv_rate, is_init)
% map_out_inc: 1D incremental map
% map_out: m-D map


if(is_init) % as initilization
    num = length(g);
    [~, tmpsidx] = sort(g(:),'descend');
    slct_vtmpidx = tmpsidx(1:kappa);
    % make the support map (binary indicator)
    map_out_1c = zeros(num,1); map_out_1c(slct_vtmpidx) = 1;
    %% support map refine????
%     map_out_1c = map_refine(map_out_1c, maprefine_thre);
%     map_out_1c = map_ext(map_out_1c, mapext_thre);
%     map_out_inc = map_out_1c;
else
    num = length(g);
    g(logical(map_prev_1c)) = -inf;
    [~, tmpsidx] = sort(g(:),'descend');
    slct_vtmpidx = tmpsidx(1:kappa);
    % make the support map (binary indicator)
    map_out_1c = map_prev_1c;
    map_out_1c(slct_vtmpidx) = 1;
    
%     map_out_inc = zeros(num,1); map_out_inc(slct_vtmpidx) = 1;
%     map_out_1c = min(map_prev_1c+map_out_inc, 1);
    
    %% refine
%     map_out_1c = suppmap_inactive_removal(map_out_1c, g, rmv_rate(1), rmv_rate(2));
%     map_out_1c = map_refine(map_out_1c, maprefine_thre); % map refine
%     map_out_1c = map_ext(map_out_1c, mapext_thre);
end
return