function [x] = vol_deg_cal_proj(a_in, D, lambda, type_D)
%g = vol_deg_cal(xi, A, D, wa_lambda, type_D); % beta
% for worest case analysis: map the current fitting error to gradients
% worest case analysis: solve \beta from A^T\xi=D^T\beta given \xi
%   min_beta 1/2||A^T\xi-D^T\beta||_2^2+lambda/2||beta||_2^2
% beta = (DD^T+ lambda I)^{-1} DA^T\alpha
a = a_in;
cg_tol = 1e-5;
cg_max_ite = 5;

switch(type_D)
    case 'graphfuse'
        s = length(a);
%         edge_num = D.edge_num;
%         node_out = D.node_out; node_in = D.node_in;
        edge_head= D.edge_head; edge_tail=D.edge_tail;
        sI = D.sI;
        xI=0; xG =0;
        ATa = a; % bI = ATa; bG = DATa;
        rI = -sI.*ATa;
        rG = ATa(edge_tail)-ATa(edge_head);
        r_prev2 = rI'*rI + rG'*rG;
        pI = -rI; 
        pG = -rG;
        DTpG = mex_graph_diffT(pG, edge_head, edge_tail, s);
        DTp = sI.*pI + DTpG;
        ApI = sI.*DTp + lambda.*pI;
%         ApG = mex_graph_diff(DTp, edge_head, edge_tail, edge_num) +lambda.*pG;
        ApG = DTp(edge_head) - DTp(edge_tail) +lambda.*pG;
        for n=1:cg_max_ite
            alpha = r_prev2 / (pI'*ApI + pG'*ApG);
            xI = xI + alpha.*pI;
            xG = xG + alpha.*pG;
            rI = rI + alpha.*ApI;
            rG = rG + alpha.*ApG;
            r2 = rI'*rI + rG'*rG;
            beta = r2 / r_prev2;
            pI = -rI + beta.*pI;
            pG = -rG + beta.*pG;
            DTpG = mex_graph_diffT(pG, edge_head, edge_tail, s);
            DTp = sI.*pI + DTpG;
            ApI = sI.*DTp + lambda.*pI;
%             ApG = mexGraphDiff(DTp, edge_head, edge_tail, edge_num) +lambda.*pG;
            ApG = DTp(edge_head) - DTp(edge_tail) +lambda.*pG;
            r_prev2 = r2;
            if(norm(r2)<cg_tol)
                break;
            end
        end
        x = [xI; xG];
    case 'flasso_1d'
        %% CG version
        x=0;
        ATa = a;
        DATa = grad_cal_1d(ATa, 'f');
        % b = xxx; r=-b;
        r = -[D.*ATa; DATa]; % here D denotes the parameter for l1 norm,
        r_prev2 = r'*r;
        p = -r;
        Ap = ADx_cal_flasso_1d(p, D, lambda);
        for n=1:cg_max_ite
            alpha = r_prev2 / (p'*Ap);
            x = x + alpha.*p;
            r = r + alpha.*Ap;
            r2 = r'*r;
            beta = r2 / r_prev2;
            p = -r + beta.*p;
            Ap = ADx_cal_flasso_1d(p, D, lambda);
            r_prev2 = r2;
            if(norm(r)<cg_tol)
                break;
            end
        end
    case 'flasso_2d'
        % to be implemented
        x1 = 0; xv = 0; xh=0;
        h = D(1); w = D(2); sI = D(3);
        ATa2d = reshape(a, h, w);
        r1 = -sI.*ATa2d;
        rv = -grad_cal_2d(ATa2d, 'v'); % b = xxx; r=-b;
        rh = -grad_cal_2d(ATa2d, 'h');
        %% acc test, to here
        r11d = r1(:); rv1d = rv(:); rh1d = rh(:);
        r_prev2 = r11d'*r11d + rv1d'*rv1d + rh1d'*rh1d;
        p1 = -r1;
        pv = -rv;
        ph = -rh;
        [Ap1, Apv, Aph] = ADx_cal_flasso_2d(p1, pv, ph, sI, lambda);
        for n=1:cg_max_ite
            alpha = r_prev2 / (p1(:)'*Ap1(:) + pv(:)'*Apv(:) + ph(:)'*Aph(:));
            x1 = x1 + alpha.*p1;
            xv = xv + alpha.*pv;
            xh = xh + alpha.*ph;
            r1 = r1 + alpha.*Ap1;
            rv = rv + alpha.*Apv;
            rh = rh + alpha.*Aph;
            r11d = r1(:); rv1d = rv(:); rh1d = rh(:);
            r2 = r11d'*r11d + rv1d'*rv1d + rh1d'*rh1d;
            beta = r2 / r_prev2;
            p1 = -r1 + beta.*p1;
            pv = -rv + beta.*pv;
            ph = -rh + beta.*ph;
            [Ap1, Apv, Aph] = ADx_cal_flasso_2d(p1, pv, ph, sI, lambda);
            r_prev2 = r2;
            if(sqrt(r2)<cg_tol)
                break;
            end
        end
        x = [x1(:); xv(:); xh(:)];
        %%
        
    case 'graph'
        node_out = D.node_out; node_in = D.node_in;
        edge_head= D.edge_head; edge_tail=D.edge_tail;
        x=0;
        ATa = a;
%         b = ATa(edge_head)-ATa(edge_tail);
        r = ATa(edge_tail)-ATa(edge_head); % b = xxx; r=-b;
        r_prev2 = r'*r;
        p = -r;
        fh = @(idx)node_deg_clc(p, idx);
        DTp = cellfun(fh, node_out) - cellfun(fh, node_in);
        Ap = DTp(edge_head) - DTp(edge_tail) + lambda.*p;
        for n=1:cg_max_ite
            alpha = r_prev2 / (p'*Ap);
            x = x + alpha.*p;
            r = r + alpha.*Ap;
            r2 = r'*r;
            beta = r2 / r_prev2;
            p = -r + beta.*p;
            Ap = DTp(edge_head) - DTp(edge_tail) + lambda.*p;
            r_prev2 = r2;
            if(norm(r)<cg_tol)
                break;
            end
        end
    
    case 'tv_1d'
        %% CG version
        x=0;
        ATa = a;
        r = -grad_cal_1d(ATa, 'f'); % b = xxx; r=-b;
        r_prev2 = r'*r;
        p = -r;
        Ap = grad_cal_1d(grad_cal_1d(p, 'b'), 'f') + lambda.*p;
        for n=1:cg_max_ite
            alpha = r_prev2 / (p'*Ap);
            x = x + alpha.*p;
            r = r + alpha.*Ap;
            r2 = r'*r;
            beta = r2 / r_prev2;
            p = -r + beta.*p;
            Ap = grad_cal_1d(grad_cal_1d(p, 'b'), 'f') + lambda.*p;
            r_prev2 = r2;
            if(norm(r)<cg_tol)
                break;
            end
        end
    case 'tv_2d'
        % to be implemented
        xv = 0; xh=0;
        h = D(1); w = D(2);
        ATa2d = reshape(a, h, w);
        rv = -grad_cal_2d(ATa2d, 'v'); % b = xxx; r=-b;
        rh = -grad_cal_2d(ATa2d, 'h');
        rv1d = rv(:); rh1d = rh(:);
        r_prev2 = rv1d'*rv1d + rh1d'*rh1d;
        pv = -rv;
        ph = -rh;
        [Apv, Aph] = ADx_cal_tv_2d(pv, ph, lambda);
        for n=1:cg_max_ite
            alpha = r_prev2 / (pv(:)'*Apv(:) + ph(:)'*Aph(:));
            xv = xv + alpha.*pv;
            xh = xh + alpha.*ph;
            rv = rv + alpha.*Apv;
            rh = rh + alpha.*Aph;
            rv1d = rv(:); rh1d = rh(:);
            r2 = rv1d'*rv1d + rh1d'*rh1d;
            beta = r2 / r_prev2;
            pv = -rv + beta.*pv;
            ph = -rh + beta.*ph;
            [Apv, Aph] = ADx_cal_tv_2d(pv, ph, lambda);
            r_prev2 = r2;
            if(sqrt(r2)<cg_tol)
                break;
            end
        end
        x = [xv(:); xh(:)];
        
        
    otherwise
        disp('Wrong type of D.\n');
        x = nan;
end
return



function [Ax] = ADx_cal(x, D, lambda, type_D)
% if(A==1)
%     ATA = 1;
% else
%     ATA = A'*A;
% end

s = length(x)/2;
% % % Ax = ATA*x + rho.*(grad_cal_1d(grad_cal_1d(x,'forward'), 'backward') + (D2).*x); % here D denotes the parameter for l1 norm
switch(type_D)
    case 'flasso_1d'
        x1 = x(1:s); x2 = x(s+1:end);
        Ax1 = (D*D).*x1 + D.*(grad_cal_1d(x2, 'b'));
        Ax2 = D.*(grad_cal_1d(x1, 'f'))+grad_cal_1d(grad_cal_1d(x2, 'b'), 'f');
        Ax = [Ax1; Ax2] + lambda.*x;
    case 'tv_1d'
        Ax = grad_cal_1d(grad_cal_1d(x, 'b'), 'f') + lambda.*x;
    otherwise
        disp('Wrong type of D.\n');
        Ax = nan;
end
return



function [Ax] = ADx_cal_flasso_1d(x, sI, lambda)
s = length(x)/2;
x1 = x(1:s); x2 = x(s+1:end);
Ax1 = (sI*sI).*x1 + sI.*(grad_cal_1d(x2, 'b'));
Ax2 = sI.*(grad_cal_1d(x1, 'f'))+grad_cal_1d(grad_cal_1d(x2, 'b'), 'f');
Ax = [Ax1; Ax2] + lambda.*x;
return

function [Ap1, Apv, Aph] = ADx_cal_flasso_2d(p1, pv, ph, sI, lambda)
tmp = sI.*p1 + grad_cal_2d(pv, 'vT') + grad_cal_2d(ph, 'hT');
Ap1 = sI.*tmp + lambda.*p1;
Apv = grad_cal_2d(tmp, 'v') + lambda.*pv;
Aph = grad_cal_2d(tmp, 'h') + lambda.*ph;
return

function [Ap] = ADx_cal_flasso_2d_comb(p, sI, lambda)
ss = size(p, 1)/3;
p1 = p(1:ss,:); pv = p(ss+1:ss*2,:); ph = p(ss*2+1:end,:);
tmp = sI.*p1 + grad_cal_2d(pv, 'vT') + grad_cal_2d(ph, 'hT');
Ap1 = sI.*tmp + lambda.*p1;
Apv = grad_cal_2d(tmp, 'v') + lambda.*pv;
Aph = grad_cal_2d(tmp, 'h') + lambda.*ph;
Ap = [Ap1; Apv; Aph];
return

function [Apv, Aph] = ADx_cal_tv_2d(pv, ph, lambda)
% Dvpv = grad_cal_2d(pv, 'vT');
% Dhph = grad_cal_2d(ph, 'hT');
% Apv = grad_cal_2d(Dvpv, 'v') + grad_cal_2d(Dhph, 'v') + lambda.*pv;
% Aph = grad_cal_2d(Dvpv, 'h') + grad_cal_2d(Dhph, 'h') + lambda.*ph;
tmp1 = grad_cal_2d(pv, 'vT') + grad_cal_2d(ph, 'hT');
Apv = grad_cal_2d(tmp1, 'v') + lambda.*pv;
Aph = grad_cal_2d(tmp1, 'h') + lambda.*ph;
return

function degr = node_deg_clc(x, idx)
degr = sum(x(idx));
return

function degr = node_deg_dif_clc(x, idx, idx2)
degr = sum(x(idx)) - sum(x(idx2));
return


% function [Ax] = ADx_cal_tv1d(x, lambda)
% Ax = grad_cal_1d(grad_cal_1d(x, 'b'), 'f') + lambda.*x;
% return



%% acc test
% %         x= 0;
% %         rcomb = [r1; rv; rh];
% %         rcomb1d = rcomb(:);
% %         r_prev2 = rcomb1d'*rcomb1d;
% %         pcomb = -rcomb;
% %         [Ap] = ADx_cal_flasso_2d_comb(pcomb, sI, lambda);
% %         for n=1:cg_max_ite
% %             alpha = r_prev2 / (pcomb(:)'*Ap(:));
% %             x = x + alpha.*pcomb;
% %             rcomb = rcomb  + alpha.*Ap;
% %             rcomb1d = rcomb(:);
% %             r2 = rcomb1d'*rcomb1d;
% %             beta = r2 / r_prev2;
% %             pcomb = -rcomb + beta.*pcomb;
% %             [Ap] = ADx_cal_flasso_2d_comb(pcomb, sI, lambda);
% %             r_prev2 = r2;
% %             if(sqrt(r2)<cg_tol)
% %                 break;
% %             end
% %         end
% %         ss = h;
% %         x1 = x(1:ss, :); xv = x(ss+1:ss*2,:); xh = x(ss*2+1:end,:);
% %         x = [x1(:); xv(:); xh(:)];




