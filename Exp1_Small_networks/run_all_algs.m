function [Bounds, LPs, LPtime,Objectives, Times] = run_all_algs(A,Lams,ExactFlag)
% This takes an input graph, a range of lambda values, and runs the
% following set of algorithms the network for all values of lambda:
%
%   1. ExactILP in Gurobi (only if ExactFlag == 1)
%   2. The 4-approximation of Charikar et al. (mainly intended for lam =
%   1/2)
%   3. ThreeLP, guaranteed to be a 3-approximation if lambda >= .5
%   4. ICM, from Bagon and Galum 2011
%   5. Grow-Cluster
%   6. Lambda-Louvain
%   7. FiveLP, the constant factor approximation we developed prior to
%       ThreeLP


addpath('../large_scale_cc-master')
addpath('../Algs')
addpath('../Algs/Lambda_louvain/')
n = size(A,1);
T = triangle_inequality(n);
r = numel(Lams);
LPs = zeros(n);
    fprintf(' %8s %8s %8s %8s %8s %8s %8s %8s  %8s\n','Lambda','Pivot','GrowCluster','CGW','LP5','LP-piv','BG11','OPT','Bound')

for i = 1:r

    lam = Lams(i);
    if ExactFlag
        tic
        [~,~,obj] = lam_cc(A,lam,T);
        ExactTime(i) = toc;
        Exact(i) = obj;
    else
        ExactTime(i) = 0;
        Exact(i) = 0;
    end
    
    tic
    [dist,bound] = lam_cc_LPrelax(A,lam,T);
    LPtime(i) = toc;
    LPs(:,:,i) = dist;      % Keeping track of the LP relaxations
    Bounds(i) = bound;      % gives us a lower bound even when we don't 
                            % solve the problem exactly.
                                  
    tic
    [~,LPCGW] = CGW_LP_round(A-lam,dist,1);
    CGWTime(i) = toc;
    CGW(i) = LPCGW;

    tic
    [~,LPlam] = LP5_round(A-lam,dist,1);
   	LP5Time(i) = toc;
    LP5(i) = LPlam;
    
    tic
    [~,LP3] = ThreeLP_round(A-lam,dist,30);
   	LP3Time(i) = toc;
    LPthree(i) = LP3;

    tic
    [~,piv] = pivots(A-lam,30);
    PivTime(i) = toc;
    Piv(i) = piv;

    tic
    [~,gc] = GrowCluster(A,lam,5);
    GCTime(i) = toc;
    GC(i) = gc;
    
    tic
    cll = lambda_louvain(A,lam,ones(n,1));
    LLtime(i) = toc;
    LL(i) = CCminDisagreeObj((A-lam),cll);
  
    % Choose to run one algorithm from Bagon & Galun 2011
    % This only operates on sparse matrices
    B = sparse(A-lam);
    tic
    %bg11 = ab_swap(B);
    bg11 = AL_ICM(B);
    BG11Time(i) = toc;
    BG11(i) = CCminDisagreeObj((A-lam),bg11);
    
    fprintf('%8.3f %8.2f %8.2f %8.2f %8.2f %8.2f  %8.2f %8.2f %8.2f \n',lam,piv,gc,LPCGW,LPlam,LP3,BG11(i),Exact(i),bound)

end
Times = [ExactTime' PivTime' GCTime' (CGWTime'+LPtime') (LP5Time'+LPtime') (LP3Time'+LPtime') BG11Time' LLtime'];
Objectives = [Exact' Piv' GC' CGW' LP5' LPthree' BG11' LL'];