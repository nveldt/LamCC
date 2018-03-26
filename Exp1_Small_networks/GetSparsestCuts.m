%% Sparsest Cuts
%
% By repeatedly running LambdaCC for smaller values of lambda, we can find
% the scaled sparsest cut. There are more efficient ways to compute scaled
% sparsest cut, but this works for small graphs such as these.
%
addpath('../Algs/')
addpath('../Networks/')
for i = 1
    
    % Load a small network
    if i ==1
        name = 'karate';
    elseif i ==2
        name = 'dolphins';
    elseif i ==3
        name = 'lesmis';
    elseif i ==4
        name = 'polbooks';
    else
        name = 'football';
    end
    
    load(strcat(name,'.mat'))
    
    A = spones(Problem.A);

    % Repeatedly call lamCC with smaller values of lambda to find the
    % minimum scaled sparsest cut
    n = size(A,1);
    T = triangle_inequality(n);
    [C,~] = lam_cc(A,.1,T);
    cut = C(:,1);
    lam = findSSC(A,cut);

    while nnz(cut) > 0 && nnz(cut) < n
        BestCut = cut;
        lam = findSSC(A,cut);
        [C,~] = lam_cc(A,lam,T);
        cut = C(:,1)';
        fprintf('%d cluster, sparsest cut %f \n',size(C,2),findSSC(A,cut))
    end
    
    ScaledSC = lam;
    save(strcat('PlotData/',name,'SparsestCut'),'BestCut','ScaledSC')
    
end
    
    