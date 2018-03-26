%% Find if there is a better sparsest cut

name = 'ca-GrQc';

load(strcat(name,'.mat'))
A = spones(Problem.A);
A = A-diag(diag(A));
A = StandardizeFully(A);

load ca-GrQc_expansion

%%

alpha = findSSC(A,bset);

[found,bset,ssc] = BetterSSC(A,alpha,bset)

%%

load BestsscSofar_cagrqc

alpha = findSSC(A,bset);

Lams = [1e-5 .0000278, .000077]

objs = LamRange_objs(A,bset+1,Lams)

n = size(A,1)
objs2 = LamRange_objs(A,ones(n,1),Lams)