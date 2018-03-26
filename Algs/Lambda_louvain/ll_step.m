function c = ll_step(A,lam,w,itlim)

n = size(A,1);
c = (1:n)';
improved = 1;
its = 0;
while improved && its < itlim
    its = its+1;
    improved = 0;
    
    for i = 1:n
        Ci = c(i);  
        Bestmove = Ci;
        Bestdecrease = 0;
        Ni = find(A(i,:));
        I = setdiff(find(c == Ci),i);
        EdgesitoI = sum(A(i,I));
        
        % w(I) is the set of node weights for I
        % At the first iteration sum(w(I)) equals |I|.
        % The number of edges from i to I is bounded above by sum(w(I))
        NonedgesitoI = w(i)*sum(w(I)) - EdgesitoI;
        assert(NonedgesitoI >= 0);
        
        % The change if i were to leave I
        deltaI = (1-lam)*EdgesitoI - lam*NonedgesitoI;
        
        for jInd = 1:numel(Ni)       % look at all neighbors of i
            j = Ni(jInd);
            Cj = c(j);
            J = find(c == Cj);
            EdgesitoJ = sum(A(i,J));
            NonedgesitoJ = w(i)*sum(w(J)) - EdgesitoJ;
            assert(NonedgesitoJ >= 0);
            
            % The change if i were to join J
            deltaJ = -(1-lam)*EdgesitoJ + lam*NonedgesitoJ;
            delta = deltaI + deltaJ;
            
            if delta < Bestdecrease
                Bestdecrease = delta;
                Bestmove = Cj;
                improved = 1;
            end
        end
        c(i) = Bestmove;
        %fprintf('Node %d moved to cluster %d \n',i,Bestmove)
        %Lccscore = Lccscore + Bestdecrease;
    end
end
c = renumber(c);

end