function C3 = C3_SomeServiceCompletion(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2)

ss = s+lmd+mu;
C3 = 0;

    if ki+li > 1
        for di = 1:(ki+li-1)
            K2 = (mu/(s+mu))^di*(1/(s+mu))*(cr*(C+mi)*(di+1)/(2*(s+mu)) - (C-mi)*p2 + costR(ki+li+1-di,C,0,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2));
            for r1 = 0:C-mi
                K1 = (mu/ss)^di*(lmd/ss)^r1*factorial(di+r1)/(2*factorial(di)*factorial(r1)*ss)*(cr*(r1+mi-C)*(di+r1+1)/(2*ss) - (r1+mi-C)*p2 + costR(ki+li+1-di,mi+r1,0,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - costR(ki+li+1-di,C,0,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2));
                C3 = C3 + K1; 
            end
            C3 = C3 + K2; 
        end
    else
        C3 = 0;
    end



