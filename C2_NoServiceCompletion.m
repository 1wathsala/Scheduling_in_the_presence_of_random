function C2 = C2_NoServiceCompletion(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2)
    
ss = s+lmd+mu;
C2 = 0;
    
    if lmd == 0
        C2=0;
    end
    
    K2 = 1/(s+mu)*(cr*(C+mi)/(2*(s+mu)) - (C-mi)*p2 + costR(ki+1,C,li,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2));
    for r1 = 0:C-mi
        K1 = (lmd/(ss))^r1*(1/ss)*(cr*(r1+mi-C)*(r1+1)/(2*ss) - (r1+mi-C)*p2 + costR(ki+1,mi+r1,li,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - costR(ki+1,C,li,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2));
        C2 = C2 + K1;
    end
    C2 = C2 + K2;
