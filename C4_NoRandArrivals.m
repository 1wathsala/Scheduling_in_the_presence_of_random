function C4 = C4_NoRandArrivals(ki,mi,li,lmd,mu,N,C,i,s,j,cI,cs,cr,p2)

ss = s+lmd+mu;
C4 = 0;

    if mi == 0
        C4 = C4 + 1/(s+lmd)*((mu/ss)^(ki+li))*(cI/(s+lmd) + costR(1,0,0,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2));
    else
        C4 = 0;
    end

