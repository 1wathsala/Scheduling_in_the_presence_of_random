function C5= C5_NoRandsAtTheBeginnig(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2)

ss = s+lmd+mu;
C5 = 0;

    if mi == 0
        if lmd > 0 && C > 0
            K2 = (mu/ss)^(ki+li)*lmd/((s+lmd)*(s+mu))*(cr*C/(2*(s+mu)) + cI/(s+lmd) - (C+1)*p2 + costR(1,C,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2));
            for r1 = 1:C
                K1 = (mu/ss)^(ki+li)*(lmd/ss)^(r1+1)/(s+lmd)*(cr*(r1-C)*(r1+1)/(2*ss)-(r1-C)*p2 + costR(1,r1,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - costR(1,C,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2));
                C5 = C5 + K1;
            end
            C5 = C5 + K2;
        end
    end
  