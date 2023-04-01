function C9= C9_NoRandsAtTheBeginningMoreAtTheEnd(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2)

ss = s+lmd+mu;
C9 = 0;

    if mi == 0
        if lmd > 0
            for q = 0:C-1
            K4 = (mu/ss)^(ki+li)*lmd/(s+lmd)*mu/(s+mu)*mu^q/(s+mu)^(q+1)*(cI/(s+lmd) + cr*C/(2*(s+mu)) + cr*(2*C-q+1)*(q+1)/(2*(s+mu)) - p2*(C+2) + costR(1,C-q,1,i+1,N,lmd,mu,s,j+1,C,cI,cs,cr,p2));
            for r3 = 0:1
                K3 = (mu/ss)^(ki+li)*lmd/(s+lmd)*mu/(s+mu)*mu^q*lmd^r3*factorial(q+r3)/(factorial(q)*factorial(r3)*(ss))^(q+r3+1)*(cI/(s+lmd) + cr*C/(2*(s+mu)) + cr*(r3+2*C-q-1)*(q+r3+1)/(2*ss) - p2*(C+r3+1) + costR(1,C+r3-q-1,1,i+1,N,lmd,mu,s,j+1,C,cI,cs,cr,p2));
                C9 = C9 + K3;
            end
            C9 = C9 + K4;
            end
            for r2 = 0:C-1
                for q = 0:r2-1
                    K2 = (mu/ss)^(ki+li)*lmd/(s+lmd)*mu*lmd^r2/ss^(r2+1)*mu^q/(s+mu)^(q+1)*(cI/(s+lmd) + cr*(r2-C)*(r2+1)/(2*ss) + cr*(2*r2-q+1)*(q+1)/(2*(s+mu)) - p2*(r2+2) + costR(1,r2-q,1,i+1,N,lmd,mu,s,j+1,C,cI,cs,cr,p2));
                    for r3 = 0:1
                        K1 = (mu/ss)^(ki+li)*lmd/(s+lmd)*mu*lmd^r2/ss^(r2+1)*mu^q*lmd^r3*factorial(q+r3)/(factorial(q)*factorial(r3)*ss^(q+r3+1))*(cI/(s+lmd) + cr*(r2-C)*(r2+1)/(2*ss) + cr*(r3+2*r2-q-1)*(q+r3+1)/(2*ss) - p2*(r2+r3+1) + costR(1,r2+r3-q-1,1,i+1,N,lmd,mu,s,j+1,C,cI,cs,cr,p2));
                        C9 = C9 + K1;
                    end
                    C9 = C9 + K2;
                end
            end
        end
    end