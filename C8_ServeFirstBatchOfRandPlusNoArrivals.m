function C8 = C8_ServeFirstBatchOfRandPlusNoArrivals(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2)

ss = s+lmd+mu;
C8 = 0;

    if lmd >0 && C >0
    v2 = costR(1,0,0,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2);
    K3 = (mu/(s+mu))^(ki+li)*(mu/ss)^C*(1/(s+lmd))*(cr*(C+mi)*(ki+li)/(2*(s+mu)) + cr*C^2/(2*ss) + cI/(s+lmd) - (C-mi)*p2 + v2);
    C8 = C8 + K3;
    for r1 = 0:C-mi
        K2 = (mu/ss)^(ki+li+C)*(lmd/ss)^r1*(1/(s+lmd))*(factorial(ki+li+r1-1)/(factorial(ki+li-1)*factorial(r1)))*(cr*(C+mi)*(ki+li+r1)/(2*ss) + cr*C^2/(2*ss) + cI/((s+lmd)) - (C-mi)*p2 + v2);
        C8 = C8 - K2;
        if mi+r1 > 0
            K1 = (mu/ss)^(ki+li+mi+r1)*(lmd/ss)^r1*(1/(s+lmd))*(factorial(ki+li+r1-1)/(factorial(ki+li-1)*factorial(r1)))*(cr*(2*mi+r1)*(ki+li+r1)/(2*ss) + cr*(mi+r1)^2/(2*ss) + cI/(s+lmd) - r1*p2 + v2);
            C8 = C8 + K1;
        end
    end
    else
        C8 = 0;
    end

    