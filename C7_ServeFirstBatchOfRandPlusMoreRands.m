function C7 = C7_ServeFirstBatchOfRandPlusMoreRands(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2)

ss = s+lmd+mu;
ss2 = s+mu;
C7 = 0;
    
    if lmd > 0 && C > 0
        v2 = costR(1,C,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2);
        a1 = (mu/ss2)^(ki+li)*mu^C/ss2;
        K18 = a1/ss2^C*(cr*(C+mi)*(ki+li)/(2*ss2) + cr*(C+2)*C/(2*ss2) - (2*C-mi+1)*p2 + v2);
        K14 = a1*C/(ss^(C+1))*(cr*(C+mi)*(ki+li)/(2*ss2) + cr*(C+1)^2/(2*ss) + cr*C/(2*ss2) - (2*C-mi+1)*p2 + v2);
        C7 = C7 + K18 + K14;
    
        for r3 = 0:C
            a2 = (mu/ss2)^(ki+li)*(lmd/ss)^r3*(1/ss);
            K17 = a2*(mu/ss2)^C*(cr*(r3-C)*(r3+1)/(2*ss) - (r3-C)*p2 + costR(1,r3,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v2);
            K13 = a2*(mu/ss)^C*lmd*C/ss*(cr*(r3-C)*(r3+1)/(2*ss) - (r3-C)*p2 + costR(1,r3,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v2);
            C7 = C7 + K17 + K13;
        end
    
        for r2 = 0:1
            a3 = (mu/ss2)^(ki+li)*(mu/ss)^C*(lmd/ss)^r2*(factorial(C+r2-1)/factorial(C-1));
            K16 = a3/ss2*(cr*(C+mi)*(ki+li)/(2*ss2) + cr*(C+1)*(C+r2)/(2*ss) + cr*C/(2*ss2) - (2*C-mi+1)*p2 + v2);
            C7 = C7 - K16;
        
            for r3 = 0:C
                K15 = a3*lmd^r3/(ss^(r3+1))*(cr*(r3-C)*(r3+1)/(2*ss) - (r3-C)*p2 + costR(1,r2+r3-1,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v2);
                C7 = C7 - K15;
            end
        end
    
        for r1 = 0:C-mi
            a4 = (mu/ss)^(ki+li)*mu^C*(lmd/ss)^r1*(1/ss2)*(factorial(ki+li+r1-1)/(factorial(ki+li-1)*factorial(r1)));
            K12 = a4/(ss2^C)*(cr*(C+mi)*(ki+li+r1)/(2*ss) + cr*(C+1)*C/(2*ss2) + cr*C/(2*ss2) - (2*C-mi+1)*p2 + v2);
            C7 = C7 - K12;
            K8 = a4*C/(ss^(C+1))*(cr*(C+mi)*(ki+li+r1)/(2*ss) + cr*(C+1)^2/(2*ss) + cr*C/(2*ss2) - (2*C-mi+1)*p2 + v2);
            C7 = C7 - K8;
        
            if mi+r1 > 0
                for q = 0:C-mi-r1
                    v3 = costR(1,C-q,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2);
                    a5 = (mu/ss)^(ki+li+r1)*(mu/ss2)^mi*mu^q*(lmd/ss2)^r1*(factorial(ki+li+r1-1)/(factorial(ki+li-1)*factorial(r1)));
                    K6 = a5/(ss2^(q+1))*(cr*(2*mi+r1)*(ki+li+r1)/(2*ss) + cr*(C+1)*(mi+r1)/(2*ss2) + cr*(2*C-mi-r1-q-1)*(q+1)/(2*ss2) - (C+r1+1)*p2 + v3);
                    C7 = C7 + K6;
                
                    for r3 = 0:mi+r1
                        K5 = a5*lmd^r3*factorial(q+r3)/(factorial(q)*factorial(r3)*ss^(q+r3+1))*(cr*(r3-mi-r1)*(q+r3+1)/(2*ss) - (r3-mi-r1)*p2 + costR(1,C+r3-mi-r1-q,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v3);
                        C7 = C7 + K5;
                    end
                
                    for r2 = 1:C-mi-r1+1
                        a6 = (mu/ss)^(ki+li+mi+r1)*mu^q*(lmd/ss)^(r1+r2)*(factorial(ki+li+r1-1)*factorial(mi+r1+r2-1)/(factorial(mi+r1-1)*factorial(ki+li-1)*factorial(r1)*factorial(r2)));
                        K4 = a6/(ss2^(q+1))*(cr*(2*mi+r1)*(ki+li+r1)/(2*ss) + cr*(C+1)*(mi+r1+r2)/(2*ss) + cr*(2*C-mi-r1-q-1)*(q+1)/(2*ss2) - (C+r1+1)*p2 + v3);
                        C7 = C7 - K4;
                    
                        for r3 = 0:mi+r1
                        K3 = a6*lmd^r3*factorial(q+r3)/(factorial(q)*factorial(r3)*ss^(q+r3+1))*(cr*(r3-mi-r1)*(q+r3+1)/(2*ss) - (r3-mi-r1)*p2 + costR(1,C+r3-mi-r1-q,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v3);
                        C7 = C7 - K3;
                        end 
                    end
                end
            
                for r2 = 1:C-mi-r1+1
                    for q = 0:r2-1
                        v3 = costR(1,C-q,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2);
                        a7 = (mu/ss)^(ki+li+mi+r1)*mu^q*(lmd/ss)^(r1+r2)*(factorial(ki+li+r1-1)*factorial(mi+r1+r2-1)/(factorial(mi+r1-1)*factorial(ki+li-1)*factorial(r1)*factorial(r2)));
                        K2 = a7/(ss2^(q+1))*(cr*(2*mi+r1)*(ki+li+r1)/(2*ss) + cr*(mi+r1+r2)^2/(2*ss) + cr*(C+r2-q-2)*(q+1)/(2*ss2) - (C+r1+1)*p2 + v3);
                        C7 = C7 + K2;
                    
                        for r3 = 0:C-r2+1
                        K1 = a7*lmd^r3*factorial(q+r3)/(factorial(q)*factorial(r3)*ss^(q+r3+1))*(cr*(r2+r3-C-1)*(q+r3+1)/(2*ss) - (r2+r3-C-1)*p2 + costR(1,r2+r3-q-1,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v3);
                        C7 = C7 + K1;
                        end 
                    end
                end
            end
        
            for r3 = 0:C
            a8 = (mu/ss)^(ki+li)*mu^C*(lmd/ss)^(r1+r3)*(1/ss)*(factorial(ki+li+r1-1)/(factorial(ki+li-1)*factorial(r1)));
            K11 = a8/(ss2^C)*(cr*(r3-C)*(r3+1)/(2*ss) - (r3-C)*p2 + costR(1,r3,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v2);
            C7 = C7 - K11;
            K7 = a8*lmd*C/(ss^(C+1))*(cr*(r3-C)*(r3+1)/(2*ss) - (r3-C)*p2 + costR(1,r3,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v2);
            C7 = C7 - K7;
            end
        
            for r2 = 0:1
                a9 = (mu/ss)^(ki+li+C)*(lmd/ss)^(r1+r2)*(factorial(ki+li+r1-1)*factorial(C+r2-1)/(factorial(C-1)*factorial(ki+li-1)*factorial(r1)));
                K10 = a9/(ss2)*(cr*(C+mi)*(ki+li+r1)/(2*ss) + cr*(C+1)*(C+r2)/(2*ss) + cr*C/(2*ss2) - (2*C-mi+1)*p2 + v2);
                C7 = C7 + K10;
                for r3 = 0:C
                    K9 = a9*lmd^r3/(ss^(r3+1))*(cr*(r3-C)*(r3+1)/(2*ss) - (r3-C)*p2 + costR(1,r2+r3-1,1,i+1,N,lmd,mu,s,N-1,C,cI,cs,cr,p2) - v2);
                    C7 = C7 + K9;
                end
            end
        end
    else
        C7 = 0;
    end