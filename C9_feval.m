function C99 = C9_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat)

ss = s+lmd+mu;
C99 = 0;



    if mi == 0
        if lmd > 0
            for q = 0:C-1
            K4 = (mu/ss)^(ki+li)*lmd/(s+lmd)*mu/(s+mu)*mu^q/(s+mu)^(q+1)*(cI/(s+lmd) + cr*C/(2*(s+mu)) + cr*(2*C-q+1)*(q+1)/(2*(s+mu)) - p2*(C+2) + Mat(find(Mat(:,1)==1 & Mat(:,2)==C-q & Mat(:,3)==1),4));
            for r3 = 0:1
                K3 = (mu/ss)^(ki+li)*lmd/(s+lmd)*mu/(s+mu)*mu^q*lmd^r3*factorial(q+r3)/(factorial(q)*factorial(r3)*(ss))^(q+r3+1)*(cI/(s+lmd) + cr*C/(2*(s+mu)) + cr*(r3+2*C-q-1)*(q+r3+1)/(2*(ss)) - p2*(C+r3+1) + Mat(find(Mat(:,1)==1 & Mat(:,2)==C+r3-q-1 & Mat(:,3)==1),4));
                C99 = C99 + K3;
            end
            C99 = C99 + K4;
            end
            for r2 = 0:C-1
                if r2 > 0
                    for q = 0:r2-1
                        K2 = (mu/ss)^(ki+li)*lmd/(s+lmd)*mu*lmd^r2/ss^(r2+1)*mu^q/(s+mu)^(q+1)*(cI/(s+lmd) + cr*(r2-C)*(r2+1)/(2*(ss)) + cr*(2*r2-q+1)*(q+1)/(2*(s+mu)) - p2*(r2+2) + Mat(find(Mat(:,1)==1 & Mat(:,2)==r2-q & Mat(:,3)==1),4));
                        for r3 = 0:1
                            K1 = (mu/ss)^(ki+li)*lmd/(s+lmd)*mu*lmd^r2/ss^(r2+1)*mu^q*lmd^r3*factorial(q+r3)/(factorial(q)*factorial(r3)*ss^(q+r3+1))*(cI/(s+lmd) + cr*(r2-C)*(r2+1)/(2*ss) + cr*(r3+2*r2-q-1)*(q+r3+1)/(2*ss) - p2*(r2+r3+1) + Mat(find(Mat(:,1)==1 & Mat(:,2)==r2+r3-q-1 & Mat(:,3)==1),4));
                            C99 = C99 + K1;
                        end
                        C99 = C99 + K2;
                    end
                end
            end
        end
    end 