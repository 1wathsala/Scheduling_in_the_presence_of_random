function C66 = C6_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat)

C66 = 0;
ss = s+lmd+mu;

    if lmd > 0 && C > 0
        for q = 0:C-1
            K6 = (mu/(s+mu))^(ki+li+q)*(1/(s+mu))*(cr*(C+mi)*(ki+li)/(2*(s+mu)) + cr*(2*C-q-1)*(q+1)/(2*(s+mu)) - (C-mi+1)*p2 + Mat(find(Mat(:,1)==1 & Mat(:,2)==C-q & Mat(:,3)==1),4));
            for r2 = 0:1
                K5 = (mu/(s+mu))^(ki+li)*(mu/ss)^q*(lmd/ss)^r2*factorial(q+r2)/(factorial(q)*factorial(r2)*ss)*(cr*(r2-1)*(q+r2+1)/(2*ss) - (r2-1)*p2 + Mat(find(Mat(:,1)==1 & Mat(:,2)==C+r2-q-1 & Mat(:,3)==1),4) - Mat(find(Mat(:,1)==1 & Mat(:,2)==C-q & Mat(:,3)==1),4));
                C66 = C66 + K5;
            end
            C66 = C66 + K6;
        end
        for r1 = 0:C-mi
            for q = 0:C-1
                K4 = (mu/ss)^(ki+li)*(mu/(s+mu))^q*(lmd/ss)^r1*factorial(ki+li+r1-1)/(factorial(ki+li-1)*factorial(r1)*(s+mu))*(cr*(C+mi)*(ki+li+r1)/(2*ss) + cr*(2*C-q-1)*(q+1)/(2*(s+mu)) - (C-mi+1)*p2 + Mat(find(Mat(:,1)==1 & Mat(:,2)==C-q & Mat(:,3)==1),4));
                for r2 = 0:1
                    K3 = (mu/ss)^(ki+li+q)*(lmd/ss)^(r1+r2)*factorial(ki+li+r1-1)*factorial(q+r2)/(factorial(q)*factorial(ki+li-1)*factorial(r1)*factorial(r2)*ss)*(cr*(r2-1)*(q+r2+1)/(2*ss) - (r2-1)*p2 + Mat(find(Mat(:,1)==1 & Mat(:,2)==C+r2-q-1 & Mat(:,3)==1),4) - Mat(find(Mat(:,1)==1 & Mat(:,2)==C-q & Mat(:,3)==1),4));
                    C66 = C66 - K3;
                end
                C66 = C66 - K4;
            end
            if mi+r1 > 0
                for q = 0:(mi+r1-1)
                    K2 = (mu/ss)^(ki+li)*(mu/(s+mu))^q*(lmd/ss)^r1*factorial(ki+li+r1-1)/(factorial(ki+li-1)*factorial(r1)*(s+mu))*(cr*(2*mi+r1)*(ki+li+r1)/(2*ss) + cr*(C+mi+r1-q-1)*(q+1)/(2*(s+mu)) - (C-mi+1)*p2 + Mat(find(Mat(:,1)==1 & Mat(:,2)==C-q & Mat(:,3)==1),4));
                    for r2 = 0:C-mi-r1+1
                        K1 = (mu/ss)^(ki+li+q)*(lmd/ss)^(r1+r2)*factorial(ki+li+r1-1)*factorial(q+r2)/(factorial(q)*factorial(ki+li-1)*factorial(r1)*factorial(r2)*ss)*(cr*(mi+r1+r2-C-1)*(q+r2+1)/(2*ss) - (mi+r1+r2-C-1)*p2 + Mat(find(Mat(:,1)==1 & Mat(:,2)==mi+r1+r2-q-1 & Mat(:,3)==1),4) - Mat(find(Mat(:,1)==1 & Mat(:,2)==C-q & Mat(:,3)==1),4));
                        C66 = C66 + K1;
                    end
                    C66 = C66 + K2;
                end
            else
                C66 = 0;
            end
        end
    else
        C66 = 0;
    end

