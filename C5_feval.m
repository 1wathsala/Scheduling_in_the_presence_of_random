function C55= C5_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat)

ss = s+lmd+mu;
C55 = 0;

    if mi == 0
        if lmd > 0 && C > 0
            K2 = (mu/ss)^(ki+li)*lmd/((s+lmd)*(s+mu))*(cr*C/(2*(s+mu)) + cI/(s+lmd) - (C+1)*p2 + Mat(find(Mat(:,1)==1 & Mat(:,2)==C & Mat(:,3)==1),4));
            for r1 = 1:C
                K1 = (mu/ss)^(ki+li)*(lmd/ss)^(r1+1)/(s+lmd)*(cr*(r1-C)*(r1+1)/(2*ss)-(r1-C)*p2 + Mat(find(Mat(:,1)==1 & Mat(:,2)==r1 & Mat(:,3)==1),4) - Mat(find(Mat(:,1)==1 & Mat(:,2)==C & Mat(:,3)==1),4));
                C55 = C55 + K1;
            end
            C55 = C55 + K2;
        end
    end
