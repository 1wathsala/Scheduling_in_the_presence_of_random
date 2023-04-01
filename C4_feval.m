function C44 = C4_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat)

ss = s+lmd+mu;
C44 = 0;

    if mi == 0
        C44 = C44 + 1/(s+lmd)*((mu/ss)^(ki+li))*(cI/(s+lmd) + Mat(find(Mat(:,1)==1 & Mat(:,2)==0 & Mat(:,3)==0),4));
    else
        C44 = 0;
    end
