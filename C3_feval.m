function C33 = C3_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat)

ss = s+lmd+mu;
C33 = 0;

    if ki+li > 1
        for di = 1:(ki+li-1)
            K2 = (mu/(s+mu))^di*(1/(s+mu))*(cr*(C+mi)*(di+1)/(2*(s+mu)) - (C-mi)*p2 + Mat(find(Mat(:,1)==ki+li+1-di & Mat(:,2)==C & Mat(:,3)==0),4));
            for r1 = 0:C-mi
                K1 = (mu/ss)^di*(lmd/ss)^r1*factorial(di+r1)/(2*factorial(di)*factorial(r1)*ss)*(cr*(r1+mi-C)*(di+r1+1)/(2*ss) - (r1+mi-C)*p2 + Mat(find(Mat(:,1)==ki+li+1-di & Mat(:,2)==mi+r1 & Mat(:,3)==0),4) - Mat(find(Mat(:,1)==ki+li+1-di & Mat(:,2)==C & Mat(:,3)==0),4));
                C33 = C33 + K1; 
            end
            C33 = C33 + K2; 
        end
    else
        C33 = 0;
    end



