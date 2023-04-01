function C22 = C2_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat)
    
ss = s+lmd+mu;
C22 = 0;
        
    if lmd == 0
        C22=0;
    end
    
    K2 = 1/(s+mu)*(cr*(C+mi)/(2*(s+mu)) - (C-mi)*p2 + Mat(find(Mat(:,1)==ki+1 & Mat(:,2)==C & Mat(:,3)==li),4));
    for r1 = 0:C-mi
        K1 = (lmd/(ss))^r1*(1/ss)*(cr*(r1+mi-C)*(r1+1)/(2*ss) - (r1+mi-C)*p2 + Mat(find(Mat(:,1)==ki+1 & Mat(:,2)==mi+r1 & Mat(:,3)==li),4) - Mat(find(Mat(:,1)==ki+1 & Mat(:,2)==C & Mat(:,3)==li),4));
        C22 = C22 + K1;
    end
    C22 = C22 + K2;
