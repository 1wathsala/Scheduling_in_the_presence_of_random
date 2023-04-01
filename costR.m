function Co = costR(ki,mi,li,i,N,lmd,mu,s,j,C,cI,cs,cr,p2)
    
    Co = 0;
    
    Ci = 0;
    
    if i < N
        
        v1 = 1;
        for c = j+1:N-1
            v1 = v1*s(c);
        end

        C1 = (cs/s)*(ki+li-1)/mu;
        C2 = C2_NoServiceCompletion(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2);
        C3 = C3_SomeServiceCompletion(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2);
        C4 = C4_NoRandArrivals(ki,mi,li,lmd,mu,N,C,i,s,j,cI,cs,cr,p2);
        C5 = C5_NoRandsAtTheBeginnig(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2);
        C6 = C6_SomeRandServiceCompletion(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2);
        C7 = C7_ServeFirstBatchOfRandPlusMoreRands(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2);
        C8 = C8_ServeFirstBatchOfRandPlusNoArrivals(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2);
        C9 = 0;%C9_NoRandsAtTheBeginningMoreAtTheEnd(ki,mi,li,lmd,mu,C,N,i,s,j,cI,cs,cr,p2);
               
            
        Ci =  Ci + C1+C2+C3+C4+C5+C6+C7+C8+C9;

            
    elseif i==N
        CN = 0;
        CN = CN + mi*(ki+li)/mu;
          
        Ch = 0;
            if mi>0
                for a = 1:mi
                    Ch = Ch + (a-1)/mu;  
                end
            else
                Ch = 0;
            end
            
        CN = CN + Ch;
        Ci = Ci + (cr*CN + cs*((ki+li-1)/mu));% + (ki+mi-(1-li))/(s(1)*mu);
    else
        Ci = 0;
    end
        Co = Co + Ci;
 