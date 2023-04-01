function Ct = costQ(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat)

L1 = (cs/s)*(ki+li-1)/mu;
L2 = C2_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat);
L3 = C3_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat);
L4 = C4_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat);
L5 = C5_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat);
L6 = C6_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat);
L7 = C7_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat);
L8 = C8_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat);
L9 = C9_feval(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat);

Ct = L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9;