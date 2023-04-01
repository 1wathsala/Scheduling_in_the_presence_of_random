function fval = cost_fun_eval(t,N,C,mu,lmd,p2,p1,cI,cs,cr)
tic;
syms s;
M1 =64;

iter = (N-1).*(6*C+1);

a=0;

i = N-1;
j = N-1;

if N==2
    fun1 = matlabFunction(costR(1,0,0,i,N,lmd,mu,s,j,C,cI,cs,cr,p2) + 100/s);
    
    %fval = talbot_inversion(fun1, t(i), M1) - 100 - N*p1;        
    fval = inv_cme(fun1, t(i), 1000) - 100 - N*p1;

else
    
Mat = zeros(iter,4);
for ki = 1:i
    for li = 0:1
        for mi = 0:C
   
            a=a+1;
            %a=ki*(6*C+1)-6*C+li*(3*C+1)+mi; 
            
            Mat(a,1) = ki;
            Mat(a,2) = mi;
            Mat(a,3) = li;
            
            fun1 = matlabFunction(costR(ki,mi,li,i,N,lmd,mu,s,j,C,cI,cs,cr,p2) + 100/s);
            
            %f1 = talbot_inversion(fun1, t(i), M1) - 100;
            f1 = inv_cme(fun1, t(i), 1000) - 100;
            
            Mat(a,4) = f1;  
  
        end
    end
end

Mat;
a=0;

if N > 3
    for z = 2:N-2
        i=N-z;
    
        iter = i*(6*C+1);
        k = zeros(iter,1);
        m = zeros(iter,1);
        l = zeros(iter,1);
        ofv = zeros(iter,1);
    
    for ki = 1:i
        for li = 0:1
            for mi = 0:C
            
                a=a+1;
                %a=ki*(6*C+1)-6*C+li*(3*C+1)+mi    
                k(a) = ki;
                m(a) = mi;
                l(a) = li;
            
                fun = matlabFunction(costQ(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat) + 100/s);
                
                %f = talbot_inversion(fun, t(i), M1) - 100;
                f = inv_cme(fun, t(i), 1000) - 100;
            
                ofv(a) = f;
  
  
            end
        end
    end
    Mat = zeros(4);
    Mat =  [k m l ofv];
    end
end

i=1;
ki=1;mi=0;li=0;

funf = matlabFunction(costQ(ki,mi,li,lmd,mu,C,s,cI,cs,cr,p2,Mat) + 100/s);

%ff = talbot_inversion(funf, t(i), M1) - 100;
ff = inv_cme(funf, t(i), 1000) - 100 ;

fval = ff - N*p1;
end
time = toc;

