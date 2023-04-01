tic



s1 = sym('s1');
t1 = sym('t1');
s = [s1];
t = [t1]; 

x = [];
N = 2;
mu = 3;
lmd = 1;
p1 = 1;
p2 = 1;
C = 5;
cI = 0;
cs = 0;
cr = 0;

%fun0 = costR(1,0,0,1,N,lmd,mu,s,1,C,cI,cs,cr,p2);
fun = costR(1,0,0,1,N,lmd,mu,s,1,C,cI,cs,cr,p2);
%idlecost = 1/(s(1)^2*s(2))+1/(s(1)*s(2)^2);



%for i = 1:length(s)
 %   idlecost = idlecost + 1/s(i)^2;
%end
%fun = fun + idlecost;

for i = 1:length(s)
    fun = ilaplace(fun,s(i),t(i));
end
fun = fun - p1 * N;
fun1 = matlabFunction(fun,'Vars',{[t1]});

x0 = [0.1];
A = [-1];
b = [0];
Aeq = [0];
beq = 0;
lb = [0];
ub = [1];

[x,fval] = fmincon(fun1,x0,A,b,Aeq,beq,lb,ub);
x
fval


f = matlabFunction(fun);

fplot(f,[0,1])

%{
figure(1)
scatter(lmd,xopt(1),'r','filled')
hold on
%scatter(p2,xopt(2),'b','filled')
%hold on
%}

%x1 = fmincon(fun1,x0,A,b)
%x2 = gamultiobj(fun1,2,A,b)
%x3 = ga(fun1,2,A,b,Aeq,beq)
%x4 = particleswarm(fun1,2,lb,ub)
%x5 = simulannealbnd(fun1,x0,lb,ub)
%x3 = paretosearch(fun1,1,A,b)


%f1 = double(subs(subs(fun,t1,max([x1(1),0])),t2,max([x1(1),0])))
%f2 = double(subs(subs(fun,t1,max([x2(1),0])),t2,max([x2(1),0])))
%f3 = double(subs(subs(fun,t1,max([x3(1),0])),t2,max([x3(1),0])))
%f4 = double(subs(subs(fun,t1,max([x4(1),0])),t2,max([x4(1),0])))
%f5 = double(subs(subs(fun,t1,max([x5(1),0])),t2,max([x5(1),0])))
%f = min([f1,f2])

%if f1 == f
 %   min_obj_fun_val = f
  %  sol = x1
%elseif f2 == f
 %   min_obj_fun_val = f
  %  sol = x2
%elseif f3 == f
    %min_obj_fun_val = f
    %sol = x3
%elseif f4 == f
    %min_obj_fun_val = f
    %sol = x4
%elseif f5 == f
    %min_obj_fun_val = f
    %sol = x2
%end

%fprintf('t=%.4f \n M=%.4f \n t1=%.4f \n t2=%.4f \n t3=%.4f \n l=%.4f \n',t,M,x(1),x(2),x(3),lmd)
%{
scatter(cs,x,'b')
hold on;
%scatter(lmd,x2(2))
%hold on;
%scatter(lmd,sum(x2))
%hold on;
grid on
xlabel('c_s (unit cost of waiting for scheduled customers)');
ylabel('optimum inter-arrival times');
hold on;
%}

%grid on
%xlabel('inter-arrival time')
%ylabel('cost')
%title(['Unit revenue (random)=p_2=' num2str(p2)])
%dim = [0.2 0.5 0.3 0.3];
%str = {['xopt = ' num2str(xopt)],['fopt = ' num2str(fopt)]};
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
%hold on


%{
nexttile
f = matlabFunction(fun);
x = 0:0.01:4;  % define range and mesh of x and y which will be shown in figure
y = 0:0.01:4;
[X, Y] = meshgrid(x, y);
surf(X, Y, f(X,Y));
grid on
xlabel('x1')
ylabel('x2')
zlabel('cost')
title(['p_2 = ' num2str(p2)])
dim = [0.2 0.5 0.3 0.3];
str = {['xopt = ' num2str(xopt)],['fopt = ' num2str(fopt)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold on
%double(subs(subs(fun,t1,0.1326),t2,0.1319))
%}

%dim = [0.2 0.5 0.3 0.3];
%str = {['N = ' num2str(N)],['C = ' num2str(C)],['\mu = ' num2str(mu)],['p1 = ' num2str(p1)],['p2 = ' num2str(p2)],['cI = ' num2str(cI)],['cs = ' num2str(cs)],['cr = ' num2str(cr)]};
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
%legend()])
%xlabel('\lambda');
%ylabel('Optimal inter-arrival time');
%title('Optimal inter-arrival times variation with \lambda');
%legend(['p2 = ' num2str(p2)]);

%hold on

%print -depsc scheduletwo_cs


