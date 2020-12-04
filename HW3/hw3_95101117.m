%% Problem 5 part a
clear
clc
n=5;
cvx_status='Infeasible';
while ~strcmp(cvx_status, 'Solved')
        clear cvx_status
    c = randn(n,1);
    A = randn(n,n);
    b = randn(n,1);
    load problem5data.mat
    d = ones(n,1);
    cvx_begin 
        variable x(n)
        minimize(c'*x)
        subject to
        A*x == b;
        x<=d;
        x>=0;
    cvx_end
end
%% Problem 5 part b
%clear
clc
n=5;
cvx_status='Infeasible';
while ~strcmp(cvx_status,'Solved')
    clear cvx_status
    %c = randn(n,1);
    %A = randn(n,n);
    %b = randn(n,1);
    d = ones(n,1);

    cvx_begin
        variables u(n)  w(n) v(n)
        maximize( -b'*u -d'*w)
        subject to
        A'*u-v+w+c==0;
        v>=0;
        w>=0;
        u>=0;
    cvx_end
end
%% Problem 5 part b (second version)
%clear
clc
n=5;
cvx_status='Infeasible';
while ~strcmp(cvx_status,'Solved')
    clear cvx_status
    %c = randn(n,1);
    %A = randn(n,n);
    %b = randn(n,1);
    d = ones(n,1);
    cvx_begin
        variable u(n)
        obj=0
        for i=1:n
            obj = obj+ min(0, c(i)+A(:,i)'*u);
        end
        maximize(-b'*u+obj)
        subject to
        u>=0
    cvx_end
end

%% Problem 6 dual
%clear
clc
n=4;
m=3;
cvx_status='Infeasible';
while ~strcmp(cvx_status,'Solved')
    clear cvx_status
    b=randn(m,1);
    y=abs(randn(n,1));
    A=randn(m,n);
    load problem6data.mat
    cvx_begin
        variable z(m)
        obj=0;
        for i =1:n
            obj=obj+y(i)*exp(A(:,i)'*z);
        end
        maximize(b'*z-log(obj));
    cvx_end
end