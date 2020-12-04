%% Question 5
clc
clear
H=zeros(7,7);
for i=1:6
    H(i,i)=1;
    H(i+1,i)=-1;
end
H=2*H*H';
f=zeros(7,1);
G=eye(7,7);
h=30*ones(7,1);
A=zeros(7,7);
A([1,49])=1;
b=[4 0 0 0 0 0 4]';
lb=[4 5 13 4 2 12 4]';
ub=[4 9 15 8 6 18 4]';
Optimum_points=quadprog(H,f,G,h,A,b,lb,ub);
X=4*[0:6]';
Optimum_points=fliplr(Optimum_points');
Y=Optimum_points';
fprintf('The walker will pass through the following points:\n')
T=table(X,Y)
