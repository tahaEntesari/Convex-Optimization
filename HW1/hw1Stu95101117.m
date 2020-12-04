%% Q6
clc
close all
clear
for n=2
    figure('units','normalized','outerposition',[0 0 1 1])
    q1=rand(n,n);
    q1=.5*(q1+q1')*5;
    q2=rand(n,n);
    q2=.5*(q2+q2')*5;
    %q1=q1/norm(q1);
    %q2=q2/norm(q2);
    tries=500000;
    parfor i=1:tries
        fprintf('%d\n',i);
        x=rand(n,1);
        x=(-1).^randi(2,n,1).*x/norm(x);
        %x=randn(n,1);
        %x=x/norm(x);
        S1(i)=x'*q1*x;
        S2(i)=x'*q2*x;
    end
    plot(S1,S2,'*','LineWidth',3);
    grid on
    title(['\textbf{$S = \{(X^TQ_1X,X^TQ_2X)| \mid\mid{X}\mid\mid_2 = 1\}, n = ', num2str(n),'$}'],'interpreter','latex','FontSize',15)
    xlabel('\textbf{$X^TQ_1X$}','interpreter','latex','FontSize',15)
    ylabel('\textbf{$X^TQ_2X$}','interpreter','latex','FontSize',15)
    txt='q1%d';
    txt=sprintf(txt,n);
    saveas(gca,txt,'epsc')
end



%% Q7
clc
clear
close all
tic
n=2000;
range=10;
y1=[rand(1,n),-rand(1,n),rand(1,n),-rand(1,n)].*randi(range,1,4*n);
y2=[rand(1,n),-rand(1,n),-rand(1,n),rand(1,n)].*randi(range,1,4*n);
x1=[rand(1,n),-rand(1,n),rand(1,n),-rand(1,n)].*randi(range,1,4*n);
x2=[rand(1,n),-rand(1,n),-rand(1,n),rand(1,n)].*randi(range,1,4*n);
a=2*x1.^4+x2.^4;
b=x1.*x2;
%
res=zeros(4*n,1);
for i=1:4*n
    
    fprintf('%d\n',i);
    res(i)=all(a+b.*(y1(i)*x2.^2+y2(i)*x1.^2)>=0);
end
res=res==1;
toc
figure('units','normalized','outerposition',[0 0 1 1])
plot(y1(res'),y2(res'),'*')
title('\textbf{$S = \{(y_1, y_2)| 2x_1^4+x_2^4+y_1x_1x_2^3+y_2x_1^3x_2\geq 0, \forall (x_1,x_2)\in R \}$}','interpreter','latex','FontSize',15)
xlabel('\textbf{$y_1$}','interpreter','latex','FontSize',15)
ylabel('\textbf{$y_2$}','interpreter','latex','FontSize',15)
grid on
saveas(gca,'q7','epsc');

%% Q8
clear
clc
close all
A=imread('HajiFirouz.jpg');
A=im2double(A) ;
A=rgb2gray(A) ;
[u,sigma,v]=svd(A);
w=1;
for k=[1:5:40,50,100,200,400]
    uk=u(:,1:k);
    vk=v(:,1:k);
    sigmak=sigma(1:k,1:k);
    Ak=uk*sigmak*vk';
    figure();
    subplot(2,1,1)
    imshow(A)
    txt='The value of the cost function for k = %d is %f\n';
    txt=sprintf(txt,k,norm(A-Ak));
    title(txt)
    subplot(2,1,2)
    imshow(Ak)
    title(['Compressed image for k = ', num2str(k)], 'FontSize',15)
    txt='ReconstructedHajik=%d.jpg';
    txt=sprintf(txt,k);
    imwrite(Ak,txt);
    err(w)=norm(Ak-A);
    w=w+1;
end