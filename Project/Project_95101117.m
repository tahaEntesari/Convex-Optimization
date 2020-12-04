%% Step 1: Choosing a Portfolio
clc
clear
pbar=load("C:\Users\marson\Dropbox\Convex Optimization\Project\data\pbar.txt");
sigma=load("C:\Users\marson\Dropbox\Convex Optimization\Project\data\sigma.txt");
m=30;
eta=logspace(-2,3,m);
n=length(pbar);
% Solving for each eta
for i=1:m
    cvx_begin
    variable x(n)
    maximize(pbar'*x-eta(i)*x'*sigma*x)
    subject to
    sum(x)==1;
    x>=0;
    cvx_end
    optimumValue(i)=cvx_optval;
    X(i).x=x;
    standardDeviation(i)=sqrt(x'*sigma*x);
    meanReturn(i)=pbar'*x;
end
% Plotting the figure('units','normalized','outerposition',[0 0 1 1])s
% meanReturn Vs SD
figure('units','normalized','outerposition',[0 0 1 1])
plot(standardDeviation,meanReturn,'LineWidth',3);
title('\textbf{Mean return - Risk }','interpreter','latex','FontSize',15);
xlabel('\textbf{$ x^T\Sigma x$}','interpreter','latex','FontSize',15);
ylabel('\textbf{$ \bar{p}^Tx$}','interpreter','latex','FontSize',15);
grid on
% meanReturn Vs eta times SD on semilogx
figure('units','normalized','outerposition',[0 0 1 1])
semilogx(eta.*standardDeviation,meanReturn,'LineWidth',3);
title('\textbf{Mean return - $\eta$ * Risk}','interpreter','latex','FontSize',15);
xlabel('\textbf{$ \eta x^T\Sigma x$}','interpreter','latex','FontSize',15);
ylabel('\textbf{$ \bar{p}^Tx$}','interpreter','latex','FontSize',15);
grid on
% Semilogx of portfolio allocation vs eta
figure('units','normalized','outerposition',[0 0 1 1])
c=[X(:).x];
for i=1:n
    semilogx(eta,c(i,:),'LineWidth',3);
    hold on
end
grid on
title('\textbf{Allocation of each stock subject to tradeoff coefficient}','interpreter','latex','FontSize',15);
xlabel('\textbf{$\mathbf{\eta}$}','interpreter','latex','FontSize',15);
ylabel('\textbf{Amount of each stock}','interpreter','latex','FontSize',15);
save('step1Data.mat','c','eta','standardDeviation','meanReturn');


%% Step 2: Short Positions
clear
clc
close all
pbarShort=load("C:\Users\marson\Dropbox\Convex Optimization\Project\data\pbar_short.txt");
sigmaShort=load("C:\Users\marson\Dropbox\Convex Optimization\Project\data\sigma_short.txt");
m=30;
eta=logspace(-2,3,m);
n=length(pbarShort);
gamma=.5;
% Optimization with short selling
for i=1:m
    cvx_begin
    variables xs(n) xl(n)
    maximize(pbarShort'*(xl-xs)-eta(i)*(xl-xs)'*sigmaShort*(xl-xs))
    subject to
    sum(xl)==1;
    xl>=0;
    xs>=0;
    sum(xs)<=gamma*sum(xl);
    cvx_end
    optimumValue(i)=cvx_optval;
    XShort(i).x=(xl-xs);
    standardDeviationShort(i)=sqrt((xl-xs)'*sigmaShort*(xl-xs));
    meanReturnShort(i)=pbarShort'*(xl-xs);
end
% Optimization problem without short selling
for i=1:m
    cvx_begin
    variable x(n)
    maximize(pbarShort'*x-eta(i)*x'*sigmaShort*x)
    subject to
    sum(x)==1;
    x>=0;
    cvx_end
    optimumValue(i)=cvx_optval;
    X(i).x=x;
    standardDeviation(i)=sqrt(x'*sigmaShort*x);
    meanReturn(i)=pbarShort'*x;
end
% Plotting the figures
figure('units','normalized','outerposition',[0 0 1 1])
plot(standardDeviation,meanReturn,'LineWidth',3);
hold on
plot(standardDeviationShort,meanReturnShort,'LineWidth',3);
title('\textbf{Mean return - Risk with new added stock}','interpreter','latex','FontSize',15);
xlabel('\textbf{$ x^T\Sigma x$}','interpreter','latex','FontSize',15);
ylabel('\textbf{$ \bar{p}^Tx$}','interpreter','latex','FontSize',15);
grid on
legend('Without short positions','With short positions','interpreter','latex','FontSize',15,'Location','NorthWest')

figure('units','normalized','outerposition',[0 0 1 1])
semilogx(standardDeviation.*eta,meanReturn,'LineWidth',3);
hold on
semilogx(eta.*standardDeviationShort,meanReturnShort,'LineWidth',3);
title('\textbf{Mean return - Risk with new added stock}','interpreter','latex','FontSize',15);
xlabel('\textbf{$\eta x^T\Sigma x$}','interpreter','latex','FontSize',15);
ylabel('\textbf{$ \bar{p}^Tx$}','interpreter','latex','FontSize',15);
grid on
legend('Without short positions','With short positions','interpreter','latex','FontSize',15,'Location','NorthWest')


%% Step 3: Portfolio Optimization with Real Data
clear
clc
close all
load("C:\Users\marson\Dropbox\Convex Optimization\Project\data\stock_dataset.mat");
pbarestimate=mean(past_returns,1)';
sigmaestimate=cov(past_returns);
m=30;
eta=logspace(-2,3,m);
n=length(pbarestimate);
gamma=.5;
% Optimization with short selling
for i=1:m
    cvx_begin
    variables xs(n) xl(n)
    maximize(pbarestimate'*(xl-xs)-eta(i)*(xl-xs)'*sigmaestimate*(xl-xs))
    subject to
    sum(xl)==1;
    xl>=0;
    xs>=0;
    sum(xs)<=gamma*sum(xl);
    cvx_end
    optimumValue(i)=cvx_optval;
    XShort(i).x=(xl-xs);
    standardDeviationShort(i)=sqrt((xl-xs)'*sigmaestimate*(xl-xs));
    meanReturnShort(i)=pbarestimate'*(xl-xs);
end
% Optimization problem without short selling
for i=1:m
    cvx_begin
    variable x(n)
    maximize(pbarestimate'*x-eta(i)*x'*sigmaestimate*x)
    subject to
    sum(x)==1;
    x>=0;
    cvx_end
    optimumValue(i)=cvx_optval;
    X(i).x=x;
    standardDeviation(i)=sqrt(x'*sigmaestimate*x);
    meanReturn(i)=pbarestimate'*x;
end
figure('units','normalized','outerposition',[0 0 1 1])
plot(standardDeviation,meanReturn,'LineWidth',3);
hold on
plot(standardDeviationShort,meanReturnShort,'LineWidth',3);
title('\textbf{Mean return - Risk based on stock data}','interpreter','latex','FontSize',15);
xlabel('\textbf{$ x^T\Sigma x$}','interpreter','latex','FontSize',15);
ylabel('\textbf{$ \bar{p}^Tx$}','interpreter','latex','FontSize',15);
grid on
legend('Without short positions','With short positions','interpreter','latex','FontSize',10,'Location','SouthEast');
% Part d: returns
c=[X(:).x];
cShort=[XShort(:).x];
for i=1:m
    noShortReturns(:,:,i)=future_returns.*c(:,i)';
    withShortReturns(:,:,i)=future_returns.*cShort(:,i)';
end
cumNoShort=cumsum(noShortReturns,1);
cumWithShort=cumsum(withShortReturns,1);
for i=1:n
    figure('units','normalized','outerposition',[0 0 1 1])
    semilogx(eta,squeeze(cumNoShort(end,i,:)),'LineWidth',3);
    hold on
    semilogx(eta,squeeze(cumWithShort(end,i,:)),'LineWidth',3);
    legend('Without short selling','With short selling','interpreter','latex','FontSize',15,'Location','NorthEast');
    grid on
    txt='Portfolio number %d';
    txt2='portfolio%d';
    txt=sprintf(txt,i);
    txt2=sprintf(txt2,i);
    txt=['\textbf{', txt, '}'];
    title(txt,'interpreter','latex','FontSize',15);
    xlabel('\textbf{$\eta$}','interpreter','latex','FontSize',15);
    ylabel('\textbf{Portfolio return}','interpreter','latex','FontSize',15);
    saveas(gca,txt2,'epsc');
end
% Total return plot
figure('units','normalized','outerposition',[0 0 1 1])
totalReturn=cumsum(future_returns*c,1);
totalReturnShort=cumsum(future_returns*cShort,1);
semilogx(eta,totalReturn(end,:),'LineWidth',3)
hold on
semilogx(eta,totalReturnShort(end,:),'LineWidth',3)
grid on
legend('Without short selling','With short selling','interpreter','latex','FontSize',15);
title('\textbf{Total return of all portfolios with respect to $ \eta$}','interpreter','latex','FontSize',15);
xlabel('\textbf{$\eta$}','interpreter','latex','FontSize',15);
ylabel('\textbf{Total return}','interpreter','latex','FontSize',15);


%% Step 4:  Improve your estimate
clear
clc
close all
load("C:\Users\marson\Dropbox\Convex Optimization\Project\data\stock_dataset.mat");
pbarestimate=mean(past_returns,1)';
S=cov(past_returns);
n=12;
e=ones(12,1);
lambda=10^(-5);
cvx_begin
    variable theta(n,n)
    obj=0;
    for i=1:n
        for j=1:n
            if j==i 
                continue
            end
            obj=obj+abs(theta(i,j));
        end
    end
    %maximize(log_det(theta)-trace(S*theta)-lambda*(e'*theta*e-trace(theta)))
    maximize(log_det(theta)-trace(S*theta)-lambda*obj)
    subject to
    theta>=0;
cvx_end
sigmaestimate=inv(theta);
%
m=30;
eta=logspace(-2,3,m);
n=length(pbarestimate);
gamma=.5;
for i=1:m
    cvx_begin
    variable x(n)
    maximize(pbarestimate'*x-eta(i)*x'*sigmaestimate*x)
    subject to
    sum(x)==1;
    x>=0;
    cvx_end
    optimumValue(i)=cvx_optval;
    X(i).x=x;
    standardDeviation(i)=sqrt(x'*sigmaestimate*x);
    meanReturn(i)=pbarestimate'*x;
end
for i=1:m
    cvx_begin
    variable x(n)
    maximize(pbarestimate'*x-eta(i)*x'*S*x)
    subject to
    sum(x)==1;
    x>=0;
    cvx_end
    optimumValue(i)=cvx_optval;
    XSimple(i).x=x;
    standardDeviationSimple(i)=sqrt(x'*S*x);
    meanReturnSimple(i)=pbarestimate'*x;
end
c=[X(:).x];
cSimple=[XSimple(:).x];
for i=1:m
    noShortReturns(:,:,i)=future_returns.*c(:,i)';
    noShortReturnsSimple(:,:,i)=future_returns.*cSimple(:,i)';
end
cumNoShort=cumsum(noShortReturns,1);
cumNoShortSimple=cumsum(noShortReturnsSimple,1);
for i=1:n
    figure('units','normalized','outerposition',[0 0 1 1])
    semilogx(eta,squeeze(cumNoShort(end,i,:)),'LineWidth',3);
    hold on
    semilogx(eta,squeeze(cumNoShortSimple(end,i,:)),'LineWidth',3);
    legend('Using new estimate','Simple previous estimate','interpreter','latex','FontSize',15,'Location','West');
    grid on
    txt='Portfolio number %d';
    txt=sprintf(txt,i);
    txt2='step4_portfolio%d';
    txt2=sprintf(txt2,i);
    txt=['\textbf{', txt, '}'];
    title(txt,'interpreter','latex','FontSize',15);
    xlabel('\textbf{$\eta$}','interpreter','latex','FontSize',15);
    ylabel('\textbf{Portfolio return}','interpreter','latex','FontSize',15);
    saveas(gca,txt2,'epsc');
end
% Total return plot
figure('units','normalized','outerposition',[0 0 1 1])
totalReturn=cumsum(future_returns*c,1);
totalReturnSimple=cumsum(future_returns*cSimple,1);
semilogx(eta,totalReturn(end,:),'LineWidth',3)
hold on
semilogx(eta,totalReturnSimple(end,:),'LineWidth',3)
grid on
legend('Using new estimate','Simple previous estimate','interpreter','latex','FontSize',15);
title('\textbf{Total return of all portfolios with respect to $ \eta$}','interpreter','latex','FontSize',15);
xlabel('\textbf{$\eta$}','interpreter','latex','FontSize',15);
ylabel('\textbf{Total return}','interpreter','latex','FontSize',15);


%% Step 5: Portfolio optimization with loss risk constraints

clear
clc
close all
load("C:\Users\marson\Dropbox\Convex Optimization\Project\data\stock_dataset.mat");
pbarestimate=mean(past_returns,1)';
S=cov(past_returns);
n=12;
e=ones(12,1);
lambda=10^(-5);
cvx_begin
    variable theta(n,n)
    obj=0;
    for i=1:n
        for j=1:n
            if j==i 
                continue
            end
            obj=obj+abs(theta(i,j));
        end
    end
    %maximize(log_det(theta)-trace(S*theta)-lambda*(e'*theta*e-trace(theta)))
    maximize(log_det(theta)-trace(S*theta)-lambda*obj)
    subject to
    theta>=0;
cvx_end
sigmaestimate=inv(theta);
m=30;
eta=logspace(-2,3,m);
n=length(pbarestimate);
gamma=.5;
q=20;
beta=logspace(-4,-1,q);
squareRootSigma=sqrtm(sigmaestimate);
%
for j=1:q
    for i=1:m
        cvx_begin
            variable x(n)
            maximize(pbarestimate'*x-eta(i)*x'*sigmaestimate*x)
            subject to
            sum(x)==1;
            x>=0;
            pbarestimate'*x+norminv(beta(j))*norm(squareRootSigma*x)>-.5;
        cvx_end
        optimumValue(i)=cvx_optval;
        X(i,j).x=x;
        standardDeviation(i,j)=sqrt(x'*sigmaestimate*x);
        meanReturn(i,j)=pbarestimate'*x;
    end
end
%
n1=1;
n2=10;
n3=20;
n4=30;
c1=[X(n1,:).x];
c2=[X(n2,:).x];
c3=[X(n3,:).x];
c4=[X(n4,:).x];
for i=1:n
    txt2='Step5_portfolio%d';
    txt2=sprintf(txt2,i);
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1)
    plot(beta,c1(i,:),'LineWidth',3);
    txt=['\textbf{Allocation of portfolio ', num2str(i), ' for $ \eta = $ ', num2str(eta(n1)),' }'];
    title(txt,'interpreter','latex','FontSize',15);
    xlabel('\textbf{$\beta$}','interpreter','latex','FontSize',15);
    ylabel('\textbf{Stock amount}','interpreter','latex','FontSize',15);
    grid on
    subplot(2,2,2)
    plot(beta,c2(i,:),'LineWidth',3);
    txt=['\textbf{Allocation of portfolio ', num2str(i), ' for $ \eta = $ ', num2str(eta(n2)),' }'];
    title(txt,'interpreter','latex','FontSize',15);
    xlabel('\textbf{$\beta$}','interpreter','latex','FontSize',15);
    ylabel('\textbf{Stock amount}','interpreter','latex','FontSize',15);
    grid on
    subplot(2,2,3)
    plot(beta,c3(i,:),'LineWidth',3);
    txt=['\textbf{Allocation of portfolio ', num2str(i), ' for $ \eta = $ ', num2str(eta(n3)),' }'];
    title(txt,'interpreter','latex','FontSize',15);
    xlabel('\textbf{$\beta$}','interpreter','latex','FontSize',15);
    ylabel('\textbf{Stock amount}','interpreter','latex','FontSize',15);
    grid on
    subplot(2,2,4)
    plot(beta,c4(i,:),'LineWidth',3);
    txt=['\textbf{Allocation of portfolio ', num2str(i), ' for $ \eta = $ ', num2str(eta(n4)),' }'];
    title(txt,'interpreter','latex','FontSize',15);
    xlabel('\textbf{$\beta$}','interpreter','latex','FontSize',15);
    ylabel('\textbf{Stock amount}','interpreter','latex','FontSize',15);
    grid on
    saveas(gca,txt2,'epsc');
end
%%
return1=cumsum(future_returns*c1,1);
return2=cumsum(future_returns*c2,1);
return3=cumsum(future_returns*c3,1);
return4=cumsum(future_returns*c4,1);
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
semilogx(beta,return1(end,:),'LineWidth',3)
title(['\textbf{Cumulative total return with respect to $ \beta$ for $ \eta = ',...
    num2str(eta(n1)),' $}'],'interpreter','latex','FontSize',15);
xlabel('\textbf{$\beta$}','interpreter','latex','FontSize',15);
ylabel('\textbf{Cumulative total return}','interpreter','latex','FontSize',15);
grid on
subplot(2,2,2)
semilogx(beta,return2(end,:),'LineWidth',3)
title(['\textbf{Cumulative total return with respect to $ \beta$ for $ \eta = ',...
    num2str(eta(n2)),' $}'],'interpreter','latex','FontSize',15);
xlabel('\textbf{$\beta$}','interpreter','latex','FontSize',15);
ylabel('\textbf{Cumulative total return}','interpreter','latex','FontSize',15);
grid on
subplot(2,2,3)
semilogx(beta,return3(end,:),'LineWidth',3)
title(['\textbf{Cumulative total return with respect to $ \beta$ for $ \eta = ',...
    num2str(eta(n3)),' $}'],'interpreter','latex','FontSize',15);
xlabel('\textbf{$\beta$}','interpreter','latex','FontSize',15);
ylabel('\textbf{Cumulative total return}','interpreter','latex','FontSize',15);
grid on
subplot(2,2,4)
semilogx(beta,return4(end,:),'LineWidth',3)
grid on
title(['\textbf{Cumulative total return with respect to $ \beta$ for $ \eta = ',...
    num2str(eta(n4)),' $}'],'interpreter','latex','FontSize',15);
xlabel('\textbf{$\beta$}','interpreter','latex','FontSize',15);
ylabel('\textbf{Cumulative total return}','interpreter','latex','FontSize',15);
grid on
%%
fr=cumsum(future_returns,1);
fr=fr(end,:);
c=zeros(n,m,q);
for i=1:n
    for j=1:n
        for k=1:q
            c(i,j,k)=X(j,k).x(i);
        end
    end
end
for j=1:m
    for k=1:q
        Returns(j,k)=fr*squeeze(c(:,j,k));
    end
end
figure('units','normalized','outerposition',[0 0 1 1])
[a,b]=meshgrid(beta,eta);
surf(a,b,Returns);
[etamax, betamax]=find(Returns==max(max(Returns)),1);
x=X(etamax,betamax).x;
etamax=eta(etamax);
betamax=beta(betamax);
save('portfolio.mat','x','etamax','betamax','sigmaestimate','pbarestimate');


%% Step 6: maximum risk for the chosen portfolio
%clear
close all
clc
%load portfolio.mat
n=12;

cvx_begin
    variable S(n,n)
    maximize(x'*S*x)
    subject to
    S>=0;
    diag(S)==diag(sigmaestimate);
    for i=1:12
        for j=1:12
            if i==j
                continue
            elseif sigmaestimate(i,j)>0
                S(i,j)>=0;
            else
                S(i,j)<=0;
            end
        end
    end
    
cvx_end


