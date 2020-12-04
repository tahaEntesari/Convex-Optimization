%% Question 6
clear
clc
close all
run('C:\Users\marson\Dropbox\Convex Optimization\HW2\treatment_planning_data.m')
% coeffs as "a" and bounds as "c"
a1=[-eye(400), Aother];
c1=Dother*ones(400,1);
a2=[zeros(100,400), -Atumor];
c2=-Dtarget*ones(100,1);
a3=[-eye(400), zeros(400,300)];
c3=zeros(400,1);
a4=[zeros(300,400), eye(300)];
c4=Bmax*ones(300,1);
a5=[zeros(300,400), -eye(300)];
c5=zeros(300,1);
A=[a1;a2;a3;a4;a5];
c=[c1;c2;c3;c4;c5];

f=[ones(n,1);zeros(mother,1)];
x=linprog(f,A,c);
tumordose=Atumor*x(401:end);
otherdose=Aother*x(401:end);
figure('units','normalized','outerposition',[0 0 1 1])
hist(tumordose)
title('\textbf{Dose histogram for voxels containing the tumor}','interpreter','latex','FontSize',25);
xlabel('\textbf{Dose}','interpreter','latex','FontSize',25);
grid on
saveas(gca,'tumordose','epsc');
figure('units','normalized','outerposition',[0 0 1 1])
hist(otherdose)
title('\textbf{Dose for voxels \textbf{not} containing the tumor}','interpreter','latex','FontSize',25);
xlabel('\textbf{Dose}','interpreter','latex','FontSize',25);
grid on
saveas(gca,'nontumordose','epsc');