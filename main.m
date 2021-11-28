%% 优化方法大作业:使用三种牛顿法解决问题
% Func1.m：标准newton方法
% Func2.m：出发点不可行newton方法
% Func3.m：对偶newton方法
% 输出图形：迭代所得函数值与迭代次数的关系
%% 初始化，随机选择矩阵A、正向量x
clear
clc

n=100;                  %A矩阵的列数
p=30;                   %A矩阵的行数
A=randn(p,n);           %随机生成矩阵A
RA=rank(A);

while RA ~= p           %判断矩阵是否满秩，直到符合要求为止
    fprintf('A不是满秩矩阵，重新生成A\n');
    A=randn(p,n);
end
fprintf('A是满秩矩阵，继续程序！\n');
x=rand(n,1);           %生成向量x在[0,1]上均匀分布
b=A*x;                 %生成b
%% 初始化参数
% α=0.01，β=0.5 误差阈值ε=10^(-8)，最大迭代次数MaxIter=100
alpha=0.01;             %设置α值
beta=0.5;               %设置β值
yita=10^(-8);           %设置阈值
MaxIter=100;         %最大迭代次数
%% 标准Newton方法
fprintf('\n')
fprintf('标准Newton方法：')
fprintf('\n')
figure(1)
[figure0,calTime0]=Func1(x,MaxIter,yita,alpha,beta,A,p,n);
calltime1=1:1:calTime0;
plot(calltime1,figure0,'b*-')
title('标准Newton方法');
%% 不可行初始点Newton方法
fprintf('\n')
fprintf('不可行初始点Newton方法：')
fprintf('\n')
figure(2)
fprintf("初始点为x0=x：\n")
[figure1,calTime1]=Func2(x,MaxIter,yita,alpha,beta,A,b,p,n);
calltime2=1:1:calTime1;
plot(calltime2,figure1,'b*-')
hold on
fprintf("初始点为x0=1：\n")
[figure2,calTime2]=Func2(ones(n,1),MaxIter,yita,alpha,beta,A,b,p,n); 
calltime3=1:1:calTime2;
plot(calltime3,figure2,'ro-')
title('不可行初始点Newton方法');
legend('初始点为X0=x','初始点为x0=1');
%% 对偶Newton方法
fprintf('\n')
fprintf('对偶Newton方法：\n')
figure(3)
[ figure3,calTime3]=Func3(MaxIter,yita,alpha,beta,A,b,p);
calltime4=1:1:calTime3;
plot(calltime4,figure3,'b*-')
title('对偶Newton方法');
hold on