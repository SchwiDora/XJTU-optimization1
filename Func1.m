function[output,IterK]=Func1(xk,MaxTime,err,alpha,beta,A,p,n)
for IterK=1:MaxTime
    Grad=log(xk)+1;                                   %计算梯度
    Hessian=diag(1./xk);                              %计算Hessian矩阵
    Rans=[Hessian,A';A,zeros(p,p)]\[-Grad;zeros(p,1)];%计算Dnt方向和对偶变量
    Xnt=Rans(1:n);                                    %计算Dnt\牛顿方向
    Lamb_squre=(Xnt'*Hessian*Xnt);                         %计算牛顿减小量的平方
    output(IterK)=xk'*log(xk);                        %计算目标函数的最优值
    
    %停止准则：
    if Lamb_squre<=2*err
        break;
    end
    
    %回溯直线搜索
    t=1;
    while (min(xk+t*Xnt)<=0)       
        t=beta*t;
    end
    while (xk+t*Xnt)'*log(xk+t*Xnt)>=(xk)'*log(xk)+alpha*t*Grad'*Xnt%回溯
        t=beta*t;
    end
    xk=xk+t*Xnt;                     
end