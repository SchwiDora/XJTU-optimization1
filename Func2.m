function[output,IterK]=Func2(xk,MaxTime,err,alpha,beta,A,b,p,n)
dk=randn(p,1);           %随机生成矩阵v0
for IterK=1:MaxTime
    Grad=log(xk)+1;                         %计算函数的梯度
    Hessian=diag(1./xk);                    %计算函数的Hessian矩阵
    r=[Grad+A'*dk;A*(xk)-b];                %原对偶残差
    output(IterK)=xk'*log(xk);              %输出目标函数值
    Rans=-[Hessian,A';A,zeros(p,p)]\r;      %求解牛顿方向和对偶变量
    xnt=Rans(1:n);                          %求解Newton方向Dnt
    dnt=Rans(n+1:n+p);                      %求解对偶变量   
    
    %停止准则：
    if norm(r)<=err
        break;
    end
    
    %对原多残差||r||2进行回溯，确定步长tk
    t=1;
    while (min(xk+t*xnt)<=0)                %使x在定义域内 
        t=beta*t;
    end
    
    %回溯直线搜索停止准则
    while norm([log(xk+t*xnt)+1+A'*(dk+t*dnt);A*(xk+t*xnt)-b])>(1-alpha*t)*norm(r)
        t=beta*t;
    end
    xk=xk+t*xnt;                       
    dk=dk+t*dnt;                       
   
end