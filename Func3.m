function[output,IterK]=Func3(MaxTime,err,alpha,beta,A,b,p)
v=randn(p,1); 
for IterK=1:MaxTime
    Grad=b-A*exp(-A'*v-1);                           %计算梯度
    Hesssian=A*diag(exp(-A'*v-1))*A';                %计算Hessian矩阵
    vnt=-Hesssian\Grad;                              %计算牛顿方向dnt
    Lamd_square=Grad'*(Hesssian^-1)*Grad;                  %计算牛顿减小量λ   
    output(IterK)=b'*v+sum(exp(-A'*v-1));            %计算目标函数值最优值
    
    %停止准则：
    if Lamd_square<=2*err
        break;
    end
    
    %回溯直线搜索
    t=1;
    while b'*(v+t*vnt)+sum(exp(-A'*(v+t*vnt)-1))>=b'*v+sum(exp(-A'*v-1))+alpha*t*Grad'*vnt
        t=beta*t;
    end
    v=v+t*vnt;
end