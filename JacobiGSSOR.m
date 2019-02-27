A=[7 3 -1 2;3 8 1 -4;-1 1 4 -1; 2 -4 -1 6];
b=[-1;0;-3;1];
error=Inf;
tol=1e-4; 
x=rand(4,1);
n=size(x,1);
JacobItr=0;

plotJacobi=[];
while error>tol
    xo=x;   
    for i=1:n
        x(i)=b(i);        
        for j=1:i-1            
            x(i)=x(i)-A(i,j)*xo(j);         
        end
        for j=i+1:n
            x(i)=x(i)-A(i,j)*xo(j);
        end
        x(i)=x(i)/A(i,i);
    end    
    JacobItr=JacobItr+1;
    error=norm(xo-x);
    plotJacobi=[plotJacobi;error];
end
fprintf('Solution of the Jacobi system is : \n%f\n%f\n%f\n%f in %d iterations',x,JacobItr);
fprintf('\n')

A=[7 3 -1 2;3 8 1 -4;-1 1 4 -1; 2 -4 -1 6];
b=[-1;0;-3;1];
error=Inf;
tol=1e-4; 
x=rand(4,1);
n=size(x,1);
GaussItr=0;
plotGauss=[];
while error>tol
    xogs=x;
    for i=1:n
        x(i)=b(i);        
        for j=1:i-1            
            x(i)=x(i)-A(i,j)*x(j);         
        end
        for j=i+1:n
            x(i)=x(i)-A(i,j)*xogs(j);
        end
        x(i)=x(i)/A(i,i);
    end        
    GaussItr=GaussItr+1;
    error=norm(xogs-x);
    plotGauss=[plotGauss;error];
end

fprintf('Solution of the GS System is : \n%f\n%f\n%f\n%f in %d iterations',x,GaussItr);
fprintf('\n')

A=[7 3 -1 2;3 8 1 -4;-1 1 4 -1; 2 -4 -1 6];
b=[-1;0;-3;1];
error=Inf;
tol=1e-4; 
x=rand(4,1);
n=size(x,1);
SORItr=0;
plotSOR=[];
w=1.4;
while error>tol
    xosor=x;
    for i=1:n
        x(i)=b(i);        
        for j=1:i-1            
            x(i)=x(i)-A(i,j)*x(j);         
        end
        for j=i+1:n
            x(i)=x(i)-A(i,j)*xosor(j);
        end
        x(i)=(1-w)*xosor(i)+w*x(i)/A(i,i);
    end        
    SORItr=SORItr+1;
    error=norm(xosor-x);
    plotSOR=[plotSOR;error];
end
fprintf('Solution of the SOR System is : \n%f\n%f\n%f\n%f in %d iterations',x,SORItr);


figure
hold on
plot(1:5:GaussItr,plotGauss(1:5:GaussItr),'LineWidth',2)
plot(1:5:JacobItr,plotJacobi(1:5:JacobItr),'LineWidth',2)
plot(1:5:SORItr,plotSOR(1:5:SORItr),'LineWidth',2)
text(GaussItr,0.2,'\downarrow')
text(GaussItr,0.3,'Gauss Seidel')
text(JacobItr,0.3,'\downarrow')
text(JacobItr,0.4,'Jacobi Method')
text(SORItr,0.2,'\downarrow')
text(SORItr,0.3,'SOR')
legend('Gauss Seidel Method','Jacobi Method', 'SOR')
ylabel('Error Value')
xlabel('Number of iterations')
title('Gauss Seidel Method Vs Jacobi Method vs SOR Method')
