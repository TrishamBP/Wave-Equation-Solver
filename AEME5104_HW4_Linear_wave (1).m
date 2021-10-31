%Solution of the one-d wave equation u_t+c*u_x=0
% for 0<x<L, t>0
%Lax, Lax-Wendroff, Euler implicit.

clear all

% 1 implies that we use Euler implicit scheme; 0 corresponds to explicit
implicit=1;
k=logical(implicit);
c = 300;
dx = 5;
L=600;
ni=0.75;%0.45;
time=1.0;

dt = ni*dx/c;
ncells=L/dx;
nx=ncells+1;
timesteps = time/dt;

syms xx

%loop to generate grid and populate vector of initial conditions between 0 and 600

for i=1:nx
    x(i) = (i-1)*dx;
    if (x(i)<0)
        uexact(1,i)=0;
    end
    if (x(i)>50&x(i)<110)
        uexact(1,i)=100*sin(pi*(x(i)-50)/60);
    else
        uexact(1,i)=0;
    end
    u(1,i)=uexact(1,i);
end

%loop to populate subsequent lines of matrix

for iter = 1:timesteps
    %current time level is n=iter and current time=(iter-1)*dt
    n=iter;
    %next time level is n+1 and next time is time=iter*dt
    %Implement boundary conditions at x=0;
    u(n+1,1) = 0;
    uexact(n+1,1)=0;
    if (logical(implicit))
        %form the Tridiagonal matrix for the nx-2 inner grid points
        at=-ni/2;
        bt=1;
        ct=ni/2;
        e=ones(nx-2,1);
        T=spdiags([at*e,bt*e,ct*e],-1:1,nx-2,nx-2);
        %fill the rhs of the equation
        rhs(1,1)=u(n,2)-(-bt)*u(n,1);
        rhs(nx-2,1)=u(n,nx-1)-(at)*u(n,nx);
        for jj=2:nx-3
            rhs(jj,1)=u(n,jj+1);
        end
%       [L,U,P]=lu(T);
%       sol=U\(L\(P'*rhs));
        sol=T\rhs;
        for jj=1:nx-2
            u(n+1,jj+1)=sol(jj);
        end
        %spy(T)
    else
        %Explicit methods- loop through the interior grid points
        for j = 2:nx-1
            % Lax-Wendroff Method
            % u(n+1,j) = u(n,j) - 0.5*ni* (u(n,j+1)- u(n,j-1)) + 0.5*ni^2 * (u(n,j+1) - 2*u(n,j) + u(n,j-1));
            % Lax Method
             u(n+1,j) = 0.5 * (u(n,j+1)+ u(n,j-1)) - .5 *ni * (u(n,j+1) -  u(n,j-1));
        end
    end
    %Exact solution
    for j = 2:nx-1
        xx=x(j)-c*(n+1)*dt;
        if (xx<0)
            uexact(n+1,j)=0;
        end
        if (xx>50&xx<110)
            uexact(n+1,j)=100*sin(pi*(xx-50)/60);
        else
            uexact(n+1,j)=0;
        end
    end
    
    %boundary conditions at x=L
    u(n+1,nx) = u(n+1,nx-1);
    uexact(n+1,nx)=uexact(n+1,nx-1);
end

%graphics

for n = 1:timesteps+1
    %hold on
    %plot(x,u(1,1:nx));
    %plot exact and numerical
    subplot, plot(x,u(n,1:nx),'b-',x,uexact(n,1:nx),'r-',x,u(1,1:nx),'r.-'),xlabel('Time (s)'),ylabel('U (m/s)');
    
    %M(:,n) = getframe;
    %plot exact solution
    %plot(x,uexact(n,1:nx),'r-',x,u(1,1:nx),'r.-');
    %plot numerical solution
    %plot(x,u(n,1:nx),'b:',x,u(1,1:nx),'r.-');
    
    axis([0 600 0 100])
    M(:,n) = getframe;
end
