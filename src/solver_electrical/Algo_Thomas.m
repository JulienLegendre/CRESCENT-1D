function X = Algo_Thomas(A,B,C,D)
%Thomas Algorithm : Resolution of a tridiagonal system

    %Initialization
    N=length(B);
    A=[0;A];
    C=[C;0];
    s=zeros(N,1);
    e=zeros(N,1);
    s(1)=-C(1)/B(1);
    e(1)=D(1)/B(1);

    %Down Step
    for i=2:N
        s(i)=-C(i)/(A(i)*s(i-1)+B(i));
        e(i)=(D(i)-A(i)*e(i-1))/(A(i)*s(i-1)+B(i));
    end

    X=zeros(N,1);
    X(N)=e(N);

    %Up step
    for i=N-1:-1:1
        X(i)=s(i)*X(i+1)+e(i);
    end

end

