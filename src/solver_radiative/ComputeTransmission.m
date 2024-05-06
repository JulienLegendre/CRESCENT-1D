function [T,Ts] = ComputeTransmission(a,b,w,kr,epsMat,z,zs,zl)
%COMPUTETRANSMISSION Compute the transmission coefficient between severals
%points, using Fluctuational Electrodynamics formalism
    
    arguments
        a       (1,1) uint32 % starting index for the S-matrix, between 0 and N (-)
        b       (1,1) uint32 % ending index for the S-matrix, between 0 and N (-)
        w       (1,1) double % angular frequency (rad.s^-1)
        kr      (1,1) double % parallel element of k (rad.m^-1)
        epsMat  (1,:) double % dielectric function at w for the different layers (N+1 terms) (-)
        z       (1,:) double % distance matrix (z=0 at the first interface; N terms) (m)
        zs      (:,1) double % spatial grid in layer a (m)
        zl      (:,1) double % spatial grid in layer b (m)
    end

    global c %#ok<GVMIS>
    
%% Section : Computation of the S-matrices

    S=zeros(2,2,2);
    
    % Error in input
    if a>b || a<0 || b>length(z) || length(z)~=length(epsMat)-1
        
        S(:,:,1)=nan*ones(2);
        S(:,:,2)=nan*ones(2);
        S0(:,:,1)=nan*ones(2);
        S0(:,:,2)=nan*ones(2);
        Sn(:,:,1)=nan*ones(2);
        Sn(:,:,2)=nan*ones(2);
      
    % Correct input    
    else
        
        S0(:,:,1)=eye(2);
        S0(:,:,2)=eye(2);
        S(:,:,1)=eye(2);
        S(:,:,2)=eye(2);
        n=sqrt(epsMat);
        kz=sqrt(epsMat.*(w/c).^2-kr^2);
        N=length(z);
        
        % Computes the reflexion/transmission coefficient depending on the
        % polarization
        r(:,1)=(kz(1:N).*epsMat(2:N+1) - kz(2:N+1).*epsMat(1:N))...
          ./(kz(1:N).*epsMat(2:N+1) + kz(2:N+1).*epsMat(1:N));
        t(:,1)=2*kz(1:N).*n(1:N).*n(2:N+1)...
          ./(kz(1:N).*epsMat(2:N+1) + kz(2:N+1).*epsMat(1:N));
        r(:,2)=(kz(1:N) - kz(2:N+1))./(kz(1:N) + kz(2:N+1));
        t(:,2)=2*kz(1:N)./(kz(1:N) + kz(2:N+1));
      
        E=[1 exp(1i*kz(2:N).*(z(2:N)-z(1:N-1)))];
        
        % Computation of S0=S(0,s)
        for L=1:a
            oldS=permute(S0,[1 3 2]);
            
            % Data from Francoeur (2009)
            S0(1,1,:)=oldS(1,:,1).*t(L,:)*E(L)...
                     ./(1-oldS(1,:,2).*r(L,:)*E(L)^2);
            S0(1,2,:)=(oldS(1,:,2)*E(L)^2-r(L,:))...
                     ./(1-oldS(1,:,2).*r(L,:)*E(L)^2);
            S0(2,1,:)=permute(S0(1,1,:),[1 3 2]).*oldS(2,:,2).*r(L,:)*E(L)...
                     ./t(L,:)+oldS(2,:,1);
            S0(2,2,:)=oldS(2,:,2).*(r(L,:).*permute(S0(1,2,:),[1 3 2])+1)...
                     *E(L)./t(L,:);
            
        end
        
        % Computation of S=S(s,l)
        for L=a+1:b
            oldS=permute(S,[1 3 2]);
            
            S(1,1,:)=oldS(1,:,1).*t(L,:)*E(L)...
                     ./(1-oldS(1,:,2).*r(L,:)*E(L)^2);
            S(1,2,:)=(oldS(1,:,2)*E(L)^2-r(L,:))...
                     ./(1-oldS(1,:,2).*r(L,:)*E(L)^2);
            S(2,1,:)=permute(S(1,1,:),[1 3 2]).*oldS(2,:,2).*r(L,:)*E(L)...
                     ./t(L,:)+oldS(2,:,1);
            S(2,2,:)=oldS(2,:,2).*(r(L,:).*permute(S(1,2,:),[1 3 2])+1)...
                     *E(L)./t(L,:);
            
        end
        
        % Computation of Sn=S(s,N)
        Sn=S;
        for L=b+1:N
            oldS=permute(Sn,[1 3 2]);
            
            Sn(1,1,:)=oldS(1,:,1).*t(L,:)*E(L)...
                     ./(1-oldS(1,:,2).*r(L,:)*E(L)^2);
            Sn(1,2,:)=(oldS(1,:,2)*E(L)^2-r(L,:))...
                     ./(1-oldS(1,:,2).*r(L,:)*E(L)^2);
            Sn(2,1,:)=permute(Sn(1,1,:),[1 3 2]).*oldS(2,:,2).*r(L,:)*E(L)...
                     ./t(L,:)+oldS(2,:,1);
            Sn(2,2,:)=oldS(2,:,2).*(r(L,:).*permute(Sn(1,2,:),[1 3 2])+1)...
                     *E(L)./t(L,:);
        end
        
        S0=permute(S0,[1 3 2]);
        S=permute(S,[1 3 2]);
        Sn=permute(Sn,[1 3 2]);
    end
    
%% Section 2 : Computation of the A,B,C,D coefficients

    param=1e-14;

    if a>0
        %Emitter: thin layer
        Es=ones(2,1);
        Bs=(Sn(2,:,1)./(1-S0(1,:,2).*Sn(2,:,1))).'.*Es;
        As=S0(1,:,2).'.*Bs;
        numB=Bs-S(2,:,1).'.*(As+Es);
        if max(abs(numB(:)))<param*max(abs(Bs(:))) || max(abs(numB))==0 || b==length(z)
            B=zeros(2,1);
        else
            B=numB./(S(2,:,2).');
        end       
        A=S(1,:,1).'.*(As+Es)+S(1,:,2).'.*B;
        Cs=(S0(1,:,2)./(1-S0(1,:,2).*Sn(2,:,1))).'./Es;
        Ds=Sn(2,:,1).'.*Cs;
        numD=Ds-S(2,:,1).'.*Cs;
        if max(abs(numD(:)))<param*max(abs(Ds(:))) || max(abs(numD))==0 || b==length(z)
            D=zeros(2,1);
        else
            D=numD./(S(2,:,2).');
        end
        C=S(1,:,1).'.*Cs+S(1,:,2).'.*D;
        
    else
        % Emitter: Semi-infinite layer
        numB=(Sn(2,:,1)-S(2,:,1));
        limit=(abs(numB(:))<param*abs(Sn(2,:,1).') | b==length(z));
        B=limit.*zeros(2,1)+~limit.*(numB./S(2,:,2)).';
        B(logical(limit.*isnan(B)))=zeros(sum(limit.*isnan(B)==1),1);
        A=S(1,:,1).'+S(1,:,2).'.*B;
        C=zeros(2,1);
        D=zeros(2,1);
        
    end

%% Section 3 : Use of the Weyl components, integration over zs

    if a>0
        Nzs=length(zs);
        gtotres=gtot(A,B,C,D,zs,zl,w,kr,sqrt(epsMat(a+1)),sqrt(epsMat(b+1)));
        gt=permute(sum((zs(2:Nzs)-zs(1:Nzs-1))'/2.*(gtotres(:,1:Nzs-1,:)...
                                                   +gtotres(:,2:Nzs,:)),2),[1 3 2]);
    else
        kzs=sqrt(epsMat(a+1)*(w/c)^2-kr^2);
        gtotres=zeros(2,1,1);
        gt=1/(2*imag(kzs))*gtot(A,B,C,D,0,zl,w,kr,sqrt(epsMat(a+1)),sqrt(epsMat(b+1)));
    end
    
%% Section 4 : Computation of the equivalent transmission coefficient

    K1=4*(w/c)^2*1i*imag(epsMat(a+1));
    K2=1/(4*pi^2)*kr;
    
    % Transmission coeff : real(K1*(gt(1,:)+gt(2,:)))
    % But we rather return T, as it is then simpler to compute the
    % radiative quantities
    T=K2*real(K1*(gt(1,:)+gt(2,:)))';
    Ts=permute(K2*real(K1*(gtotres(1,:,:)+gtotres(2,:,:))),[2 3 1]);
    
end