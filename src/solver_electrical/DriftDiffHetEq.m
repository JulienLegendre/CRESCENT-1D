function d = DriftDiffHetEq(d,tol,iterVmax)
%DRIFTDIFFEQ Compute the equilibrium quantities of component d.
    arguments
        d (1,1) activeComponent
        tol (1,1) double
        iterVmax (1,1) uint32
    end

    global eps0 kb e

    %% Initialisation and definition of parameters
    
    % Initialisation
    d.Nm=zeros(d.Nzet,1);
    d.zeq=zeros(d.Nzet-1,1);
    d.Dgv=zeros(d.Nzet,1);
    d.xinv=zeros(d.Nzet,1);
    d.xipv=zeros(d.Nzet,1);
    d.Dg=zeros(d.Nzet+d.Nle-1,1);
    d.xin=zeros(d.Nzet+d.Nle-1,1);
    d.xip=zeros(d.Nzet+d.Nle-1,1);
    d.Ncm=zeros(d.Nzet+d.Nle-1,1);
    d.Nvm=zeros(d.Nzet+d.Nle-1,1);

    d.CTm.mu_min=zeros(d.Nzet+d.Nle-1,2);
    d.CTm.DmuT=zeros(d.Nzet+d.Nle-1,2);
    d.CTm.NrefT=zeros(d.Nzet+d.Nle-1,2);
    d.CTm.alpha=zeros(d.Nzet+d.Nle-1,2);

    d.nim=zeros(d.Nzet+d.Nle-1,1);
    d.taum=zeros(d.Nzet+d.Nle-1,2);
    d.Bm=zeros(d.Nzet+d.Nle-1,1);
    d.Cm=zeros(d.Nzet+d.Nle-1,2);

    % signed doping level
    Ndops=d.DopSign(d.NleRange(1:end)).*d.Ndop(d.NleRange(1:end));

    % Definition
    % essentially, defining quantities at each point of the spatial grid
    % instead of each layer
    for i=1:d.Nle
        irange=sum(d.Nze(1:i-1))+1:sum(d.Nze(1:i));

        d.Nm(sum(d.Nze(1:i-1))+3-i-(i==1):sum(d.Nze(1:i))-i+(i==d.Nle))=repmat(Ndops(i)/d.ni(1),d.Nze(i)-2+(i==1 || i==d.Nle),1);
        d.zeq(sum(d.Nze(1:i-1)-1)+1:sum(d.Nze(1:i)-1))=sqrt(d.epsr(i)*eps0*kb*d.T/(e^2*d.ni(1)));
        d.Dgv(sum(d.Nze(1:i-1))+3-i-(i==1):sum(d.Nze(1:i))-i+(i==d.Nle))=repmat(-(d.Eg(i)-d.Eg(1))/kb/d.T+log(d.Nv(i)*d.Nc(i)/d.Nv(1)/d.Nc(1)),d.Nze(i)-2+(i==1 || i==d.Nle),1);
        d.xinv(sum(d.Nze(1:i-1))+3-i-(i==1):sum(d.Nze(1:i))-i+(i==d.Nle))=repmat((d.chi(i)-d.chi(1))/kb/d.T+log(d.Nc(i)/d.Nc(1)),d.Nze(i)-2+(i==1 || i==d.Nle),1);
        d.xipv(sum(d.Nze(1:i-1))+3-i-(i==1):sum(d.Nze(1:i))-i+(i==d.Nle))=repmat(-(d.Eg(i)-d.Eg(1)+d.chi(i)-d.chi(1))/kb/d.T+log(d.Nv(i)/d.Nv(1)),d.Nze(i)-2+(i==1 || i==d.Nle),1);
        d.Dg(irange)=repmat(-(d.Eg(i)-d.Eg(1))/kb/d.T+log(d.Nv(i)*d.Nc(i)/d.Nv(1)/d.Nc(1)),d.Nze(i),1);
        d.xin(irange)=repmat((d.chi(i)-d.chi(1))/kb/d.T+log(d.Nc(i)/d.Nc(1)),d.Nze(i),1);
        d.xip(irange)=repmat(-(d.Eg(i)-d.Eg(1)+d.chi(i)-d.chi(1))/kb/d.T+log(d.Nv(i)/d.Nv(1)),d.Nze(i),1);
        d.Ncm(irange)=d.Nc(i);
        d.Nvm(irange)=d.Nv(i);

        d.CTm.mu_min(irange,:)=repmat(d.CTmodel(i).mu_min,length(irange),1);
        d.CTm.DmuT(irange,:)=repmat(d.CTmodel(i).mu_max.*(300./d.T).^d.CTmodel(i).theta1-d.CTmodel(i).mu_min,length(irange),1);
        d.CTm.NrefT(irange,:)=repmat(d.CTmodel(i).Nref.*(d.T/300).^d.CTmodel(i).theta2,length(irange),1);
        d.CTm.alpha(irange,:)=repmat(d.CTmodel(i).alpha,length(irange),1);

        d.nim(irange)=repmat(d.ni(i),length(irange),1);
        d.taum(irange,:)=repmat([d.tauN(i) d.tauP(i)],length(irange),1);
        d.Bm(irange)=repmat(d.Brad(i),length(irange),1);
        d.Cm(irange,:)=repmat([d.Cn(i) d.Cp(i)],length(irange),1);
    end
    d.dzb=d.dz./d.zeq;

    % Parameters for finite difference method
    d.a=d.dz(2:d.Nzet-1)./d.dz(1:d.Nzet-2);
    d.a=d.a-(d.a>=1).*(abs(d.a-round(d.a))./d.a<1e-15).*(d.a-round(d.a))... % reduce error if a should be an integer
         -(d.a<1).*(abs(1./d.a-round(1./d.a)).*d.a<1e-15).*(d.a-1./max(1,round(1./d.a))); 
    d.K=2*[1./(d.a+1) -1./d.a 1./(d.a.*(d.a+1))];               % 2nd derivative, 2nd order, centred, irregular
    d.KL=[d.a.^2./(1+d.a) -(1+d.a) (1+2*d.a)./(1+d.a)];         % 1st derivative, 2nd order, left, irregular
    d.KR=[-(2+d.a)./(1+d.a) (1+d.a)./d.a -1./(d.a.*(1+d.a))];   % 1st derivative, 2nd order, right, irregular

    %% Resolution of Poisson equation

    % Poisson - First guess : zero electrical field in each layer
    d.cp.Vb=zeros(d.Nzet,1);
    for i=1:d.Nle
        rg=sum(d.Nze(1:i-1)-1)+2;
        d.cp.Vb(sum(d.Nze(1:i-1)-1)+1:sum(d.Nze(1:i)-1)+1)=repmat(eqValue(d.Nm(rg),d.xinv(rg),d.xipv(rg)),[d.Nze(i),1]);
    end

    % Poisson - Iterative process (due to linearization)
    % solution of M * dVb = Y
    iterV=0;
    err=1;
    A=[d.K(1:d.Nzet-2,1);0]; % subdiagonal terms of M
    C=[0;d.K(1:d.Nzet-2,3)]; % superdiagonal terms of M
    zint=cumsum(d.Nze(1:end-1));
    while err>tol && iterV<iterVmax
        % Normalised equilibrium charge carrier concentration
        neq=exp(d.xinv).*exp(d.cp.Vb);
        peq=exp(d.xipv).*exp(-d.cp.Vb);

        % Volumetric terms
        B=[1;d.K(1:d.Nzet-2,2)-d.dzb(1:d.Nzet-2).^2.*(neq(2:d.Nzet-1)+peq(2:d.Nzet-1));1]; % main diagonal of M
        Y=[0;-(d.K(1:d.Nzet-2,1).*d.cp.Vb(1:d.Nzet-2)+d.K(1:d.Nzet-2,2).*d.cp.Vb(2:d.Nzet-1)+d.K(1:d.Nzet-2,3).*d.cp.Vb(3:d.Nzet))+...
              d.dzb(1:d.Nzet-2).^2.*(neq(2:d.Nzet-1)-peq(2:d.Nzet-1)+d.Nm(2:d.Nzet-1));0];

        % At interfaces: electric displacement is constant. The use of
        % non-centred finite differences causes the equation to involve 5
        % different nodes. This equations must thus be combined with the
        % one obtained at the two adjacent nodes to eliminate the influence
        % of nodes k-2 and k+2.
        A(zint-1)=C(zint+1).*d.a(zint-1).*d.epsr(1:end-1)'.*(A(zint-2).*d.KL(zint-2,2)-B(zint-1).*d.KL(zint-2,1));
        B(zint)=A(zint-2).*C(zint+1).*(d.a(zint-1).*d.epsr(1:end-1)'.*d.KL(zint-2,3)-d.epsr(2:end)'.*d.KR(zint,1))-...
                C(zint-1).*C(zint+1).*d.a(zint-1).*d.epsr(1:end-1)'.*d.KL(zint-2,1)+...
                A(zint-2).*A(zint).*d.epsr(2:end)'.*d.KR(zint,3);
        C(zint)=-A(zint-2).*d.epsr(2:end)'.*(C(zint+1).*d.KR(zint,2)-B(zint+1).*d.KR(zint,3));
        Y(zint)=A(zint-2).*C(zint+1).*(-d.a(zint-1).*d.epsr(1:end-1)'.*(d.KL(zint-2,1).*d.cp.Vb(zint-2)+d.KL(zint-2,2).*d.cp.Vb(zint-1)+d.KL(zint-2,3).*d.cp.Vb(zint))+...
                                                     d.epsr(2:end)'.*(d.KR(zint,1).*d.cp.Vb(zint)+d.KR(zint,2).*d.cp.Vb(zint+1)+d.KR(zint,3).*d.cp.Vb(zint+2)))-...
                C(zint+1).*Y(zint-1).*d.a(zint-1).*d.epsr(1:end-1)'.*d.KL(zint-2,1)+...
                A(zint-2).*Y(zint+1).*d.epsr(2:end)'.*d.KR(zint,3);
        
        % Thomas algorithm
        dVb=Algo_Thomas_mex(A,B,C,Y);
        d.cp.Vb=d.cp.Vb+dVb;
        iterV=iterV+1;
        err=max(abs(dVb));
    end

    % compute additional quantities related to current point
    d.cp.U=0;    
    d.cp.J=0;
    d.cp.Efnb=ones(sum(d.Nze),1);
    d.cp.Efpb=ones(sum(d.Nze),1);
    d=d.setnp;
    d=d.setEcv;

    % save values of spatially-resolved quantities at equilibrium
    d.PointEq=d.cp;

end

function V=eqValue(Nm,gam,kap)
%EQVALUE Returns the electrostatic potential of a layer with no electrical
%field.
    D=gam+kap;
    x=exp(D)*(2/Nm)^2;
    if abs(x)<1e-15 && Nm>0
        V=kap-log(abs(Nm));
    else
        V=-gam+log(abs(Nm/2))+log(-sign(Nm)+sqrt(1+x));
    end
end

