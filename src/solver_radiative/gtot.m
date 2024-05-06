function gt = gtot(A,B,C,D,zs,zl,w,kr,ns,nl)
%GTOT Computes the product/sum of the Weyl components of the Green tensors
    arguments
        A   (2,1) double % coefficients found using the S-matrix (A(1):p-polarized; A(2):s-polarized). A: forward emitted - forward received
        B   (2,1) double % B: forward emitted - backward received
        C   (2,1) double % C: backward emitted - forward received
        D   (2,1) double % D: backward emitted - backward received
        zs  (:,1) double % spatial grid for emission (m)
        zl  (:,1) double % spatial grid for reception (m)
        w   (1,1) double % angular frequency (rad.s^-1)
        kr  (1,1) double % parallel component of the angular wavevector (rad.m^-1)
        ns  (1,1) double % complex refractive index of the emitter (-)
        nl  (1,1) double % complex refractive index of the absorber (-)
    end

    global c %#ok<GVMIS>
    ks=w*ns/c;
    kl=w*nl/c;
    kzs=sqrt(ks^2-kr^2);
    kzl=sqrt(kl^2-kr^2);
    Nzs=length(zs);
    Nzl=length(zl);

    % Don't do the calculation if all coefficients are zero
    gt=zeros(Nzs,Nzl,2);
    if all(A==0) && all(B==0) && all(C==0) && all(D==0)
        gt=permute(gt,[3 1 2]);
        return
    end

    % Definition of the exponential terms
    % Being used several times, defining this term here allows to reduce
    % computational time
    rm1=exp(log(A(1))/2+1i*(-kzs*zs))*exp(log(A(1))/2+1i*(+kzl*zl'));
    rm2=exp(log(B(1))/2+1i*(-kzs*zs))*exp(log(B(1))/2+1i*(-kzl*zl'));
    rm3=exp(log(C(1))/2+1i*(+kzs*zs))*exp(log(C(1))/2+1i*(+kzl*zl'));
    rm4=exp(log(D(1))/2+1i*(+kzs*zs))*exp(log(D(1))/2+1i*(-kzl*zl'));
    rt1=exp(log(A(2))/2+1i*(-kzs*zs))*exp(log(A(2))/2+1i*(+kzl*zl'));
    rt2=exp(log(B(2))/2+1i*(-kzs*zs))*exp(log(B(2))/2+1i*(-kzl*zl'));
    rt3=exp(log(C(2))/2+1i*(+kzs*zs))*exp(log(C(2))/2+1i*(+kzl*zl'));
    rt4=exp(log(D(2))/2+1i*(+kzs*zs))*exp(log(D(2))/2+1i*(-kzl*zl'));
    
    % Computation of the Weyl components of Green's tensors
    % (e/h: electric/magnetic field, r/t(heta)/z: direction)
    gerr=1i*kzl/(2*ks*kl)*          (+rm1-rm2-rm3+rm4);
    ghtr=conj(kl/(2*ks)*            (-rm1-rm2+rm3+rm4));
    gerz=1i*kzl*kr/(2*kzs*ks*kl)*   (-rm1+rm2-rm3+rm4);
    ghtz=conj(kl*kr/(2*ks*kzs)*     (+rm1+rm2+rm3+rm4));
    gett=1i/(2*kzs)*                (+rt1+rt2+rt3+rt4);
    ghrt=conj(kzl/(2*kzs)*          (+rt1-rt2+rt3-rt4));
    
    % Compute the total TE and TM terms
    gt(:,:,1)=gerr.*ghtr+gerz.*ghtz;
    gt(:,:,2)=-gett.*ghrt;
    gt=permute(gt,[3 1 2]);

end


