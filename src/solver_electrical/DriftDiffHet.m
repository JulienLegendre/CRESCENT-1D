function [d,calc_state,varargout] = DriftDiffHet(d,Vfb,G,id,options)
%DRIFTDIFFHET Computes the state of heterojunctions under non-equilibrium 
%conditions. Can account for thermionic emission and intraband tunnelling.

    arguments
        d (1,1) activeComponent % current component of interest
        Vfb (1,1) double        % applied voltage (V)
        G (:,1) double          % radiative generation rate (m^-2.s^-1)
        id (1,1) uint64         % location in the voltage grid (-)     
        options (1,1) struct    % options from crescent1D
    end

    global e kb hb %#ok<GVMIS>

    calc_state="over - convergence"; % stays so until the end only if it does not diverge during the resolution
    weightFermiChange="both";
    iterE=0;
    err=10; % arbitrarily large error, must be larger than tolerance to enter the while loop
    dOld=d; % save component before starting calculation
    mu=zeros(d.Nzet+d.Nle-1,2);
    idint=cumsum(d.Nze(1:end)); % index of interfaces

    % compute the chemical potential used for computing the radiative terms
    % it will be compared with the updated value mumean
    murad=zeros(d.Nzet+d.Nle-1,1);
    for k=1:d.Nle
        murad(sum(d.Nze(1:k-1))+1:sum(d.Nze(1:k)))=d.mu(d.NleRange(k),id);
    end
    mumean=murad;

    Vc=d.cp.Vb;
    Gem=d.cp.GaEm; % radiative recombination rate

    while err>options.tol*d.wnext && iterE<options.iterEmax        
        %% Computation of required quantities : mobility, recombinations, susceptibilities
        
        % Setting the previous iteration results
        EfnbOld=d.cp.Efnb;
        EfpbOld=d.cp.Efpb;
        nOld=d.cp.n;
        pOld=d.cp.p;
        Vb=d.z2za(d.cp.Vb); % changing mesh to account for n/p discontinuities

        iterEf=0;
        iterEfd=Inf;
        if options.iterEf % authorise to iterate on quasi Fermi levels
            iterEfmax=5e3;
        else % forbid to iterate on quasi Fermi levels
            iterEfmax=2;
        end
        clear errm2n errm2p
        errm2=1;errm2n=1;errm2p=1;
        tol2=1e-4;

        % iterations on drift-diffusion and continuity equations stop when
        %   - the solution found after the first iteration (without weight) is physical
        %   - the number of iterations exceeds the maximum authorised
        %   - the errors on quasi Fermi levels become smaller that the provided tolerance
        while (~d.physEf(id) || iterEf==0) && iterEf<iterEfmax && (errm2n(end)>=tol2*d.wefn || errm2p(end)>=tol2*d.wefp)
            % mobility, diffusion coefficient and recombination rate
            if all(d.cp.n>0) && all(d.cp.p>0)
                mu=(d.CTm.mu_min+d.CTm.DmuT./(1+(repmat(d.cp.n+d.cp.p,1,2)/1e6./d.CTm.NrefT).^d.CTm.alpha))/1e4;
            end
            DtoEq=(d.cp.n.*d.cp.p-d.nim.^2).*(abs(d.cp.n.*d.cp.p./d.nim.^2-1)>=5e-15);
            d.cp.R=repmat(DtoEq,1,3).*[1./(d.taum(:,2).*(d.cp.n+d.nim)+d.taum(:,1).*(d.cp.p+d.nim)) d.Bm d.Cm(:,1).*d.cp.n+d.Cp(:,2).*d.cp.p];
            d.cp.D=kb*d.T*mu/e;
            R=sum(d.cp.R,2);
            d.cp.Rs=d.locSurfRec*prod(d.Sr).*DtoEq([1;end])./(d.Sr(1)*(d.cp.n([1;end])+d.ni([1;end])')+d.Sr(2)*(d.cp.p([1;end])+d.ni([1;end])'));
            
            % Slight correction of the emission rate with the new 
            % electrochemical potential, using Boltzmann approximation.
            % Sometimes, the significant change of Geq at each iteration 
            % prevent reaching convergence. Then, a factor lower than 1
            % can be put in front of (mumean-murad) to reduce the change in 
            % Geq between iteration. Since mumean must converge towards 
            % murad (because of the iterative process on mu in the main 
            % code), the fact that the correction does not follow the
            % Bose-Einstein distribution is not of concern.
            if isempty(Gem) || ~options.IterChemPotential
                Geq=G;
            else  
                Geq=G+Gem.*(1-exp((mumean-murad)*e/kb/d.T)); 
            end
            RG=e*(R-Geq);

            EfnbOld2=d.cp.Efnb;
            EfpbOld2=d.cp.Efpb;

            % conduction/valence band
            d=d.setEcv;
            
            % Definition of the susceptibility terms, relating the current
            % to the difference of Slotboom variable Efb
            Gn=e*d.ni(1)./d.dza.*exp(d.xin(2:end)+Vb(2:end)).*(d.cp.D(1:end-1,1)+d.cp.D(2:end,1))/2.*Ber(Vb(2:end)-Vb(1:end-1));
            Gnth=d.Ar(1:end-1,1)*d.T^2.*min((exp([-d.cp.Ec(idint(1:end-1)) -d.cp.Ec(idint(1:end-1)+1)]*e/kb/d.T)),[],2);
            Gp=-e*d.ni(1)./d.dza.*exp(d.xip(1:end-1)-Vb(1:end-1)).*(d.cp.D(1:end-1,2)+d.cp.D(2:end,2))/2.*Ber(Vb(2:end)-Vb(1:end-1));
            Gpth=-d.Ar(1:end-1,2)*d.T^2.*min((exp([d.cp.Ev(idint(1:end-1)) d.cp.Ev(idint(1:end-1)+1)]*e/kb/d.T)),[],2);
    
            %% Tunnelling
    
            phin=zeros(d.Nle-1,2);
            phip=zeros(d.Nle-1,2);
            % Jtunnel represent the fraction of the total current related
            % to tunnelling, obtained as phi(z)/(eta+phi)*J.
            % Here, phi(z) corresponds to the integral from the minimum energy 
            % to E(z) in the case of the conduction band, and the opposite in 
            % the valence band. This is the opposite of the notation in the
            % manuscript, where phi(z) is the integral from E(z) to Ecmax,
            % and the fraction of current related to tunnelling is thus
            % (phi-phi(z))/(eta+phi)*J.
            Jntunnel=zeros(idint(end)-1,d.Nle-1);
            Jptunnel=zeros(idint(end)-1,d.Nle-1);

            % For each interface i, tunnelling of electrons and holes is computed
            % The side on which tunnelling happens is checked first
            % Treatment is different depending on the barrier shape
            if options.TunnelDepthMax>0
                % loop on the interface index
                for i=1:d.Nle-1
                    % Electron tunnelling
                    % tunnelling on the right of interface i
                    if (d.chi(i+1)<d.chi(i) && d.cp.Ec(idint(i)+2)<d.cp.Ec(idint(i)+1))
                        Ecmin=max(d.cp.Ec(idint(i)),min(d.cp.Ec(idint(i)+1:idint(i+1))));
                        if Ecmin==d.cp.Ec(idint(i)) % minimum reached at the interface
                            zl=interp1(d.cp.Ec(idint(i)+1:idint(i)+find(d.cp.Ec(idint(i)+1:end)<Ecmin)),d.zea(idint(i)+1:idint(i)+find(d.cp.Ec(idint(i)+1:end)<Ecmin)),Ecmin);
                        else % minimum reached in the bulk
                            zl=d.zea(idint(i)+find(d.cp.Ec(idint(i)+1:end)==Ecmin,1));
                        end
                        zlim=min(zl,d.zea(find(d.zea>=min(d.zea(idint(i)+1)+options.TunnelDepthMax,d.zea(end)),1)));
                        if zlim==zl && Ecmin==d.cp.Ec(idint(i))
                            Ecint=e*[d.cp.Ec(idint(i)+1:find(d.zea<zlim,1,'last'));Ecmin];
                            zint=[d.zea(idint(i)+1:find(d.zea<zlim,1,'last'));zlim];
                        else
                            Ecint=e*d.cp.Ec(idint(i)+1:find(d.zea==zlim,1));
                            zint=d.zea(idint(i)+1:find(d.zea==zlim,1));
                        end
                        logTz=sqrt(2*d.mx(i+1,1)*max((Ecint-Ecint'),0));
                        logT=-2/hb*trapz(zint,logTz,1).';
                        intTunnel=exp((Ecint(1)-Ecint)/kb/d.T+logT);
                        phin(i,:)=d.etat(i,1)*[0 -1/kb/d.T*trapz(Ecint,intTunnel)];
                        Jntunnel(idint(i)+1:idint(i+1)-1,i)=interp1(zint,d.etat(i,1)*-1/kb/d.T*flip(cumtrapz(-flip(Ecint),flip(intTunnel)))/(d.eta(i,1)+phin(i,2)),d.zea(idint(i)+1:idint(i+1)-1)/2+d.zea(idint(i)+2:idint(i+1))/2);
                        lastDer=2*abs((Jntunnel(idint(i+1)-1,i)-Jntunnel(idint(i+1)-2,i))/(d.zea(idint(i+1))-d.zea(idint(i+1)-2))/Jntunnel(idint(i)+1,i));
                        if Ecmin~=d.cp.Ec(idint(i)) && lastDer>=1e5 % relative variation of 1e-4 of max value per nm
                            warning("Tunnelling not well-calculated due to too thin layer (n, layer "+(i+1)+", interface "+i+").");
                        end
                    % tunnelling on the left of interface i
                    elseif (d.chi(i)<d.chi(i+1) && d.cp.Ec(idint(i)-1)<d.cp.Ec(idint(i)))
                        Ecmin=max(d.cp.Ec(idint(i)+1),min(d.cp.Ec((i>1)*idint(max(1,i-1))+1:idint(i))));
                        if Ecmin==d.cp.Ec(idint(i)+1)
                            zl=interp1(d.cp.Ec((i>1)*idint(max(1,i-1))+find(d.cp.Ec((i>1)*idint(max(1,i-1))+1:idint(i))<Ecmin,1,'last'):idint(i)),d.zea((i>1)*idint(max(1,i-1))+find(d.cp.Ec((i>1)*idint(max(1,i-1))+1:idint(i))<Ecmin,1,'last'):idint(i)),Ecmin);
                        else
                            zl=d.zea(find(d.cp.Ec==Ecmin,1));
                        end
                        zlim=max(zl,d.zea(find(d.zea>max(d.zea(idint(i))-options.TunnelDepthMax,0),1)-1));
                        if zlim==zl && Ecmin==d.cp.Ec(idint(i)+1)
                            Ecint=e*[Ecmin;d.cp.Ec(find(d.zea>zlim,1):idint(i))];
                            zint=[zlim;d.zea(find(d.zea>zlim,1):idint(i))];
                        else
                            Ecint=e*d.cp.Ec(find(d.zea==zlim,1):idint(i));
                            zint=d.zea(find(d.zea==zlim,1):idint(i));
                        end
                        logTz=sqrt(2*d.mx(i,1)*max((Ecint-Ecint'),0));
                        logT=-2/hb*trapz(zint,logTz).';
                        intTunnel=exp((Ecint(end)-Ecint)/kb/d.T+logT);
                        phin(i,:)=d.etat(i,1)*[1/kb/d.T*trapz(Ecint,intTunnel) 0];
                        Jntunnel((i>1)*idint(max(1,i-1))+1:idint(i)-1,i)=interp1(zint,d.etat(i,1)*1/kb/d.T*cumtrapz(Ecint,intTunnel)/(d.eta(i,1)+phin(i,1)),d.zea((i>1)*idint(max(1,i-1))+1:idint(i)-1)/2+d.zea((i>1)*idint(max(1,i-1))+2:idint(i))/2);
                        lastDer=2*abs((Jntunnel((i>1)*idint(max(1,i-1))+2,i)-Jntunnel((i>1)*idint(max(1,i-1))+1,i))/(d.zea((i>1)*idint(max(1,i-1))+3)-d.zea((i>1)*idint(max(1,i-1))+1))/Jntunnel(idint(i)-1,i));
                        if Ecmin~=d.cp.Ec(idint(i)+1) && lastDer>=1e5 % relative variation of 1e-4 of max value per nm
                            warning("Tunnelling not well-calculated due to too thin layer (n, layer "+i+", interface "+i+").");
                        end
                    % no tunnelling at interface i
                    else
                        phin(i,:)=0;
                    end
                    % Hole tunnelling
                    % tunnelling on the right of interface i
                    if (d.chi(i)+d.Eg(i)<d.chi(i+1)+d.Eg(i+1) && d.cp.Ev(idint(i)+1)<d.cp.Ev(idint(i)+2))
                        Evmax=min(d.cp.Ev(idint(i)),max(d.cp.Ev(idint(i)+1:idint(i+1))));
                        if Evmax==d.cp.Ev(idint(i))
                            zl=interp1(d.cp.Ev(idint(i)+1:idint(i+1)),d.zea(idint(i)+1:idint(i+1)),Evmax);
                        else
                            zl=d.zea(find(d.cp.Ev==Evmax,1));
                        end
                        zlim=min(zl,d.zea(find(d.zea>=min(d.zea(idint(i)+1)+options.TunnelDepthMax,d.zea(end)),1)));
                        if zlim==zl && Evmax==d.cp.Ev(idint(i))
                            Evint=e*[d.cp.Ev(idint(i)+1:find(d.zea<zlim,1,'last'));Evmax];
                            zint=[d.zea(idint(i)+1:find(d.zea<zlim,1,'last'));zlim];
                        else
                            Evint=e*d.cp.Ev(idint(i)+1:find(d.zea==zlim,1));
                            zint=d.zea(idint(i)+1:find(d.zea==zlim,1));
                        end
                        logTz=sqrt(2*d.mx(i+1,2)*max((Evint'-Evint),0));
                        logT=-2/hb*trapz(zint,logTz).';
                        intTunnel=exp((Evint-Evint(1))/kb/d.T+logT);
                        phip(i,:)=d.etat(i,2)*[0 1/kb/d.T*trapz(Evint,intTunnel)];
                        Jptunnel(idint(i)+1:idint(i+1)-1,i)=interp1(zint,d.etat(i,2)*1/kb/d.T*flip(cumtrapz(-flip(Evint),flip(intTunnel)))/(d.eta(i,2)+phip(i,2)),d.zea(idint(i)+1:idint(i+1)-1)/2+d.zea(idint(i)+2:idint(i+1))/2);
                        lastDer=2*abs((Jptunnel(idint(i+1)-1,i)-Jptunnel(idint(i+1)-2,i))/(d.zea(idint(i+1))-d.zea(idint(i+1)-2))/Jptunnel(idint(i)+1,i));
                        if Evmax~=d.cp.Ev(idint(i)) && lastDer>=1e5 % relative variation of 1e-4 of max value per nm
                            warning("Tunnelling not well-calculated due to too thin layer (p, layer "+(i+1)+", interface "+i+").");
                        end
                    % tunnelling on the left of interface i
                    elseif (d.chi(i+1)+d.Eg(i+1)<d.chi(i)+d.Eg(i) && d.cp.Ev(idint(i))<d.cp.Ev(idint(i)-1))
                        Evmax=min(d.cp.Ev(idint(i)+1),max(d.cp.Ev((i>1)*idint(max(1,i-1))+1:idint(i))));
                        if Evmax==d.cp.Ev(idint(i)+1)
                            zl=interp1(d.cp.Ev((i>1)*idint(max(1,i-1))+1:idint(i)),d.zea((i>1)*idint(max(1,i-1))+1:idint(i)),Evmax);
                        else
                            zl=d.zea(find(d.cp.Ev==Evmax,1));
                        end
                        zlim=max(zl,d.zea(find(d.zea>max(d.zea(idint(i))-options.TunnelDepthMax,0),1)-1));
                        if zlim==zl && Evmax==d.cp.Ev(idint(i)+1)
                            Evint=e*[Evmax;d.cp.Ev(find(d.zea>zlim,1):idint(i))];
                            zint=[zlim;d.zea(find(d.zea>zlim,1):idint(i))];
                        else
                            Evint=e*d.cp.Ev(find(d.zea==zlim,1,'last'):idint(i));
                            zint=d.zea(find(d.zea==zlim,1,'last'):idint(i));
                        end
                        logTz=sqrt(2*d.mx(i,2)*max((Evint'-Evint),0));
                        logT=-2/hb*trapz(zint,logTz).';
                        intTunnel=exp((Evint-Evint(end))/kb/d.T+logT);
                        phip(i,:)=d.etat(i,2)*[-1/kb/d.T*trapz(Evint,intTunnel) 0];
                        Jptunnel((i>1)*idint(max(1,i-1))+1:idint(i)-1,i)=interp1(zint,d.etat(i,2)*-1/kb/d.T*cumtrapz(Evint,intTunnel)/(d.eta(i,2)+phip(i,1)),d.zea((i>1)*idint(max(1,i-1))+1:idint(i)-1)/2+d.zea((i>1)*idint(max(1,i-1))+2:idint(i))/2);
                        lastDer=2*abs((Jptunnel((i>1)*idint(max(1,i-1))+2,i)-Jptunnel((i>1)*idint(max(1,i-1))+1,i))/(d.zea((i>1)*idint(max(1,i-1))+3)-d.zea((i>1)*idint(max(1,i-1))+1))/Jptunnel(idint(i)-1,i));
                        if Evmax~=d.cp.Ev(idint(i)+1) && lastDer>=1e5 % relative variation of 1e-4 of max value per nm
                            warning("Tunnelling not well-calculated due to too thin layer (p, layer "+i+", interface "+i+").");
                        end
                    % no tunnelling
                    else
                        phip(i,:)=0;
                    end
                end
                Jntunnel(isnan(Jntunnel(:)))=0;
                Jptunnel(isnan(Jptunnel(:)))=0;
            end
    
            if exist('zint','var')
                if any(zint~=sort(zint))
                    error('Mistake in tunnelling grid.');
                end
            end
    
            %% Current density (drift-diffusion,continuity,thermionic,tunneling)
            
            % initialisation
            d.cp.Jn=zeros(idint(end)-1,1);
            d.cp.Jp=zeros(idint(end)-1,1);
            
            % Resolution of continuity equation
            % Different resolution depending on where carriers are in majority or minority
            if d.DopSign(d.NleRange(1))<0
                % Electrons
                d.cp.Jn(end)=-e*d.cp.Rs(2)-d.dza(end)/2*(.25*RG(end-1)+.75*RG(end));
                for j=d.Nle-1:-1:1
                    d.cp.Jn(idint(j+1)-2:-1:idint(j)+1)=d.cp.Jn(idint(j+1)-1)-cumsum((d.dza(idint(j+1)-2:-1:idint(j)+1)+d.dza(idint(j+1)-1:-1:idint(j)+2))/2.*(.25*RG(idint(j+1)-2:-1:idint(j)+1)+.5*RG(idint(j+1)-1:-1:idint(j)+2)+.25*RG(idint(j+1):-1:idint(j)+3)));
                    d.cp.Jn(idint(j))=d.cp.Jn(idint(j)+1)-d.dza(idint(j)+1)/2*(.75*RG(idint(j)+1)+.25*RG(idint(j)+2));
                    d.cp.Jn(idint(j)-1)=d.cp.Jn(idint(j))-d.dza(idint(j)-1)/2*(.25*RG(idint(j)-1)+.75*RG(idint(j)));
                end
                d.cp.Jn(idint(1)-2:-1:1)=d.cp.Jn(idint(1)-1)-cumsum((d.dza(idint(1)-2:-1:1)+d.dza(idint(1)-1:-1:2))/2.*(.25*RG(idint(1)-2:-1:1)+.5*RG(idint(1)-1:-1:2)+.25*RG(idint(1):-1:3)));
    
                % Holes
                d.cp.Jp(1)=-e*d.cp.Rs(1)-d.dza(1)/2*(.75*RG(1)+.25*RG(2));
                d.cp.Jp(2:idint(1)-1)=d.cp.Jp(1)-cumsum((d.dza(1:idint(1)-2)+d.dza(2:idint(1)-1))/2.*(.25*RG(1:idint(1)-2)+.5*RG(2:idint(1)-1)+.25*RG(3:idint(1))));
                for j=1:d.Nle-1
                    d.cp.Jp(idint(j))=d.cp.Jp(idint(j)-1)-d.dza(idint(j)-1)/2*(.25*RG(idint(j)-1)+.75*RG(idint(j)));
                    d.cp.Jp(idint(j)+1)=d.cp.Jp(idint(j))-d.dza(idint(j)+1)/2*(.75*RG(idint(j)+1)+.25*RG(idint(j)+2));
                    d.cp.Jp(idint(j)+2:idint(j+1)-1)=d.cp.Jp(idint(j)+1)-cumsum((d.dza(idint(j)+1:idint(j+1)-2)+d.dza(idint(j)+2:idint(j+1)-1))/2.*(.25*RG(idint(j)+1:idint(j+1)-2)+.5*RG(idint(j)+2:idint(j+1)-1)+.25*RG(idint(j)+3:idint(j+1))));
                end
            else 
                % Electrons
                d.cp.Jn(1)=e*d.cp.Rs(1)+d.dza(1)/2*(.75*RG(1)+.25*RG(2));
                d.cp.Jn(2:idint(1)-1)=d.cp.Jn(1)+cumsum((d.dza(1:idint(1)-2)+d.dza(2:idint(1)-1))/2.*(.25*RG(1:idint(1)-2)+.5*RG(2:idint(1)-1)+.25*RG(3:idint(1))));
                for j=1:d.Nle-1
                    d.cp.Jn(idint(j))=d.cp.Jn(idint(j)-1)+d.dza(idint(j)-1)/2*(.25*RG(idint(j)-1)+.75*RG(idint(j)));
                    d.cp.Jn(idint(j)+1)=d.cp.Jn(idint(j))+d.dza(idint(j)+1)/2*(.75*RG(idint(j)+1)+.25*RG(idint(j)+2));
                    d.cp.Jn(idint(j)+2:idint(j+1)-1)=d.cp.Jn(idint(j)+1)+cumsum((d.dza(idint(j)+1:idint(j+1)-2)+d.dza(idint(j)+2:idint(j+1)-1))/2.*(.25*RG(idint(j)+1:idint(j+1)-2)+.5*RG(idint(j)+2:idint(j+1)-1)+.25*RG(idint(j)+3:idint(j+1))));
                end
    
                % Holes
                d.cp.Jp(end)=e*d.cp.Rs(2)+d.dza(end)/2*(.25*RG(end-1)+.75*RG(end));
                for j=d.Nle-1:-1:1
                    d.cp.Jp(idint(j+1)-2:-1:idint(j)+1)=d.cp.Jp(idint(j+1)-1)+cumsum((d.dza(idint(j+1)-2:-1:idint(j)+1)+d.dza(idint(j+1)-1:-1:idint(j)+2))/2.*(.25*RG(idint(j+1)-2:-1:idint(j)+1)+.5*RG(idint(j+1)-1:-1:idint(j)+2)+.25*RG(idint(j+1):-1:idint(j)+3)));
                    d.cp.Jp(idint(j))=d.cp.Jp(idint(j)+1)+d.dza(idint(j)+1)/2*(.75*RG(idint(j)+1)+.25*RG(idint(j)+2));
                    d.cp.Jp(idint(j)-1)=d.cp.Jp(idint(j))+d.dza(idint(j)-1)/2*(.25*RG(idint(j)-1)+.75*RG(idint(j)));
                end
                d.cp.Jp(idint(1)-2:-1:1)=d.cp.Jp(idint(1)-1)+cumsum((d.dza(idint(1)-2:-1:1)+d.dza(idint(1)-1:-1:2))/2.*(.25*RG(idint(1)-2:-1:1)+.5*RG(idint(1)-1:-1:2)+.25*RG(idint(1):-1:3)));
            end
            % Correction due to tunneling
            d.cp.Jn=d.cp.Jn-sum(Jntunnel.*repmat(d.cp.Jn(idint(1:end-1))',length(d.cp.Jn),1),2);
            d.cp.Jp=d.cp.Jp-sum(Jptunnel.*repmat(d.cp.Jp(idint(1:end-1))',length(d.cp.Jp),1),2);
    
            %% Quasi-Fermi levels (from current density and susceptibilities)
            
            % General expression of quasi-Fermi levels variations inside layers
            DEfn=d.cp.Jn./Gn;
            DEfp=d.cp.Jp./Gp;
            
            % Boundary conditions at interfaces
            for j=1:d.Nle-1
                je=d.NleRange(j);
                if (strcmp(d.mat(je),d.mat(je+1)) && d.x(je)==d.x(je+1) && d.ni(j)==d.ni(j+1)) || ~options.TIatInterfaces
                    DEfn(idint(j))=0;
                    DEfp(idint(j))=0;
                else
                    DEfn(idint(j))=d.cp.Jn(idint(j))./Gnth(j)./(d.eta(j,1)+sum(phin(j,:)));
                    DEfp(idint(j))=d.cp.Jp(idint(j))./Gpth(j)./(d.eta(j,2)+sum(phip(j,:)));
                end
            end
            % Complete quasi-Fermi levels
            if d.DopSign(d.NleRange(1))<0
                d.cp.Efnb=exp(Vfb)+[0;cumsum(DEfn)];
                d.cp.Efpb=1-[flip(cumsum(flip(DEfp)));0];
            else
                d.cp.Efnb=exp(Vfb)-[flip(cumsum(flip(DEfn)));0];
                d.cp.Efpb=1+[0;cumsum(DEfp)];
            end

            % for nonphysical results at the first iteration, the iterative process is set
            % convergence is largely improved by alterning between standard iterations with 
            % small weights (to prevent divergence) and punctual leaps with high weights given at 
            % the new result (to reduce comp. time)
            if (any(d.cp.Efnb<0) || any(d.cp.Efpb<0)) || ~d.physEf(id)
                d.physEf(id)=false;
                % Leap
                if (iterEf>=2 && (abs(errm2n(iterEf)-errm2n(iterEf-1))/errm2n(iterEf)<2e-3 || errm2n(iterEf)<tol2*d.wefn/2) && (abs(errm2p(iterEf)-errm2p(iterEf-1))/errm2p(iterEf)<2e-3 || errm2p(iterEf)<tol2*d.wefp/2) && errm2(iterEf)<1e-2) %|| iterEf>iterEfd+100 %&& id~=10%iterE>=2 && abs(errM(iterE)-errM(iterE-1))/errM(iterE)<1e-2 && errM(iterE)<1e-2
                    iterEfd=iterEf;
                    d.cp.Efnb=EfnbOld2+d.wefleap*(d.cp.Efnb-EfnbOld2);
                    d.cp.Efpb=EfpbOld2+d.wefleap*(d.cp.Efpb-EfpbOld2);               
                % Standard iteration
                else
                    d.cp.Efnb=d.wefn*d.cp.Efnb+(1-d.wefn)*EfnbOld2;
                    d.cp.Efpb=d.wefp*d.cp.Efpb+(1-d.wefp)*EfpbOld2;
                    % if last leap was performed long ago, convergence is not good. 
                    % we change the weights to improve it (starting again the whole computation).
                    if iterEf>iterEfd+100
                        if abs(errm2n(iterEf)-errm2n(iterEf-1))/errm2n(iterEf)>=2e-3 && abs(sum(sign(diff(errm2n(end-10:end)))))<2
                            disp('Reducing weight for convergence of electron quasi-Fermi level (oscillations).');
                            weightFermiChange="n";
                        elseif abs(errm2p(iterEf)-errm2p(iterEf-1))/errm2p(iterEf)>=2e-3 && abs(sum(sign(diff(errm2p(end-10:end)))))<2
                            disp('Reducing weight for convergence of hole quasi-Fermi level (oscillations).');
                            weightFermiChange="p";
                        else
                            disp('Reducing weight for convergence of quasi-Fermi levels (leaps).');
                            weightFermiChange="leap";
                        end
                        calc_state="on";                        
                        break;
                    end
                end
            end
            
            % If some NaN value is obtained, we stop the process.
            if any(isnan(d.cp.Efnb)) || any(isnan(d.cp.Efpb))
        	    warning('Complex or NaN quasi-Fermi level. Decreasing factor w.');
                calc_state="on";
                break;
            end
            
            % New chemical potential
            d=d.setnp;
            for k=1:d.Nle
                mumean(sum(d.Nze(1:k-1))+1:sum(d.Nze(1:k)))=kb*d.T/e*log(mean(d.cp.n(sum(d.Nze(1:k-1))+1:sum(d.Nze(1:k))).*d.cp.p(sum(d.Nze(1:k-1))+1:sum(d.Nze(1:k)))/d.ni(k)^2));
            end

            % Update iteration number and error.
            iterEf=iterEf+1;
            errm2n(iterEf)=max(abs(d.cp.Efnb-EfnbOld2)./abs(EfnbOld2)); %#ok<AGROW> 
            errm2p(iterEf)=max(abs(d.cp.Efpb-EfpbOld2)./abs(EfpbOld2)); %#ok<AGROW> 
            errm2(iterEf)=max([errm2n(iterEf) errm2p(iterEf)]); %#ok<AGROW> 

            % if after few iterations (at start, or after a leap), error is large, it means that it 
            % diverges. Weights are reduced accordingly.
            if iterEf>=5 && errm2n(iterEf)>1 && (iterEf<iterEfd || iterEf>iterEfd+20)
                disp('Reducing weight for convergence of electron quasi-Fermi level (no convergence).');
                calc_state="on";
                weightFermiChange="n";
                break;
            end
            if iterEf>=5 && errm2p(iterEf)>1 && (iterEf<iterEfd || iterEf>iterEfd+20)
                disp('Reducing weight for convergence of hole quasi-Fermi level (no convergence).');
                calc_state="on";
                weightFermiChange="p";
                break;
            end

        end

        if calc_state=="on"
            break;
        end

        if iterEf==iterEfmax
            warning('Exceeded max number of iterations (quasi-Fermi levels).');
            calc_state="on";
            weightFermiChange="leap";
            break
        end

        %% Electrostatic potential (Poisson equation)

        % Poisson - new boundary conditions
        if d.DopSign(d.NleRange(1))==-1 && d.DopSign(d.NleRange(end))==1
            d.cp.Vb(1)=-d.xinv(1)-log(d.cp.Efnb(1))+log(sqrt((d.Ndop(d.NleRange(1))/(2*d.ni(1)))^2+exp(d.Dgv(1))*d.cp.Efpb(1)*d.cp.Efnb(1))+(d.Ndop(d.NleRange(1))/(2*d.ni(1))));
            d.cp.Vb(d.Nzet)=d.xipv(end)+log(d.cp.Efpb(d.Nzet))-log(sqrt((d.Ndop(d.NleRange(end))/(2*d.ni(1)))^2+exp(d.Dgv(end))*d.cp.Efpb(d.Nzet)*d.cp.Efnb(d.Nzet))+(d.Ndop(d.NleRange(end))/(2*d.ni(1))));
        elseif d.DopSign(d.NleRange(1))==1 && d.DopSign(d.NleRange(end))==-1
            d.cp.Vb(1)=d.xipv(1)+log(d.cp.Efpb(1))-log(sqrt((d.Ndop(d.NleRange(1))/(2*d.ni(1)))^2+exp(d.Dgv(1))*d.cp.Efpb(1)*d.cp.Efnb(1))+(d.Ndop(d.NleRange(1))/(2*d.ni(1))));
            d.cp.Vb(d.Nzet)=-d.xinv(end)-log(d.cp.Efnb(d.Nzet))+log(sqrt((d.Ndop(d.NleRange(end))/(2*d.ni(1)))^2+exp(d.Dgv(end))*d.cp.Efpb(d.Nzet)*d.cp.Efnb(d.Nzet))+(d.Ndop(d.NleRange(end))/(2*d.ni(1))));
        end
        
        % Warning if V reaches complex values
        if any(imag(d.cp.Vb)) || any(isnan(d.cp.Vb))
            if calc_state~="on"
        	    warning('Complex or NaN potential. Decreasing factor w.');
                calc_state="on";
            end
            break;
        end

        % Warning if exponential quantities get negative
        if any(d.cp.Efnb<=0) || any(d.cp.Efpb<=0)
            warning('Negative exponential quantity. Decreasing factor w.');
            if all(d.cp.Efnb>0)
                weightFermiChange="p";
            elseif all(d.cp.Efpb>0)
                weightFermiChange="n";
            end
            calc_state="on";
            break;
        end
        
        % Poisson - iterative process
        iterV=0;
        err=1;
        A=[d.K(1:d.Nzet-2,1);0];
        C=[0;d.K(1:d.Nzet-2,3)];
        zint=cumsum(d.Nze(1:end-1));
        Efnb=d.za2z(d.cp.Efnb);
        Efpb=d.za2z(d.cp.Efpb);
        Kn=exp(d.xinv).*Efnb;
        Kp=exp(d.xipv).*Efpb;
        while err>options.tol && iterV<options.iterVmax
            nnorm=Kn.*exp(d.cp.Vb);
            pnorm=Kp.*exp(-d.cp.Vb);
    
            % Direct terms
            B=[1;d.K(1:d.Nzet-2,2)-d.dzb(1:d.Nzet-2).^2.*(nnorm(2:d.Nzet-1)+pnorm(2:d.Nzet-1));1];
            Y=[0;-(d.K(1:d.Nzet-2,1).*d.cp.Vb(1:d.Nzet-2)+d.K(1:d.Nzet-2,2).*d.cp.Vb(2:d.Nzet-1)+d.K(1:d.Nzet-2,3).*d.cp.Vb(3:d.Nzet))+...
                  d.dzb(1:d.Nzet-2).^2.*(nnorm(2:d.Nzet-1)-pnorm(2:d.Nzet-1)+d.Nm(2:d.Nzet-1));0];
    
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

        % Warning if max iteration reached
        if iterV==options.iterVmax
           warning(join(['Max number of iterations reached for V (err=',num2str(err),'). Decreasing factor w.']));
           calc_state="on";
           break;
        end
        
        % Warning if V reaches complex values
        if any(imag(d.cp.Vb)) || any(isnan(d.cp.Vb))
        	warning('Complex or NaN potential. Decreasing factor w.');
            calc_state="on";
            break;
        end

        %% Finishing iteration

        % Ponderation to limit oscillation
        if iterE>=1
            d.cp.Efnb=d.wnext*d.cp.Efnb+(1-d.wnext)*EfnbOld;
            d.cp.Efpb=d.wnext*d.cp.Efpb+(1-d.wnext)*EfpbOld;
            d.cp.Vb=d.wnext*d.cp.Vb+(1-d.wnext)*Vc;
        end

        % New solutions
        Vc=d.cp.Vb;
        d=d.setnp;
        iterE=iterE+1;
        
        errn=max(abs(d.cp.n-nOld)./abs(nOld));
        errp=max(abs(d.cp.p-pOld)./abs(pOld));
        err=max(errn,errp);
        errM(iterE)=err; %#ok<AGROW,NASGU> 

        if iterE>=10 && err>=10
            warning('Not converging. Decreasing factor w.');
            calc_state="on";
            weightFermiChange="next";
            break;
        end

    end

    %% Wrapping up

    d.J(id)=d.cp.Jn(1)+d.cp.Jp(1);
    % if reversed polarization, need to take the opposite value
    if d.DopSign(d.NleRange(1))>0
        d.J(id)=-d.J(id);
    end

    % Warning if max iteration reached
    if iterE==options.iterEmax && err>options.tol
        warning(join(['Max number of iterations reached for E (err=',num2str(err),'). Decreasing factor w.']));
        calc_state="on";
        weightFermiChange="next";
    end
    
    % Stop process if changing weights is forbidden by option
    if calc_state=="on" && ~options.VariableWeights
        calc_state="over - divergence";
    end

    % Change the weights in case of divergence
    if calc_state=="on"
        physEf=d.physEf(id);
        d=dOld;
        if physEf % physical quasi Fermi level -> change global weight
            switch d.wnext
                case .8, d.wnext=.1;
                case .1, calc_state="over - divergence";
                otherwise, d.wnext=.8;
            end
        else % non-physical quasi Fermi level -> change related weights
            switch weightFermiChange
                case "n"
                    d.wefn=d.wefn/2;
                case "p"
                    d.wefp=d.wefp/2;
                case "leap"
                    d.wefleap=d.wefleap/2;
                case "next"
                    switch d.wnext
                        case .8, d.wnext=.1;
                        case .1, calc_state="over - divergence";
                        otherwise, d.wnext=.8;
                    end
                otherwise
                    d.wefn=d.wefn/2;
                    d.wefp=d.wefp/2;
                    d.wefleap=d.wefleap/2;
            end
            % if weights are too small, stop the calculation
            if d.wefn<1e-10 || d.wefp<1e-10 || d.wefleap<.1e-1
                calc_state="over - divergence";
            end
        end
        d.cp.R=zeros(length(d.ze),3);
        if ~isempty(Gem)
            varargout{1}=G;
            varargout{2}=Gem;
        end
        return
    elseif calc_state=="over - divergence"
        if ~isempty(Gem)
            varargout{1}=G;
            varargout{2}=Gem;
        end
        d=dOld;
    else
        if ~isempty(Gem)
            varargout{1}=Geq;
            varargout{2}=Gem+G-Geq;
        end
    end

end

function y=Ber(x)
    %BER Compute the Bernoulli function, using a Taylor expansion if x is
    %too small.
    y=zeros(size(x));
    x1=x(abs(x)>1e-2);
    y(abs(x)>1e-2)=x1./(exp(x1)-1+(abs(x1)<=1e-5));
    x2=x(abs(x)<=1e-2);
    y(abs(x)<=1e-2)=1-x2/2+x2.^2/12;%-x2.^4/720;
end

