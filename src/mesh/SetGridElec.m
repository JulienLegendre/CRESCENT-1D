function d = SetGridElec(d)
%SETGRIDELEC Set the electrical grid
    arguments
        d (1,1) activeComponent % basic structure defining the component
    end

    % Initialisation
    d.Nzet=sum(d.Nze)-(d.Nle-1);
    d.zel=cell(d.Nle,1);
    d.ze=zeros(d.Nzet,1);
    d.zea=zeros(sum(d.Nze),1);

    dz=1e-12; % spatial step on each side of an interface, should be small to improve convergence
    k0=d.NleRange(1); % layer index in the component, of the first layer in the electrical system

    % Set grid for each layer
    for k=1:d.Nle
        ke=d.NleRange(k); % layer index in the component, of the kth layer in the electrical system
        d.zel{k}=[0;linspace(dz,d.tj(ke)-dz,d.Nze(k)-2)';d.tj(ke)];
        d.ze(sum(d.Nze(1:k-1))+3-k:sum(d.Nze(1:k))+1-k)=(k>1)*sum(d.tj(k0:ke-1))+d.zel{k}(2:d.Nze(k));
        d.zea(sum(d.Nze(1:k-1))+1:sum(d.Nze(1:k)))=sum(d.tj(k0:ke-1))+d.zel{k}(1:d.Nze(k));
    end

    % Additional quantities (spatial step)
    d.dz=d.ze(2:d.Nzet)-d.ze(1:d.Nzet-1);
    d.dza=d.zea(2:end)-d.zea(1:end-1);

end
