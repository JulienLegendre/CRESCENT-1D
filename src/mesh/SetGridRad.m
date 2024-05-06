function d = SetGridRad(d)
%SETGRIDRAD Set the radiative grid
    arguments
        d (:,1) component % basic structure defining the component
    end
    
    % Convention used to treat semi-infinite layers
    temp_tj=d.tj;
    for k=[1 d.Nl]
        if isinf(d.tj(k))
            d.Nzr(k)=2;
            d.tj(k)=0;
        end
    end
    
    % Initialisation
    d.Nzrt=sum(d.Nzr)-(d.Nl-1);
    d.zrl=cell(d.Nl,1);
    d.zr=zeros(d.Nzrt,1);
    d.zra=zeros(sum(d.Nzr),1);
    d.zrange=cell(d.Nl,1);

    % For each layer
    for k=1:d.Nl
        d.zrl{k}=linspace(0,d.tj(k),d.Nzr(k))';
        d.zr(sum(d.Nzr(1:k-1))+3-k:sum(d.Nzr(1:k))+1-k)=(k>1)*sum(d.tj(1:k-1))+d.zrl{k}(2:d.Nzr(k));
        d.zra(sum(d.Nzr(1:k-1))+1:sum(d.Nzr(1:k)))=sum(d.tj(1:k-1))+d.zrl{k}(1:d.Nzr(k));
        d.zrange{k}=sum(d.Nzr(1:k-1))+2-k:sum(d.Nzr(1:k))-k;
    end

    % Additional quantities (spatial step)
    d.dzr=d.zr(2:end)-d.zr(1:end-1);

    % Reset layer thicknesses
    d.tj=temp_tj;

end
