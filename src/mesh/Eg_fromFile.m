function EgMat = Eg_fromFile(d)
%EG_FROMFILE Gives the bandgap of each material used in the device given as
%input (structure d), using the field fmat
% It uses the structure d, rather than the activeComponent instance, as
% the current function is used to mesh the frequency which is later needed
% to generate the activeComponent.

    arguments
        d (1,:) struct % basic structure defining the component
    end

    EgMat=zeros(size(d));

    for k=1:length(d)
        name=func2str(d(k).fmat);
        switch name
            case 'Properties_AlGaAs'
                [~,EgMat(k)]=epsAlGaAsSb_GonzalezCuevas(d(k).x,1,1e15,d(k).T);
            case 'Properties_InGaP'
                [~,EgMat(k)]=epsInGaP_Schubert(d(k).x,1e15,d(k).T);
            case 'Properties_InGaAs'
                [~,EgMat(k)]=epsInGaAsSb_GonzalezCuevas(d(k).x,1,1e15,d(k).T);
            otherwise
                error(join(['Properties function named "',name,'" not found.']));
        end
    end

end