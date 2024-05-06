function mu = mobilityCTstruct(N,T,CTm)
%MOBILITYCAUGHEYTHOMAS Computes the particle mobility, using Caughey-Thomas model
   
    mu=CTm.mu_min+(CTm.mu_max.*(300./T).^CTm.theta1-CTm.mu_min)./(1+(N./(CTm.Nref.*(T/300).^CTm.theta2)).^CTm.alpha);
    
end

