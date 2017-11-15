function U = wPCA(data)

X = data';
    means = mean(X,2);
    
    %Centering data
    
    meanX = repmat(means,1,size(data,1));
    X_cent = X-meanX;
    [V,E] = eig(X_cent'*X_cent);
    disp(diag(E))
    
    % Account for reduction in rank of one due to removal
    % of mean previously
    
    E = E(2:size(E),2:size(E));
    V = V(:,2:size(V,2));
        
    U = X*V*E^(-0.5)*E^(-0.5);
end