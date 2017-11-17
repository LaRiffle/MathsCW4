function W = LDA(data, labels)
    
    X = data.';
    X= X - repmat(mean(X, 2), 1, size(X,2));
    C = size(unique(labels), 1); %number of labels
    y = zeros(C);
    pos = 1;
    for i = 1:C
        y(i) = sum(labels==labels(pos));
        pos = pos + y(i);
    end
    M = 1./y(1) * ones(y(1), y(1));
    for i = 2:C
        M2 =  1./y(1) * ones(y(i), y(i));
        M = blkdiag(M, M2);
    end
    Nt = size(M,1);
    [V,D] = eig((eye(Nt, Nt) - M) * (X.' * X) * (eye(Nt, Nt) - M));
    [D,I] = sort(diag(D));
    V = V(:, I);
    D = real(diag(D(C+2:Nt)));
    V = real(V(:,1:Nt-(C+1)));
    U =X*V*D^(-0.5);
    Xtilde = U.' * X * M;
    Q = pca(Xtilde.');
    W = U * Q;
end
