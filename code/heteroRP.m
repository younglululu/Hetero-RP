function [X_reweigh, weightVec] = heteroRP(X, attrMat, repulMat)
    [m,n]=size(X); eta = 1;
    
    Dcol_attr = full(sum(attrMat,2));
    L_attr = spdiags(Dcol_attr,0,n,n)-attrMat; clear Dcol_attr;
    Y_attr = getY(X, L_attr); clear L_attr; %Y_attr = diag(X*L_attr*X');
    Y = Y_attr;
    fprintf('finish loading attrMat \n');
    
    if ~isempty(repulMat)
        Dcol_repul = full(sum(repulMat,2));
        L_repul = spdiags(Dcol_repul,0,n,n)-repulMat; clear Dcol_repul;
        Y_repul = getY(X, L_repul); clear L_repul; %Y_repul = diag(X*L_repul*X'); 
        Y = Y_attr - eta*Y_repul;
        fprintf('finish loading repulMat \n');
    end
    
    selectedIdx = find(Y > 0); selectedLen = length(selectedIdx)
    selectedY = Y(selectedIdx); 
    
    r = min(m,n); currSigma = 10; 
    B = tinv(1-sqrt(m)/(2*r*log(r)),m-1);
    lambda_0 = B / sqrt(m-1+B*B); 

    tol = 1e-6; 
    iter = 1; maxIter = 20; diffVal = inf;
    currCoeff = zeros(selectedLen,1);
    
    while diffVal > tol && iter < maxIter
        iter = iter + 1;
        cvx_begin quiet
           cvx_precision high
           variable x(selectedLen);
           minimize(selectedY'*square(x+1) + 2*lambda_0*currSigma*sum_square(x))
           subject to
                x >= -1;
                sum(x) == 0;
        cvx_end
        newCoeff = x;
        if isnan(newCoeff), continue; end
        
        diffVal = sum(abs((newCoeff-currCoeff))); fprintf('iter= %d \t diff= %f \n', iter-1,diffVal);
        currSigma = sqrt(selectedY'*square(newCoeff+1))*sqrt(m); 
        currCoeff = newCoeff;
    end
    weightVec = ones(m,1); weightVec(selectedIdx) = currCoeff + 1;
    X_reweigh = X .* repmat(weightVec, 1, n);
    
    clear Y_attr Y_repul Y currCoeff newCoeff;
end

function Y = getY(matX, matL)
    [m,n]=size(matX); 
    if m < n,
        Y = diag((matX*matL)*matX');
    else
        tmpMat = matX*matL;
        Y = sum(tmpMat .* matX, 2); clear tmpMat;
    end
end