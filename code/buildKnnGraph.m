function A_unsup = buildKnnGraph(X, outputAdjMatURL, top)
    batch = 1000;
    [~,n]=size(X);
    
    if exist(outputAdjMatURL, 'file') ~= 2,
    
    stdArr1 = []; 
    for i=1: ceil(n/batch)
        idx_sub = (i-1)*batch+1:1:min(i*batch,n);
        subX = X(:,idx_sub);
        tmpDist = bsxfun(@plus,dot(subX,subX,1)',dot(X,X,1))-2*(subX'*X);         
        stdArr1(end+1) = quantile(tmpDist(:),0.05);
    end
    bandwidth = median(stdArr1); 
    
    if bandwidth == 0,
        stdArr1 = zeros(1,n); 
        for i=1: ceil(n/batch)
            idx_sub = (i-1)*batch+1:1:min(i*batch,n);
            subX = X(:,idx_sub);
            tmpDist = sqrt(bsxfun(@plus,dot(subX,subX,1)',dot(X,X,1))-2*(subX'*X)); 
            for rowIdx = 1: size(tmpDist,1)
                realRow = (i-1)*batch + rowIdx;
                tmpRow = tmpDist(rowIdx, :); stdArr1(realRow) = 1.06*std(tmpRow)/(length(tmpRow) ^ 0.2);
            end
        end
        bandwidth = mean(stdArr1);
    end
    
    fid = fopen(outputAdjMatURL,'w');
    for i=1: ceil(n/batch)
        idx_sub = (i-1)*batch+1:1:min(i*batch,n);
        subX = X(:,idx_sub);
        tmpDist = sqrt(bsxfun(@plus,dot(subX,subX,1)',dot(X,X,1))-2*(subX'*X)); 

        for rowIdx = 1: size(tmpDist,1)
            realRow = (i-1)*batch + rowIdx;
            tmpRow = tmpDist(rowIdx, :); tmpRow(tmpRow==0) = Inf; 
            [B,I] = sort(tmpRow);
            for j=1:top, realCol = I(j); fprintf(fid,'%d\t%d\t%f\n', realRow, realCol, exp(-B(j)/bandwidth)); end
        end
    end
    fclose(fid);
   
    end
    
    weightMat = spconvert(load(outputAdjMatURL));
    if size(weightMat,1) ~= n || size(weightMat,2) ~= n, weightMat(n,n) = 0; end
    
    A_unsup_bin = (weightMat > 0); A_unsup_bin = A_unsup_bin + A_unsup_bin';
    A_unsup = weightMat + weightMat';
    A_unsup(A_unsup>0) = A_unsup(A_unsup>0) ./ A_unsup_bin(A_unsup_bin>0); clear weightMat A_unsup_bin stdArr1;
end