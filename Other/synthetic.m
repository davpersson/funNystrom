function A = synthetic(matrix_size,parameter,multiple,decay)

if strcmp(decay,'algebraic')
    
    Q = orth(randn(matrix_size,matrix_size));
    D = multiple*diag((1:matrix_size).^(-parameter));
    A = Q*D*Q';
    
elseif strcmp(decay,'exponential')
    
    Q = orth(randn(matrix_size,matrix_size));
    D = multiple*diag(exp(-(1:matrix_size)/parameter));
    A = Q*D*Q';
    
else
    
    error('Input decay not valid')
    
end

end