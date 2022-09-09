function err = error_measure_fun(A,B,error_measure)

if strcmp(error_measure,'nuclear_norm')
    
    %Requires A >= B
    err = trace(A)-trace(B);
    
    if err < 0
        
        err = inf;
        
    end
    
elseif strcmp(error_measure,'frobenius_norm')
    
    err = norm(A-B,'fro');
    
else
    
    error('error_measure must be nuclear_norm or frobenius_norm')

end

end