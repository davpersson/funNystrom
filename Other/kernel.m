function A = kernel(points,type,parameter_cell)

%Function to create kernel matrices

matrix_size = length(points(1,:));
A = zeros(matrix_size,matrix_size);

if strcmp(type,'square_exponential')
    
    sigma2 = parameter_cell{1};
    K = @(x,y) exp(-norm(x-y)^2/(2*sigma2));
    
elseif strcmp(type,'matern')
    
    nu = parameter_cell{1}; alpha = parameter_cell{2};
    K = @(x,y) sqrt(pi)*((alpha*norm(x-y))^(nu)*besselk(nu,alpha*norm(x-y)))/(2^(nu-1)*alpha^(2*nu)*gamma(nu+0.5));
    
elseif strcmp(type,'polynomial')
    
    g = parameter_cell{1}; c = parameter_cell{2}; p = parameter_cell{3};
    K = @(x,y) (g*x'*y + c)^p;
    
else
    
    error('type not valid')
    
end

if strcmp(type,'matern')

    for row = 1:matrix_size

        for column = 1:matrix_size

            if row == column
                
                A(row,column) = sqrt(pi)*gamma(nu)/(gamma(nu + 1/2)*alpha^(2*nu));

            else

                A(row,column) = K(points(:,row),points(:,column));

            end

        end

    end

else

    for row = 1:matrix_size

        for column = 1:matrix_size

            A(row,column) = K(points(:,row),points(:,column));

        end

    end

end

end