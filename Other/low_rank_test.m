function low_rank_test(U,S,fscalar,repeats,rank_list,filename)

disp(filename)

%Allocate tensors for saving results
trace_error = zeros(length(rank_list),repeats,2);
frobenius_error = zeros(length(rank_list),repeats,2);
trace_error_exact = zeros(length(rank_list),repeats,2);
frobenius_error_exact = zeros(length(rank_list),repeats,2);
trace_optimal = zeros(length(rank_list),1);
frobenius_optimal = zeros(length(rank_list),1);

%Compute important quantities
feigvals = fscalar(diag(S));
fA = U*diag(feigvals)*U';
Afun = @(X) (U*(S*(U'*X)));
fAfun = @(X) fA*X;
matrix_size = size(U,1);
trfA = sum(feigvals);
frobfA = norm(fA,'fro');

%Run tests
outer_iteration = 0;
for r = rank_list
    
    outer_iteration = outer_iteration + 1;
    r
    
    trace_optimal(outer_iteration) = sum(feigvals((r+1):end));
    frobenius_optimal(outer_iteration) = norm(feigvals((r+1):end));
    
    for repetition = 1:repeats
        
        %funNystrom q = 1
        [U,S] = nystrom(matrix_size,Afun,r,1);
        fS = diag(fscalar(diag(S)));
        trace_error(outer_iteration,repetition,1) = trfA - trace(fS);
        frobenius_error(outer_iteration,repetition,1) = norm(fA-U*fS*U','fro');
        
        %funNystrom q = 2
        [U,S] = nystrom(matrix_size,Afun,r,2);
        fS = diag(fscalar(diag(S)));
        trace_error(outer_iteration,repetition,2) = trfA - trace(fS);
        frobenius_error(outer_iteration,repetition,2) = norm(fA-U*fS*U','fro');
        
        %Nystrom q = 1
        [U,S] = nystrom(matrix_size,fAfun,r,1);
        trace_error_exact(outer_iteration,repetition,1) = trfA - trace(S);
        frobenius_error_exact(outer_iteration,repetition,1) = norm(fA-U*S*U','fro');
        
        %Nystrom q = 1
        [U,S] = nystrom(matrix_size,fAfun,r,2);
        trace_error_exact(outer_iteration,repetition,2) = trfA - trace(S);
        frobenius_error_exact(outer_iteration,repetition,2) = norm(fA-U*S*U','fro');
        
    end
    
end

%Save results
save(filename,'rank_list','trace_error','frobenius_error','trace_error_exact',...
    'frobenius_error_exact','trace_optimal','frobenius_optimal','trfA','frobfA');

end