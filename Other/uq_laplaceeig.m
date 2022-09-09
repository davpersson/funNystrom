function H0 = uq_laplaceeig(kappa,sigma,lambda,T_list)

%Function to create A_pde

N = 40;

h = 1/N;

T = diag(-2*ones(N-1,1)) + diag(ones(N-2,1),-1) + diag(ones(N-2,1),1);
J = zeros(N,N); J(end,end) = 1/2;
That = T - 2*eye(size(T));
Ttilde = diag(-2*ones(N,1)) + diag(ones(N-1,1),-1) + diag(ones(N-1,1),1);
L = kappa*(kron(eye(N),T) + kron(Ttilde,eye(N-1)) - kron(J,That))/h^2 + ...
    lambda*(eye(N*(N-1)) - kron(J,eye(N-1)));

sqrtC0 = -inv((kron(eye(N),T) + kron(Ttilde,eye(N-1)) - kron(J,That)));
%theta = sqrtC0*randn(N*(N-1),1);
I = zeros(1,49);
for k = 1:7
    
    I((7*(k-1)+1):(7*k)) = (5:5:35)+40*(k-1);
    
end

[U,Lambda] = eig(L);

F = [];

for t = T_list
    
    F = [F; U(I,:)*expm(kappa*t*Lambda)*U'];
    
end

% H = F'*F/sigma^2; 
F0 = F*sqrtC0;
H0 = F0'*F0/sigma^2;

% N_obs = length(I)*length(T_list);
% s = svd(H);
% s0 = svd(H0);
% semilogy(s(1:(N_obs)),'b');
% hold on
% semilogy(s0(1:(N_obs)),'r-')
% hold off

end

