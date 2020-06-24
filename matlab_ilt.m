function ilt = matlab_ilt(fun, T, maxFnEvals, method)

global cmeParams;
if isempty(cmeParams)
    cmeParams = jsondecode(fileread('iltcme.json'));
end
if ~exist('method','var')
    method = 'cme';
end

if strcmp(method, 'cme')

    % find the most steep CME satisfying maxFnEvals
    params = cmeParams(1);
    for i=2:length(cmeParams)
        if cmeParams(i).cv2<params.cv2 && cmeParams(i).n+1<=maxFnEvals
            params = cmeParams(i);
        end
    end

    % compute eta and beta parameters
    eta = [params.c*params.mu1, params.a'*params.mu1 + 1i*params.b'*params.mu1];
    beta = [1, 1 + 1i*(1:params.n)*params.omega] * params.mu1;

elseif strcmp(method,'euler')
    
    n_euler = floor((maxFnEvals-1)/2);
    eta = [0.5, ones(1, n_euler), zeros(1, n_euler-1), 2^-n_euler];
    for k = 1:n_euler-1
%        eta(2*n_euler-k + 1) = eta(2*n_euler-k + 2) + 2^-n_euler * nchoosek(n_euler, k);
       eta(2*n_euler-k + 1) = eta(2*n_euler-k + 2) + exp(sum(log(1:n_euler)) - n_euler*log(2) - sum(log(1:k)) - sum(log(1:(n_euler-k))));
    end
    k = 0:2*n_euler;
    beta = n_euler*log(10)/3 + 1i*pi*k;
    eta  = (10^((n_euler)/3))*(1-mod(k, 2)*2) .* eta;  
    
elseif strcmp(method,'gaver')

    if mod(maxFnEvals,2)==1
        maxFnEvals = maxFnEvals - 1;
    end
    ndiv2 = maxFnEvals/2;
    eta = zeros(1,maxFnEvals);
    beta = zeros(1,maxFnEvals);
    for k = 1:maxFnEvals % itration index
        inside_sum = 0.0;
        for j = floor((k+1)/2):min(k,ndiv2) %  eta summation index
%            inside_sum=inside_sum+((j^((ndiv2+1))/factorial(ndiv2))*(nchoosek(ndiv2, j)*nchoosek(2*j, j)*nchoosek(j, k-j)));           
           inside_sum=inside_sum+exp((ndiv2+1)*log(j) - sum(log(1:(ndiv2-j))) + sum(log(1:2*j)) - 2*sum(log(1:j)) - sum(log(1:(k-j))) - sum(log(1:(2*j-k))));
        end
        eta(k)=log(2.0)*(-1)^(k+ndiv2)*inside_sum;
        beta(k) = k * log(2.0);
    end    
    
else
    error('This inverse laplace transform method is unknown. Supported ones: cme, euler, gaver');
end

% common part for all abate-whitt variants
[eta_mesh,T_mesh]= meshgrid(eta, T);
beta_mesh = meshgrid(beta, T);
ilt = 1./reshape(T,[],1) .* sum(real(eta_mesh .* arrayfun(fun, beta_mesh./T_mesh)), 2);

end
