function ilt = inv_cme(fun, T, maxFnEvals)

global cmeParams;
if isempty(cmeParams)
    cmeParams = jsondecode(fileread('iltcme.json'));
end


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

% common part for all abate-whitt variants
[eta_mesh,T_mesh]= meshgrid(eta, T);
beta_mesh = meshgrid(beta, T);
ilt = 1./reshape(T,[],1) .* sum(real(eta_mesh .* arrayfun(fun, beta_mesh./T_mesh)), 2);

end
