function [z,w,ind,stats,stats_T] = sketchAAA(GZ,Z,deg,opts)
% Input (mandatory)
%   GZ: An ell x M matrix that contains ell sketched surrogates of F(Z) 
%       where F: C -> C^N is the vector-valued function to be approximated
%       on the sampling set Z of M points in C.
%       Construct GZ by XXX.m
%    Z: Discretisation of Omega, used to evaluate approximation errors on
%  deg: max degree of rational approximant to construct
%
% Optional fields in opts:
%  tol:        stopping tolerance (default is 1e-10) 
%  max_norm:   tolerance is max norm, otherwise Frob norm (default is true) 
%  svd_update: update the svd (default is true) 
%
% Output
%  z,w: the nodes and weights of the rational approximant of F
%  ind: the indices in Z that correspond to z
%  stats.err_g: errors of approx on the surrogates g
%
% Part of sketchAAA by Guettel, Kressner, and Vandereycken (2022).
% Code based on SV-AAA from autoCORK [1] by Lietaert et.al.,
% Automatic rational approximation and linearization of nonlinear eigenvalue 
% problems, IMA J. Numer. Anal., (2021). See %%% MODIFICATION for details.
%
% [1] https://people.cs.kuleuven.be/karl.meerbergen/files/aaa/autoCORK.zip


% BV - 11.2022


if ~exist('opts','var')
    opts = [];
end
[max_norm, tol, svd_update] = parse_options(opts);

stats_T.A = 0; stats_T.svd = 0; stats_T.err_F = 0; stats.rest = 0;
T_total = tic();

[nG,M] = size(GZ);
GG = GZ.';
G = GG(:);

R = mean(G); 
mmax = deg; 

errvec = zeros(mmax,1);
z = zeros(mmax,1);
g = zeros(mmax,nG);
ind = zeros(mmax,nG);
C = zeros(M,0);
% Set-valued AAA on surrogates gives support pts 'zval' and weights 'wval'
%%%%%%%%%%%%%%%%%%%%%%%%%%%% AAA SV 

if svd_update
    % Left scaling matrix: not needed if updating svd, use repmat instead.
    SF = spdiags(GG, 0:-M:-M*(nG-1), M*nG, M);
end

% Initialize values    
H = zeros(mmax,mmax-1);
S = zeros(mmax,mmax-1);    
Q = zeros(M*nG,0);        

% Relative tolerance
normG = compute_norm(GG, max_norm);

%%% From AAA code
J = 1:nG*M;
[~, loc] = compute_norm(R-GG, max_norm);

% AAA iteration: stacked version of AAA, one block under each other.
for m = 1:mmax
    loc = mod(loc,M);
    ind(m,:) =  loc + (M*(loc==0):M:(nG-1+(loc==0))*M);  % Get indices of the z_i
    z(m) = Z(ind(m,1));                           % Add interpolation point
    g(m,:) = G(ind(m,:));                         % Add function values
    % no updating
    for k=1:nG
        J(J == ind(m,k)) = []; 
    end

    C(:,end+1) = 1./(Z - z(m));                   % Get column of the Cauchy matrix.
    C(ind(1:m,1),m) = 0;                          % Set the selected elements to 0

    % Compute weights:
    if svd_update
        v = C(:,m)*g(m,:);                            % Compute the next vector of the basis.
        v = G .* repmat(C(:,m),[nG, 1]) - v(:);  % instead of v = SF*C(:,m)-v(:);
                
        % Update H and S to compensate for the removal of the rows
        q = Q(ind(m,:),1:m-1);
        q = q*S(1:m-1,1:m-1);
        ee = eye(m-1,m-1)-q'*q;
        ee(1:size(ee,1)+1:end) = real(diag(ee));
        Si = chol(ee);
        H(1:m-1,1:m-1) = Si*H(1:m-1,1:m-1);
        S(1:m-1,1:m-1) = S(1:m-1,1:m-1)/Si;
        S(m,m) = 1;
        Q(ind(1:m,:),:) = 0;
        
        nv = norm(v);
        H(1:m-1,m) = Q'*v;
        H(1:m-1,m) = S(1:m-1,1:m-1)'*H(1:m-1,m);
        HH = S(1:m-1,1:m-1)*H(1:m-1,m);
        v = v-(Q*HH);
        H(m,m) = norm(v);
        % Reorthoganlization is necessary for higher precision
        it = 0;
        while (it < 3) && (H(m,m) < 1/sqrt(2)*nv)
            h_new = S(1:m-1,1:m-1)'*(Q'*v);
            v = v - Q*(S(1:m-1,1:m-1)*h_new);
            H(1:m-1,m) = H(1:m-1,m) + h_new;
            nv = H(m,m);
            H(m,m) = norm(v);
            it = it+1;            
        end        
        v = v/H(m,m);
        
        % Add v
        Q(:,end+1) = v;
    
        % Solve small least squares problem with H
        [~,~,V] = svd(H(1:m,1:m));
        w = V(:,end);
    else
        kk = 0;    
        for k=1:nG
            Sf = diag(g(1:m,k)); 
            SF = spdiags(GG(:,k), 0, M, M);
            A(1+kk:M+kk,1:m) = SF*C - C*Sf; 
            kk = kk+M;
        end    
        [~, ~, V] = svd(A(J,:), 0);         % Reduced SVD.
        w = V(:,m);
    end
        
    % Get the rational approximation
    N = C*(w.*g(1:m,:));       % Numerator
    D = C*(w.*ones(m,nG));     % Denominator
    R = N./D;
    R(ind(1:m,1),:) = GG(ind(1:m,1),:);            

    % Error in the sample points:
    E = GG-R;
    [maxerr, loc] = compute_norm(E, max_norm);    
    errvec(m) =  maxerr;

    % Check if converged:
    if ( errvec(m) <= tol * normG )
        break
    end
       
end

g = g(1:m,:); w = w(1:m); z = z(1:m); ind = ind(1:m,1);
errvec = errvec(1:m)/normG;

% Remove support points with zero weight:
I = find(w == 0);
z(I) = []; w(I) = []; g(I,:) = []; ind(I) = [];


stats.err_g = errvec(1:m);
stats_T.total = toc(T_total);
stats_T.rest = stats_T.total - stats_T.svd - stats_T.A - stats_T.err_F;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxerr, loc] = compute_norm(E, max_norm)
% E is matrix with each row corresponding to evaluations in z_1, ...
    if max_norm 
        [maxerr, loc] = max(abs(E(:)));             
    else % L2 norm over each sample        
        [maxerr, loc] = max(sum(conj(E).*E,2));     
        maxerr = sqrt(maxerr);        
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [max_norm, tol, svd_update] = parse_options(opts)
% default values
max_norm = true; 
tol = 1e-10;
svd_update = true;

if isfield(opts, 'max_norm')
    max_norm = opts.max_norm;
end
if isfield(opts, 'tol')
    tol = opts.tol;
end
if isfield(opts, 'svd_update')
    svd_update = opts.svd_update;
end
end


