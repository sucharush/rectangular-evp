function [vec_tls_all, lam_tls_all] = tls_pencil_eigs(A, B)
%TLS_PENCIL_EIGS  TLS-type eigenvalue extraction for a rectangular pencil A - lambda B.
%
% INPUT
%   A, B : n-by-r matrices
%
% OUTPUT
%   lam_tls_all : r-by-1 vector of generalized eigenvalues
%
% Method:
%   Form C = [B, A], compute its SVD, partition the first r right singular
%   vectors as
%       Vc(:,1:r) = [V11; V21],   V11, V21 are r-by-r,
%   then solve
%       V21' * x = lambda * V11' * x.

    [nA, rA] = size(A);
    [nB, rB] = size(B);

    if nA ~= nB || rA ~= rB
        error('A and B must have the same size n-by-r.');
    end

    r = rA;

    C = [B, A];              % n x 2r
    [~, Sigma, Vc] = svd(C, 0);  % economy SVD

    if size(Vc, 2) < r
        error('SVD did not return enough right singular vectors.');
    end

    V11 = Vc(1:r,     1:r);
    V21 = Vc(r+1:2*r, 1:r);

    [vec_tls_all, lam_tls_all] = eig(V21', V11', 'vector');
    % lam_tls_all = diag(D);
end