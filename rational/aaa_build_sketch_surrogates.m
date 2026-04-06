function GZ = aaa_build_sketch_surrogates(F_op, Z, U, V)
% GZ(i,j) = U(:,i)^* F_op(Z(j)) V(:,i)
% size(GZ) = ell x numel(Z)

    ell = size(U, 2);
    M = numel(Z);
    GZ = zeros(ell, M);

    for j = 1:M
        Fj = F_op(Z(j));
        for i = 1:ell
            GZ(i,j) = U(:,i)' * Fj * V(:,i);
        end
    end
end
