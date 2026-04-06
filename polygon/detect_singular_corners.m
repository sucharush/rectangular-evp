function corner_list = detect_singular_corners(V)
% V is assumed closed: V(end,:) = V(1,:)

    corners = V(1:end-1,:);
    nC = size(corners,1);

    corner_list = {};
    tol = 1e-10;

    for i = 1:nC
        c = corners(i,:);
        p = corners(mod(i-2,nC)+1,:);   % previous
        q = corners(mod(i,nC)+1,:);     % next

        a = p - c;
        b = q - c;

        a = a / norm(a);
        b = b / norm(b);

        detBA = b(1)*a(2) - b(2)*a(1);
        dotBA = b(1)*a(1) + b(2)*a(2);

        omega = atan2(detBA, dotBA);
        if omega <= 0
            omega = omega + 2*pi;
        end

        alpha = pi / omega;
        is_regular = abs(alpha - round(alpha)) < tol;

        if ~is_regular
            s = struct();
            s.c = c;
            s.e0 = b;
            s.omega = omega;
            s.alpha = alpha;
            corner_list{end+1} = s; %#ok<AGROW>
        end
    end

    if isempty(corner_list)
        error('No singular corners detected.');
    end
end