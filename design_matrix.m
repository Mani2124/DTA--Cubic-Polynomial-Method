function A = design_matrix(nodes, t)
    interval = 0.5;
    x = zeros(length(nodes) * 2, 1);
    A = zeros(length(t), length(x));

    n = 1;
    k = 1;
    for i = 1:length(nodes) - 1
        for j = 1:length(t)
            if (nodes(i) <= t(j)) && (nodes(i+1) > t(j))
                values = (t(j) - nodes(i)) / interval;
                A(n, k) = 1 - 3 * values^2 + 2 * values^3;
                A(n, k+1) = (t(j) - nodes(i)) * (1 - 2 * values + values^2);
                A(n, k+2) = 3 * values^2 - 2 * values^3;
                A(n, k+3) = (t(j) - nodes(i)) * (-values + values^2);
                n = n + 1;
            end
        end
        k = k + 2;
    end
end
