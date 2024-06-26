function stderr_cluster = se_cluster(X, eps, coefinds, clusterids)
    Xeps = eps .* X;
    [cl, co] = ndgrid(clusterids, 1:size(Xeps,2));
    Cu = accumarray([cl(:) co(:)], Xeps(:));
    V_cluster = (X' * X)\(Cu' * Cu)/(X' * X);
    stderr_cluster = diag(V_cluster);
    stderr_cluster = sqrt(stderr_cluster(coefinds));
end