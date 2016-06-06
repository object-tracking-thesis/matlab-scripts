function hasGating = ellipGating(ell_mu, ell_cov, meas_mu, Pg)
    % ell_mu    - the ellipse center
    % ell_cov   - the ellipse covariance
    % meas_mu   - measurement cluster center to check
    % Pg        - desired confidence intervall in %
    %
    r = length(meas_mu);
    threshold = chi2inv(Pg, r);

    d = meas_mu - ell_mu(1:r);
    cal = d'*inv(ell_cov)*d;

%     if (cal < threshold) & (abs(d) < [3;3])
    if (abs(d) < [2;2]) & (abs(ell_mu(5)) < 5)
        hasGating = 1;
    else
        hasGating = 0;
    end
end