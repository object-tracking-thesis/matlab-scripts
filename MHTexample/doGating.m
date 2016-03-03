function okMeasurements = doGating(recMeas,predMeas, Sk, Pg)
% Performs gating of measurements according to elliptical gating. Only
% measurements which are within the threshold are returned. The threshold
% is gotten from chi2inv(Pg, n), where Pg is the confidence interval and n is
% the dimension of a measurement. Gating initself is performed according to
% (y_{k} - y_{k|k-1})'*inv(S_{k|k-1})*(y_{k} - y_{k|k-1}) < threshold 
% 
% okMeasurements = doGating(recMeas,predMeas, Sk, Pg)
%
%     recMeas  - Received Measurements. Corresponds to y_{k}
%     predMeas - The predicted measurement. Corresponds to y_{k|k-1}
%     Sk       - The innovation gain matrix, gotten from KF 
%     Pg       - The confidence interval to be used for the gate
%
    [r, c] = size(recMeas); % Rows = state, columns = measurement
    threshold = chi2inv(Pg, r);
    okMeasurements = [];

    for j = 1:c
        d = recMeas(:,j) - predMeas;
        cal = d'*inv(Sk)*d;

        if cal < threshold
            okMeasurements = [okMeasurements recMeas(:,j)];
        end
    end

end