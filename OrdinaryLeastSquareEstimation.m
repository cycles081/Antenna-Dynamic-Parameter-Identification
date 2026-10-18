function [pib_OLS, pif_OLS] = OrdinaryLeastSquareEstimation(Tau, Wb, numberOfBaseParameters)
% -------------------------------------------------------------------------
% 最小二乘法求解
% -------------------------------------------------------------------------
pi_OLS = Wb\Tau;
pib_OLS = pi_OLS(1:numberOfBaseParameters);
pif_OLS = pi_OLS(numberOfBaseParameters+1:end);
end
