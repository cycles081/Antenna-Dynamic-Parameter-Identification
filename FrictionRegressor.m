function Yf = FrictionRegressor(dq)
% -------------------------------------------------------------------------
% 根据速度生成摩擦模型的回归矩阵
% 采用线性摩擦模型
% tau_f = Fv*qd + Fc*sign(qd) + F0
% dq - 1x2 vector
% -------------------------------------------------------------------------
noJoints = size(dq,2);
Yf = zeros(noJoints, noJoints*2);
for i = 1:noJoints     
    Yf(i,2*i-1:2*i) = [dq(i), tanh(50*dq(i))]; 
end
