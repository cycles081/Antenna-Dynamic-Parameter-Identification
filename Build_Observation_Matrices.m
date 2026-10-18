function [Tau, Wb] = Build_Observation_Matrices(identraj,baseQR, a)
% -----------------------------------------------------------------------
% 根据带有电机惯量的动力学回归矩阵和摩擦模型的回归矩阵构造观测矩阵
% inputs:
%     identraj - 用于辨识的结构体，包括角位移、角速度、角加速度、力矩和时间
%     baseQR - 通过对动力学回归矩阵QR分解得到的基础参数集、基础参数集数目和对应的置换矩阵
%     a - 两个关节的减速比
% outputs:
%     Tau - 力矩
%     Wb - 观测矩阵
% 2023.12.07
% -----------------------------------------------------------------------
q = identraj.q;
dq = identraj.dq;
ddq = identraj.ddq;
tau = identraj.tau;
T = identraj.t;

Wb = [];
Tau = [];
E1 = baseQR.permutationMatrix(:,1:baseQR.numberOfBaseParameters);

    for i = 1:length(T)
        Yi = RegressorWithoutMotor(q(i,:), dq(i,:), ddq(i,:));
        Yi_frc = FrictionRegressor(dq(i,:));
        Ybi = [Yi*E1, Yi_frc];
        Wb = vertcat(Wb, Ybi);
    
        Tau = vertcat(Tau, diag(a)*tau(i,:)');
    end

end




