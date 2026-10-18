function [pib_SDP, pifrctn_SDP, pi_full,sol2] = PhysicallyConsistentEstimation(Tau, Wb, baseQR)
% -------------------------------------------------------------------------
% 在物理一致性的约束下，对待辨识参数进行辨识
% 约束条件为半正定
% -------------------------------------------------------------------------
physicalConsistency = 1;
pi_frctn = sdpvar(4,1);                             % 摩擦参数
pi_b = sdpvar(baseQR.numberOfBaseParameters, 1);    % 基础参数
numberOfDependentParameters = size(baseQR.permutationMatrix, 1) - baseQR.numberOfBaseParameters;
pi_d = sdpvar(numberOfDependentParameters, 1);    % 线性相关参数
pii = baseQR.permutationMatrix*[eye(baseQR.numberOfBaseParameters), -baseQR.beta;
                                zeros(numberOfDependentParameters, baseQR.numberOfBaseParameters), eye(numberOfDependentParameters)]*[pi_b; pi_d];
 
    % 约束
    mass_indexes = 10:10:20;
    cnstr = [];

    % 质量
    for i = 1:2
        cnstr = [cnstr, pii(mass_indexes(i)) >= 0];    
    end

    % 惯量
    if physicalConsistency
        for i = 1:10:20
            link_inertia_i = [pii(i), pii(i+1), pii(i+2); ...
                              pii(i+1), pii(i+3), pii(i+4); ...
                              pii(i+2), pii(i+4), pii(i+5)];          
            frst_mmnt_i = pii(i+6:i+8);
        
            Di = [0.5*trace(link_inertia_i)*eye(3) - link_inertia_i, frst_mmnt_i;
                                                       frst_mmnt_i',    pii(i+9)];
        
            cnstr = [cnstr, Di>=0, pii(i)>=0, pii(i+3)>=0, pii(i+5)>=0];
        end
    else
        for i = 1:10:20
            link_inertia_i = [pii(i), pii(i+1), pii(i+2); ...
                              pii(i+1), pii(i+3), pii(i+4); ...
                              pii(i+2), pii(i+4), pii(i+5)];
        
            frst_mmnt_i = vecS(pii(i+6:i+8));
        
            Di = [link_inertia_i, frst_mmnt_i'; frst_mmnt_i, pii(i+9)*eye(3)];
            cnstr = [cnstr, Di>=0, pii(i)>=0, pii(i+3)>=0, pii(i+5)>=0,pii(i+10)>=0];
        end
    end

    % 对摩擦参数的约束
    for i = 1:2
        cnstr = [cnstr, pi_frctn(2*i-1)>=0, pi_frctn(2*i)>=0];  
    end

    % 目标函数
    obj = norm(Tau - Wb*[pi_b;pi_frctn]);

    sol2 = optimize(cnstr, obj, sdpsettings('solver','sdpt3'));

    pib_SDP = value(pi_b);   % 基础惯性参数变量 
    pifrctn_SDP = value(pi_frctn);

    pi_full = baseQR.permutationMatrix*[eye(baseQR.numberOfBaseParameters), -baseQR.beta; ...
                                        zeros(numberOfDependentParameters,baseQR.numberOfBaseParameters), eye(numberOfDependentParameters)]*[value(pi_b); value(pi_d)];

end






function m = vecS(v) 
% -------------------------------------------------------------------------
% 用于叉积 
% -------------------------------------------------------------------------
    m = [0, -v(3), v(2);
        v(3), 0, -v(1);
        -v(2), v(1) 0]; 
end
     






