function Y = RegressorWithoutMotor(q, dq, ddq)
% -------------------------------------------------------------------------
% 牛顿欧拉法建立方位及俯仰关节的动力学方程和观测矩阵,包括摩擦
% Inputs:
%     q - 关节角度, 1x2 向量;
%     dq - 关节角速度, 1x2 向量;
%     ddq - 关节角加速度, 1x2 向量;
% Outputs:
%     Y - 回归矩阵，2x20 矩阵。
% 2023.10.21
% -------------------------------------------------------------------------
%% ------------------------------------------------------------------------
q1 = q(1); q2 = q(2);
dq1 = dq(1); dq2 = dq(2);
ddq1 = ddq(1); ddq2 = ddq(2);
L1 = 0.132;
g = 9.8;

%% ---------------------坐标变换矩阵---------------------------------------
T01 = [cos(q1), -sin(q1), 0, 0;
       sin(q1), cos(q1),  0, 0;
       0,       0,        1, 0;
       0,       0,        0, 1];
T12 = [cos(q2), -sin(q2), 0, L1;
       0,      0,        -1,  0;
       sin(q2), cos(q2), 0,  0 ;
       0,      0,        0,    1];
% T02 = T01*T12;
R01=T01(1:3,1:3);P01=T01(1:3,4);
R12=T12(1:3,1:3);P12=T12(1:3,4);


%% --------------------计算方位以及俯仰环节加速度-------------------------
% 方位及俯仰关节旋转方向
axis1 = [0 0 1]';
axis2 = [0 0 1]';

% % 机架及反射面质心位置
% rc1 = [hx(1) hy(1) hz(1)]';
% rc2 = [hx(2) hy(2) hz(2)]';

% 方位及俯仰关节角速度
omega_base=[0 0 0]';
omega0=omega_base;
omega1=R01'*omega0+dq1*axis1; 
omega2=R12'*omega1+dq2*axis2;

% 方位及俯仰关节速度
vel_base=[0 0 0]';
vel0=vel_base;
vel1=R01'*(vel0+getS(omega0)*P01);
vel2=R12'*(vel1+getS(omega1)*P12);

% 方位及俯仰关节角加速度
ang_acc_base=[0 0 0]';
ang_acc0=ang_acc_base;
ang_acc1=R01'*ang_acc0+R01'*getS(omega1)*axis1*dq1+axis1*ddq1;
ang_acc2=R12'*ang_acc1+R12'*getS(omega2)*axis2*dq2+axis2*ddq2;

% 方位及俯仰关节加速度
acc_base=[0 0 g]';
acc0=acc_base;
acc1=R01'*(acc0+getS(ang_acc0)*P01+getS(omega0)*getS(omega0)*P01);
acc2=R12'*(acc1+getS(ang_acc1)*P12+getS(omega1)*getS(omega1)*P12);
% % 机架及反射面质心加速度
% ac_base = [0 0 0]';
% ac1 = getS(ang_acc1)*rc1 + getS(omega1)*(getS(omega1)*rc1) + acc1;
% ac2 = getS(ang_acc2)*rc2 + getS(omega2)*(getS(omega2)*rc2) + acc2;
% ac1 = simplify(ac1); ac2 = simplify(ac2);

%% ---------------------线性化动力学方程-----------------------------------
A1 = getAi(ang_acc1,omega1,acc1);
A2 = getAi(ang_acc2,omega2,acc2);

Q12 = getTi(R12, P12);

U22 = A2;
U11 = A1; U12 = Q12*A2;

Y22 = [0 0 0 0 0 1]*U22;
Y11 = [0 0 0 0 0 1]*U11; Y12 = [0 0 0 0 0 1]*U12;

% Y = [Y11,        Y12;
%     zeros(1,10), Y22];              % Regressor Matrix, 2x20

%% ---------------------摩擦模型-------------------------------------------
%linear friction model:  tau_fri = f1*dq + f2*sign(dq) + f3
% Y = [Y11,        Y12, dq(1), sign(dq(1)), 1,     0,           0, 0;
%     zeros(1,10), Y22,     0,           0, 0, dq(2), sign(dq(2)), 1];

Y = [        Y11,  Y12;
     zeros(1,10),  Y22];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------------------主函数结束----------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----------------------旋量矩阵算子--------------------------------------  
function screw=getS(vec)
[row,col]=size(vec);
if row~=3 && col~=1
    error("worng,input of getS must be 3*1 vector");
end
screw=[0       -vec(3)  vec(2);
       vec(3)   0      -vec(1);
      -vec(2)   vec(1)  0     ];
end

%% ----------------------速度矩阵算子K-------------------------------------
function K=getK(vec3)
K=[vec3(1) vec3(2) vec3(3) 0       0       0;
    0      vec3(1) 0       vec3(2) vec3(3) 0;
    0      0       vec3(1) 0       vec3(2) vec3(3)];
end

%% ----------------------变换算子T-----------------------------------------
function Ti = getTi(R, P)
Ti(1:3,1:3) = R;
Ti(4:6,1:3) = getS(P)*R;
Ti(4:6,4:6) = R;
end

%% ----------------------力矩简化算子A-------------------------------------
function Ai = getAi(dotwi,wi,ai)
% Ai(1:3,1:6) = zeros(3,6);
Ai(1:3,7:9) = getS(dotwi) + getS(wi)*getS(wi);
Ai(1:3,10)  = ai;
Ai(4:6,1:6) = getK(dotwi)+getS(wi)*getK(wi);
Ai(4:6,7:9) = -getS(ai);
% Ai(4:6,10)  = zeros(3,1);
end

