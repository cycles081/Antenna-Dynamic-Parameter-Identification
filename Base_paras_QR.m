function [pi_lgr_base, baseQR] = Base_paras_QR(seed)
% -------------------------------------------------------------------------
% 通过QR分解计算基础参数集；当前20列回归矩阵未单独建模电机转动惯量或摩擦
% 
% 
% -------------------------------------------------------------------------
%% --------------------- 符号参数 -----------------------------------------
% syms q1 q2 dq1 dq2 ddq1 ddq2 g  L1 L2 real 
m = sym('m',[2 1],'real');
hx = sym('hx',[2 1],'real');
hy = sym('hy',[2 1],'real');
hz = sym('hz',[2 1],'real');
Ixx = sym('Ixx',[2 1],'real');
Ixy = sym('Ixy',[2 1],'real');
Ixz = sym('Ixz',[2 1],'real');
Iyz = sym('Iyz',[2 1],'real');
Iyy = sym('Iyy',[2 1],'real');
Izz = sym('Izz',[2 1],'real');
J = sym('J',[1 2],'real');

Pi_1=[Ixx(1),Ixy(1),Ixz(1),Iyy(1),Iyz(1),Izz(1),m(1)*hx(1),m(1)*hy(1),m(1)*hz(1),m(1)]'; 
Pi_2=[Ixx(2),Ixy(2),Ixz(2),Iyy(2),Iyz(2),Izz(2),m(2)*hx(2),m(2)*hy(2),m(2)*hz(2),m(2)]';
pi_lgr_sym = [Pi_1, Pi_2];

[nLnkPrms, nLnks] = size(pi_lgr_sym);

pi_lgr_sym = reshape(pi_lgr_sym, [nLnkPrms*nLnks, 1]);

%% --------------------- 随机生成观测矩阵 ----------------------------------
if nargin < 1
    seed = 0;
end
rng(seed);                    % 固定种子，便于复现 QR 参数选择

q_min = 0;      % 运动范围
q_max = pi*ones(2,1);
qd_max = pi/18*ones(2,1);
q2d_max = pi/18*ones(2,1);
% n = 30;
% q1set = -1+2*rand(n, 1); q2set = -1+2*rand(n, 1);
% dq1set = -1+2*rand(n, 1); dq2set = -1+2*rand(n, 1);
% ddq1set = -1+2*rand(n, 1); ddq2set = -1+2*rand(n, 1);

W = [];
for i = 1:50
%      q = [q1set(i), q2set(i)];
%     dq = [dq1set(i), dq2set(i)];
%     ddq = [ddq1set(i), ddq2set(i)];
    q_rnd = q_min + (q_max - q_min).*rand(2,1);
    qd_rnd = -qd_max + 2*qd_max.*rand(2,1);
    q2d_rnd = -q2d_max + 2*q2d_max.*rand(2,1);

    Y = RegressorWithMotor(q_rnd, qd_rnd, q2d_rnd);
    W = vertcat(W,Y);
end

%% --------------------- QR分解 -------------------------------------------
[~, R, E] = qr(W);
bb = rank(W);             % 观测矩阵的秩
R1 = R(1:bb,1:bb);
R2 = R(1:bb,bb+1:end);
beta = R1\R2;             % 相关系数
beta(abs(beta)<sqrt(eps)) = 0; 

W1 = W*E(:,1:bb);
W2 = W*E(:,bb+1:end);
assert(norm(W2 - W1*beta) < 1e-6,...
    'Found realationship between W1 and W2 is not correct\n')

%% --------------------- 基础参数集 ---------------------------------------
pi1 = E(:,1:bb)'*pi_lgr_sym;         % 线性无关参数
pi2 = E(:,bb+1:end)'*pi_lgr_sym;     % 线性相关和全零列对应的参数

pi_lgr_base = pi1 + beta*pi2;        % 重构后的基础参数集

%% --------------------- 最终结果 -----------------------------------------
baseQR = struct;
baseQR.numberOfBaseParameters = bb;
baseQR.permutationMatrix = E;
baseQR.beta = beta;

% end




