function identraj = Filtered_data(data,kaishi,jieshu,fc,fs,r,q0)
% --------------------------------------------------------------
% 根据提取到的数据滤波并返回相应的数据
% inputs:
%     data - 原始数据集
%     kaishi - 截取数列左端点
%     jieshu - 截取数据右端点
%     fc - 截止频率
%     fs - 采样频率
%     q0 - 初始角位移
% outputs:
%     identraj - struct
% --------------------------------------------------------------
t = data(:,1)/fs;
step = t(2) - t(1);
time = t(kaishi:jieshu);

fangwei_pos = data(kaishi:jieshu,5);               % 方位关节位置数据
fangwei_vel = diff(fangwei_pos)/step;
fangwei_vel = [0;fangwei_vel];
fangwei_acc = diff(fangwei_vel)/step;
fangwei_acc = [0;fangwei_acc];

fuyang_pos = 85 - data(kaishi:jieshu,2);           % 俯仰关节位置数据
fuyang_vel = diff(fuyang_pos)/step;
fuyang_vel = [0;fuyang_vel];
fuyang_acc = diff(fuyang_vel)/step;
fuyang_acc = [0;fuyang_acc];

fangwei_tau = data(kaishi:jieshu,7)*0.147/1000;    % 方位关节力矩
fuyang_tau = -data(kaishi:jieshu,4)*0.147/1000;    % 俯仰关节力矩

% mean_fuyangtau = mean(fuyang_tau(1:2*T*fs));       % 俯仰力矩均值
mean_fuyangtau = 0;

[b,a] = butter(6, fc/(fs/2));                      % 滤波器

fangwei_pos = filtfilt(b,a,fangwei_pos - q0(1));
fangwei_vel = filtfilt(b,a,fangwei_vel);
fangwei_acc = filtfilt(b,a,fangwei_acc);
fangwei_tau = filtfilt(b,a,fangwei_tau);

fuyang_pos = filtfilt(b,a,fuyang_pos - q0(2));
fuyang_vel = filtfilt(b,a,fuyang_vel);
fuyang_acc = filtfilt(b,a,fuyang_acc);
fuyang_tau = filtfilt(b,a,fuyang_tau - mean_fuyangtau);

q = [fangwei_pos + q0(1), fuyang_pos + q0(2)];
dq = [fangwei_vel, fuyang_vel];
ddq = [fangwei_acc, fuyang_acc];
tau = [fangwei_tau, fuyang_tau + mean_fuyangtau];

figure
subplot(2,1,1)
plot(time,q(:,1),'--r','LineWidth',1)
hold on
plot(time,dq(:,1),'--b','LineWidth',1)
hold on 
plot(time,ddq(:,1),'--m','LineWidth',1)
hold on
plot(time,r(1)*tau(:,1),'--g','LineWidth',1)
grid on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.25);
set(gca,'LineWidth',1.25,'FontName','Times New Roman','FontSize',12);
xlabel('','FontName','Times New Roman','FontSize',12);
ylabel('','FontName','Times New Roman','FontSize',12);
legend('方位位置','方位速度','方位加速度','方位力矩','FontName','Songti');

subplot(2,1,2)
plot(time,q(:,2),'--r','LineWidth',1)
hold on
plot(time,dq(:,2),'--b','LineWidth',1)
hold on 
plot(time,ddq(:,2),'--m','LineWidth',1)
hold on
plot(time,r(2)*tau(:,2),'--g','LineWidth',1)
grid on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.25);
set(gca,'LineWidth',1.25,'FontName','Times New Roman','FontSize',18);
xlabel('','FontName','Times New Roman','FontSize',18);
ylabel('','FontName','Times New Roman','FontSize',18);
legend('俯仰位置','俯仰速度','俯仰加速度','俯仰力矩','FontName','Songti');

identraj = struct;
identraj.t = t(kaishi:jieshu);
identraj.q = deg2rad(q);
identraj.dq = deg2rad(dq);
identraj.ddq = deg2rad(ddq);
identraj.tau = tau;

figure
plotCount = min(4000, numel(time));
plot(time(1:plotCount),r(1)*tau(1:plotCount,1),'-b','LineWidth',1.5)
hold on 
plot(time(1:plotCount),r(2)*tau(1:plotCount,2),'-r','LineWidth',1.5)
grid on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.25);
set(gca,'LineWidth',1.25,'FontName','Times New Roman','FontSize',12);
xlabel('time','FontName','Times New Roman','FontSize',12);
ylabel('torque','FontName','Times New Roman','FontSize',12);
legend('方位力矩','俯仰力矩','FontName','Songti');

end

