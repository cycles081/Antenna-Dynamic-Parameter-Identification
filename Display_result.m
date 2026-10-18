function [X,error_fangwei,error_fuyang] = Display_result(sol,identraj,Wb,gearRatio)
% -------------------------------------------------------------------
% 根据辨识结果及采集到的数据，计算辨识误差
% inputs:
%     sol - 
%     identraj - 
% outputs:
%     X - 
%     error_fangwei - 
%     error_fuyang - 
% 2023.12.16
% -------------------------------------------------------------------
X = [sol.pi_b;sol.pi_fr];
if nargin < 4
    gearRatio = [560;700];
end
tau2 = reshape(Wb*X, 2, []);
fangwei_tau2 = tau2(1,:);
fuyang_tau2 = tau2(2,:);
error_fangwei = fangwei_tau2' - gearRatio(1)*identraj.tau(:,1);
error_fuyang =  fuyang_tau2' - gearRatio(2)*identraj.tau(:,2);
t = identraj.t;

% figure
% subplot(2,1,1)
tu = tiledlayout(2,1,'TileSpacing','compact');
nexttile
plot(t,gearRatio(1)*identraj.tau(:,1), '-', 'LineWidth',1)
hold on 
plot(t,fangwei_tau2', '--', 'LineWidth',1)
hold on
plot(t,error_fangwei,'-.', 'LineWidth',1)
grid on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.25);
set(gca,'LineWidth',1.25,'FontName','Times New Roman','FontSize',12);
% xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
% ylabel('Torque(N·m)','FontName','Times New Roman','FontSize',12);

% subplot(2,1,2)
nexttile
plot(t,gearRatio(2)*identraj.tau(:,2), '-', 'LineWidth',1)
hold on 
plot(t,fuyang_tau2', '--', 'LineWidth',1)
hold on 
plot(t,error_fuyang, '-.', 'LineWidth',1)
grid on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.25);
set(gca,'LineWidth',1.25,'FontName','Times New Roman','FontSize',12);
% xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
% ylabel('Torque(N·m)','FontName','Times New Roman','FontSize',12);

lgd = legend({'测试'; '仿真'; '误差'}, 'FontName', 'Songti', 'Location', 'east', 'Orientation', 'horizontal');
lgd.Layout.Tile = 'north';
xlabel(tu, 'Time(s)','FontName','Times New Roman','FontSize',12)
ylabel(tu, 'Torque(N·m)','FontName','Times New Roman','FontSize',12)

end

