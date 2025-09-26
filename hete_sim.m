%% design the controller
hete_synthesis
% AgentsData
% load("controlGain.mat")
%% simulate the leader estimator
n_0 = length(A_0);
p = height (C1_0);
r_0 = width(B1_0);
t = 0:0.001:40;
omega= zeros(1,length(t));
omega(1:floor(0.1*length(t)))=1;
omega(floor(0.5*length(t)):floor(0.7*length(t)))=-1;
estimateSys = ss(kron(eye(N+1),A_0)+kron(diag(sum(A,2))-A,K*C1_0),[B1_0;zeros(N*n_0,r_0)],kron(eye(N+1),eye(n_0)),zeros((N+1)*n_0,r_0));
eta_and_x0 = lsim(estimateSys,omega,t);
figure(1)
clf
for i=1:n_0
    subplot(n_0,1,i);
    plot(t,eta_and_x0(:,i:n_0:end));
end
ylim padded
%% simulation for DOF
figure(2)
clf
% plot leader's z_0
z_0=C1_0*eta_and_x0(:,1:n_0)';
for j=1:n_0
    subplot(n_0,1,j);
    plot(t,z_0(j,:));
    hold on
end
% simulate for followers
for i = 1:max(size(A_i_data))
    A_i=A_i_data{i};
    B1_i =B1_i_data {i};
    B2_i =B2_i_data {i};
    C1_i =C1_i_data {i};
    C2_i =C2_i_data {i};
    D11_i=D11_i_data{i};
    D21_i=D21_i_data{i};
    % dimensions
    n_i = length(A_i);
    r_i = width(B1_i);
    m_i = width(B2_i);
    q_i = height(C2_i);
    Ac_i=Ac_i_record{i};
    Bc_i=Bc_i_record{i};
    Cc_i=Cc_i_record{i};
    Dc_i=Dc_i_record{i};
    Aceta_i=Aceta_i_record{i};
    Cceta_i=Cceta_i_record{i};
    DOFsys = ss([A_i+B2_i*Dc_i*C2_i,B2_i*Cc_i; ...
        Bc_i*C2_i,Ac_i],[B2_i*Cceta_i,B1_i+B2_i*Dc_i*D21_i; ...
        Aceta_i,Bc_i*D21_i],[C1_i,zeros(p,n_i)],[zeros(p,n_0),D11_i]);
    % xi_i = exp(-1*t).*randn(r_i,length(t));
    xi_i = exp(-1*t).*sin(i*length(t)).*ones(r_i,length(t));
    % xi_i = zeros(r_i,length(t));
    [z_i, tOut, xandxc_i]=lsim(DOFsys,[eta_and_x0(:,i*n_0+1:(i+1)*n_0),xi_i'],t);
    for j=1:p
        subplot(p,1,j);
        plot(tOut,z_i(:,j:p:end));
        hold on
    end
end
% figure marks
ylim padded
subplot(n_0,1,1);
legend('leader','Agent 1','Agent 2','Agent 3','Agent 4', ...
        'Location','northoutside','NumColumns',N+1);
for i= 1:p
    subplot(p,1,i);
    ylabel(['$z_' num2str(i) '^{(i)}$'],...
        'interpreter','latex', 'FontSize',15)
    xlabel('$t$','interpreter','latex','FontSize',15)
end