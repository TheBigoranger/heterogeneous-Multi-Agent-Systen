%% load agents data
AgentsData
%% graph
A=[0,0,0,0,0; ...
    1,0,1,1,1; ...
    1,1,0,1,0; ...
    0,1,1,0,1; ...
    0,1,0,1,0];
G = digraph(A);
Gbar =  graph (A(2:end,2:end));
N = length(A) -1;
Lbar = full(Gbar.laplacian);
H = diag(A(2:end,1)) + Lbar;
% figure(3)
% plot(G)
%% design for the estimator
n_0 = length(A_0);
p = height (C1_0);
r_0 = width(B1_0);
P1 = sdpvar (n_0);
lambda_min = min(eig(H));
lambda_max = max(eig(H));
hatK = sdpvar (n_0,p,'full');
gamma1 = sdpvar (1);
lowerbound = sdpvar(1);
constraint1 = [P1*A_0+lambda_min*hatK*C1_0+(P1*A_0+lambda_min*hatK*C1_0)',-P1*B1_0,eye(n_0);
    -B1_0'*P1,-gamma1*eye(r_0),zeros(r_0,n_0);
    eye(n_0),zeros(n_0,r_0),-gamma1*eye(n_0)];
constraint2 = [P1*A_0+lambda_max*hatK*C1_0+(P1*A_0+lambda_max*hatK*C1_0)',-P1*B1_0,eye(n_0);
    -B1_0'*P1,-gamma1*eye(r_0),zeros(r_0,n_0);
    eye(n_0),zeros(n_0,r_0),-gamma1*eye(n_0)];
constraint = [constraint1<=0,constraint2<=0,P1>=lowerbound*eye(n_0);lowerbound>=0];
obj = 10*gamma1 + norm(hatK,'fro') - lowerbound;
% obj = gamma1;
% ops = sdpsettings('solver','bmibnb','bmibnb.relgaptol',1e-4,'verbose',0);
% ops = sdpsettings('solver','lmilab','verbose',1,'debug',0,'shift',1e-12);
ops = sdpsettings('solver','mosek','verbose',0,'debug',0,'shift',1e-12);
% ops = sdpsettings('solver','sdpt3','verbose',1,'debug',0,'shift',1e-12);
% ops = sdpsettings('solver','sedumi','verbose',1,'debug',0,'shift',1e-12);
optimize(constraint,obj,ops);
hatK = value(hatK);
P1 = value(P1);
K = P1^-1*hatK;
%% go through all agents
Ac_i_record = {};
Bc_i_record = {};
Cc_i_record = {};
Dc_i_record = {};
Aceta_i_record = {};
Cceta_i_record = {};
for i = 1:max(size(A_i_data))
    %% follower's dynamics
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
    %% solve for Francis equation
    Pi_i = sdpvar(n_i,n_0,'full');
    Gamma_i = sdpvar(m_i,n_0,'full');
    LME_i = [A_i*Pi_i+B2_i*Gamma_i == Pi_i*A_0; C1_i*Pi_i==C1_0 ];
    ops = sdpsettings('solver','mosek','verbose',0,'debug',0,'shift',1e-12);
    optimize(LME_i,[],ops);
    Pi_i = value(Pi_i);
    Gamma_i = value(Gamma_i);
    %% solve for DOF controller
    Y_i = sdpvar(n_i);
    X_i = sdpvar(n_i);
    Wo_i = sdpvar(n_i,q_i);
    Wc_i = sdpvar(m_i,n_i);
    Dc_i = sdpvar(m_i,q_i);
    J_i = sdpvar(n_i,n_i,'full');
    P_lowerbound =  sdpvar;
    gamma_2_i = sdpvar;

    constraint1 = ...
        [(Y_i*A_i+Wo_i*C2_i)+(Y_i*A_i+Wo_i*C2_i)',        (A_i+B2_i*Dc_i*C2_i+J_i')'       , Y_i*B2_i*Gamma_i-Wo_i*C2_i*Pi_i ,   -Y_i*Pi_i*B1_0   , Y_i*B1_i+Wo_i*D21_i ,       C1_i'      ;
         A_i+B2_i*Dc_i*C2_i+J_i'         , (A_i*X_i+B2_i*Wc_i)+(A_i*X_i+B2_i*Wc_i)', B2_i*Gamma_i-B2_i*Dc_i*C2_i*Pi_i,     -Pi_i*B1_0     , B1_i+B2_i*Dc_i*D21_i,    (C1_i*X_i)'   ;
    (Y_i*B2_i*Gamma_i-Wo_i*C2_i*Pi_i)'   ,   (B2_i*Gamma_i-B2_i*Dc_i*C2_i*Pi_i)'   ,       -gamma_2_i*eye(n_0)       ,   zeros(n_0,r_0)   ,    zeros(n_0,r_i)   ,   zeros(n_0,p)   ;
            (-Y_i*Pi_i*B1_0)'            ,              (-Pi_i*B1_0)'              ,          zeros(r_0,n_0)         , -gamma_2_i*eye(r_0),    zeros(r_0,r_i)   ,   zeros(r_0,p)   ;
          (Y_i*B1_i+Wo_i*D21_i)'         ,         (B1_i+B2_i*Dc_i*D21_i)'         ,          zeros(r_i,n_0)         ,   zeros(r_i,r_0)   , -gamma_2_i*eye(r_i) ,      D11_i'      ;
                   C1_i                  ,                 C1_i*X_i                ,           zeros(p,n_0)          ,    zeros(p,r_0)    ,        D11_i        , -gamma_2_i*eye(p)];
    constraint2  = [Y_i,eye(n_i);
        eye(n_i),X_i];
    constraint = [constraint1<=0,constraint2>=P_lowerbound*eye(2*n_i),...
        P_lowerbound>=0];
    % obj = gamma_2_i;
    obj = gamma_2_i-P_lowerbound;
    % ops = sdpsettings('solver','bmibnb','bmibnb.relgaptol',1e-4,'verbose',0);
    % ops = sdpsettings('solver','lmilab','verbose',0,'debug',0,'shift',1e-12);
    ops = sdpsettings('solver','mosek','verbose',0,'debug',0,'shift',1e-12);
    % ops = sdpsettings('solver','sdpt3','verbose',1,'debug',0,'shift',1e-12);
    % ops = sdpsettings('solver','sedumi','verbose',1,'debug',0,'shift',1e-12);
    optimize(constraint,obj,ops);
    value(gamma_2_i)
    Y_i =value(Y_i );
    X_i =value(X_i );
    Wo_i=value(Wo_i);
    Wc_i=value(Wc_i);
    Dc_i=value(Dc_i);
    J_i =value(J_i );
    S_i = Y_i-X_i^-1;
    Bc_i = S_i^-1 * (Y_i * B2_i *Dc_i-Wo_i);
    Cc_i = Wc_i*X_i^-1-Dc_i*C2_i;
    % Ac_i = S_i^-1*(Y_i*A_i*X_i+Y_i*B2_i*Wc_i-S_i*Bc_i*C2_i*X_i-J_i)*X_i^-1;
    Ac_i = S_i^-1*Y_i*(A_i+B2_i*Dc_i*C2_i+B2_i*Cc_i)-Bc_i*C2_i-S_i^-1*J_i*X_i^-1;
    Cceta_i = Gamma_i-Dc_i*C2_i*Pi_i;
    Aceta_i = -Bc_i*C2_i*Pi_i;
    %% record 
    Ac_i_record = [Ac_i_record, Ac_i];
    Bc_i_record = [Bc_i_record, Bc_i];
    Cc_i_record = [Cc_i_record, Cc_i];
    Dc_i_record = [Dc_i_record, Dc_i];
    Aceta_i_record = [Aceta_i_record, Aceta_i];
    Cceta_i_record = [Cceta_i_record, Cceta_i];
end

%% save to the workspace
save('controlGain.mat', 'K','Ac_i_record','Bc_i_record','Aceta_i_record',...
    'Cc_i_record','Dc_i_record','Cceta_i_record');