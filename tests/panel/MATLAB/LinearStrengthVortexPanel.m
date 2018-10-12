% Alwin Wang
%% Set Up
clear
clc
%% Airfoil
% Sample diamond airfoil
tc = 0.1; p = 0.3;
x = [1   p   0  p   1];
y = [0 -tc/2 0 tc/2 0];
figure()
plot(x, y)
% Matrices describing the airfoil
PT = [x' y'];
DX = x(2:end) - x(1:(end-1));
DY = y(2:end) - y(1:(end-1));
CO = [x(1:(end-1))'+DX'/2, y(1:(end-1))'+DY'/2];
TH = atan2(DY, DX)';
DL = sqrt(DX.^2 + DY.^2)';
M  = length(TH);
N  = M + 1;
%% Influence Matrix
% Preallocate for speed
A = zeros(M, N);
B = zeros(M, N);
alpha = 5;
AL    = alpha*pi/180; 

for i = 1:M % For each collocation point
    %% Convert to Local Panel Coord System
    XG = CO(i,1)   - PT(1:M,1);
    YG = CO(i,2)   - PT(1:M,2);
    X2 = PT(2:N,1) - PT(1:M,1);
    Y2 = PT(2:N,2) - PT(1:M,2);
    X =  XG.*cos(TH) + YG.*sin(TH);
    Y = -XG.*sin(TH) + YG.*cos(TH);
    X1 = zeros(M,1);
    Y1 = zeros(M,1);
    X2 = X2.*cos(TH) + Y2.*sin(TH);
    Y2 = zeros(M,1);
    % Global = [PT(1:M,1:2) PT(2:N,1:2) ones(M,1)*CO(i,1:2)];
    % Local  = [X1 Y1 X2 Y2 X Y];
    %% Determine Induced Velocity Components
    % Collocation to panel
    R1  = sqrt((X-X1).^2 + (Y-Y1).^2);
    R2  = sqrt((X-X2).^2 + (Y-Y2).^2);
    TH1 = atan2(Y-Y1,X-X1);
    TH2 = atan2(Y-Y2,X-X2);
    % ColPan = [R1 R2 TH1 TH2];
    UAL = -(Y.*log(R2./R1) + X.*(TH2-TH1)-X2.*(TH2-TH1))         ./(2*pi*X2);
    UBL =  (Y.*log(R2./R1) + X.*(TH2-TH1))                       ./(2*pi*X2);
    WAL = -((X2-Y.*(TH2-TH1)) - X.*log(R1./R2) + X2.*log(R1./R2))./(2*pi*X2);
    WBL =  ((X2-Y.*(TH2-TH1)) - X.*log(R1./R2))                  ./(2*pi*X2);
    % UWLocal = [UAL WAL UBL WBL];
    %% Convert to Global Coord System
    UA =  UAL.*cos(-TH) + WAL.*sin(-TH);
    UB =  UBL.*cos(-TH) + WBL.*sin(-TH);
    WA = -UAL.*sin(-TH) + WAL.*cos(-TH);
    WB = -UBL.*sin(-TH) + WBL.*cos(-TH);
    % UWGlobal = [UA WA UB WB];
    %% Normal and Tangential Velocities
    NorA = -UA.*sin(TH(i)) + WA.*cos(TH(i));
    NorB = -UB.*sin(TH(i)) + WB.*cos(TH(i));
    TanA =  UA.*cos(TH(i)) + WA.*sin(TH(i));
    TanB =  UB.*cos(TH(i)) + WB.*sin(TH(i));
    %% Add to Matrix System
    A(i,:) = [NorA' 0] + [0 NorB'];
    B(i,:) = [TanA' 0] + [0 NorB'];    
end

%% Add Kutta Condition and BC
A(N,:)   = [1 zeros(1,N-2) 1];
A(:,N+1) = [cos(AL)*sin(TH)-sin(AL)*cos(TH); 0];

%% Solve
R = rref(A);
G = R(:,N+1);

%% Variables of Interest
% Assumed there was a unit velocity
V  = B*G - (cos(AL)*cos(TH) + sin(AL)*sin(TH)); % I think this line is wrong...
CP = 1-V.^2;
CL = -1*CP.*(cos(AL)*cos(TH) + sin(AL)*sin(TH)).*DL;