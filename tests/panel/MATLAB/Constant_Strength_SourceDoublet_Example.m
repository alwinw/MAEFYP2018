%% Constant strength source-doublet example:

clear
close all

%Half the mesh resolution of each element:
Npanel = 200;

%Define the free stream speed:
u_infty = 1;

alpha = 2*pi/180;
%alpha = 5*pi/180;

%Aerofoil of main element (must be expressed in XXXX form):
Aerofoil = 4521;
%Aerofoil_main = 0015;

%maximum camber of main element:
m = 1/100*floor(Aerofoil/1000);

%location of the maximum camber of main element:
p = floor((Aerofoil-100*m*1000)/100)/10;

%chord of main element:
chord = 1;

%associated thickness:
thickness = chord*(Aerofoil - p*10*100 - ...
    m*100*1000)/100;

%Define the panel geometry of the main element:
[x,z,x_c,z_c,n_hat,t_hat,sin_theta,cos_theta,X2] = BodyPanelDefinition(...
    chord,thickness,Npanel,m,p);

%Move the collocation points inwards:
CP = 0;

%Make sure the relocation of the collocation point ensures that the angles
%are preserved:
TestCp = 0;

if CP == 1

    %Internalise the collocation points:
    [x_c,z_c] = CollocationInternalise(X2,sin_theta,cos_theta,x_c,z_c,...
        TestCp);
    
end

%Influence coefficients (Doublet):
[a] = DoubletInfluenceCoeffDirichlet(cos_theta,sin_theta,x_c,z_c,X2,x,z,CP);

%Influence coefficients (Source):
[b] = SourceInfluenceCoeffDirichlet(cos_theta,sin_theta,x_c,...
     z_c,X2,x,z,CP); 

%RHS (equations 11.61 and 11.72):
[sigma] = SourceStrengths(n_hat,u_infty,alpha);
RHS = -b*sigma;

%Determine the strengths of the doublets:
mu = a\RHS;

Q_t = zeros(length(mu)-1,1);
zeta = zeros(length(mu)-1,1);
Cp = zeros(length(mu)-1,1);

%Equation 11.76:
for i = 1:length(mu)
    
    zeta(i) = u_infty*(x_c(i)*cos(alpha) + z_c(i)*sin(alpha)) + mu(i);
    
end

for i = 1:length(mu)-1
        
    Q_t(i) = 2*(zeta(i) - zeta(i+1))/(X2(i)+X2(i+1));
    Cp(i) = 1 - Q_t(i)^2/u_infty^2;        
    
end

figure
plot(x_c(2:end),Cp,'k')
set(gca,'Ydir','reverse')
xlabel('x/c')
ylabel('C_P')
%hold on