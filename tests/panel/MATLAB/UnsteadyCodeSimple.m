%% Code for predicting performance of a thin aerofoil in unsteady flow

%This example is a simple 2D aerfoil in unsteady flow modelled using lumped
%vortex elements. The singularities associated with ideal irrotational
%vortices are bypassed by the introduction of the Lamb-Oseen vortex core
%model (see Wiki for more details).

clear
close all

%Number of time steps:
Nstep = 250*2;

%Declare array to store wake data:
vortic = zeros(Nstep,3);

%Density of air:
rho = 1.0;

%X-component of the velocity of the aerofoil:
ut = 50;

%chord:
c = 1.56/ut;

%Heaving amplitude:
h0 = c*0.019;

%Heaving frequency:
omega = 2*ut/c*8.5;

%Orientation of blade with respect to global coordinate system:
beta = 0*5.0*pi/180;
%If omega = 0, this is just the aerofoil

%Sine alpha:
sn = sin(beta);

%Cosine alpha:
cs = cos(beta);

%Time step:
dt = 0.0009*c/ut;

%Time vector:
t = 0:dt:(Nstep-1)*dt;

%Z-component of the velocity of the aerofoil:
wt = h0*omega*cos(omega*t);

%Location of the origin of the body coordinate system:
sx = -ut*t;
sz = -h0*sin(omega*t);

%Location of most recently shed wake element:
Q = sqrt(ut^2 + wt.^2);
dxw = 0.3*Q*dt;

%Trailing vortex locations:
vortic(:,1) = (c + dxw)*cs + sx;
vortic(:,2) = -(c + dxw)*sn + sz;

%Influence matrix terms:
a = -1/(pi*c);
b = 1./(2*pi*(c/4 + dxw));
%Refer to chapter 9 (page 225) for more details.

xx1 = (0.75*c*cs + sx(2:Nstep))';
zz1 = (0.75*c*sn + sz(2:Nstep))';

%Location of the bound vortex:
xx1_1 = 0.25*c*cs + sx;
zz1_1 = 0.25*c*sn + sz;

%RHS1_a:
RHS1_a = -(ut*sn + wt*sn);

%Vortex core (simple representation):
cutoff = 0.001;

%Some arrays:
t1 = zeros(Nstep,1);
L = t1;
D = t1;    
cl = t1;
cd = t1;
clt = t1;
sx1 = t1;

%Run through the time steps:
for it = 1:Nstep

    rhs2 = 0;
    wwake = 0;

    if it == 1
        
    else
        
        %Calculate the induced velocity arising from all the wake elements.
        
        %Index marker:
        it1 = it - 1;
        
        %3/4 chord point. According to Anderson (see Wind Energy Handbook
        %also), the angle of attack should be taken at this point:
        AoA_x_pos = 0.75*c*cs + sx(it);
        AoA_z_pos = 0.75*c*sn + sz(it);
        
        %Distance from wake vortex elements to 3/4 chord point:
        rx = AoA_x_pos - vortic(1:it1,1);
        rz = AoA_z_pos - vortic(1:it1,2);
        r2 = rx.^2 + rz.^2 + 0*cutoff^2;
        
        %Induced velocity by each vortex wake element.
        u1 = rz/(2*pi).*vortic(1:it1,3)./r2.*(1 - exp(-r2/cutoff^2));
        w1 = -rx/(2*pi).*vortic(1:it1,3)./r2.*(1 - exp(-r2/cutoff^2));        

        %Get rid of NaN values due to singularity at core:
        u1(isnan(u1)) = 0;
        w1(isnan(w1)) = 0;   
        %This needs to be done because in spite of the vortex core model,
        %the code will still get confused at the core's centre. The
        %programme is trying to evaluate 0*0/0, which returns NaN. This
        %must be removed.
        
        %Sum all velocity contributions from the wake vortex elements:
        u = sum(u1);
        w = sum(w1);
        
        %Dot induced velocity with normal vector (of aerofoil) as per
        %equation 13.117:
        wwake = u*sn + w*cs;
        
        rhs2 = -sum(vortic(1:it1,3));
        
    end
    
    %Equations 13.114 - 13.119:
    rhs1 = -ut*sn - wwake - wt(it)*cs;
    vortic(it,3) = 1/(b(it)/a - 1)*(rhs1/a-rhs2);
    gammat = rhs2 - vortic(it,3);
    %For a single lumped vortex representation of an aerofoil, equation
    %13.118 reduces to a simple expression in which the 'a' matrix is just
    %a 2-by-2 matrix of the following form:
    
    % a = a_11 a_1W
    %       1    1
    
    %The column vector of unknowns is just [gamma_1,gamma_w_i]^T. The RHS
    %is jused [RHS1,RHS2]^T. Matrix inversion of a 2-by-2 matrix does not
    %even require Matlab \ command. The inverse of the 'a' matrix is just
    
    % a^{-1} = 1/(a_11 - a_1W)*[1, -a_1W; -1, a_11]
    
    %Note - although this is a 'pure' model, it is quite possible to
    %upgrade it to a LLT model. To do so, you simply need to recognise that
    %the effective angle of attack (taking into account both kinematic and
    %wake-induced velocities) is known, and then you simply go to a lookup
    %table to get the lift and drag. Then by the Kutta-Joukowski theorem,
    %you calculate the strength of the shed vortex and feed that into the
    %wake. Using lookup tables does, however, have its own issues. For
    %example, the use of lookup tables may be a source of error in
    %curvilinear flow or when a wing/blade is of finite length.
    
    if it < 1
        
    else        

        %Velocity induced by the bound vortex on the aerofoil at the
        %positions of all the wake elements:
        Bound_x = 0.25*c*cs + sx(it);
        Bound_z = 0.25*c*sn + sz(it);        
        
        %Distance from wake vortex elements to aerofoil bound vortex:
        rx = vortic(1:it,1) - Bound_x;
        rz = vortic(1:it,2) - Bound_z;        
        r2 = rx.^2 + rz.^2 + 0*cutoff^2;
        
        %Induced velocity on each vortex wake element by bound vortex.
        u = rz/(2*pi)*gammat./r2.*(1 - exp(-r2/cutoff^2));
        w = -rx/(2*pi)*gammat./r2.*(1 - exp(-r2/cutoff^2));

        %Get rid of NaN values due to singularity at core:
        u(isnan(u)) = 0;
        w(isnan(w)) = 0;        
        
        %Determine the influence of the wake elements on each other:
        for i = 1:it

            %Distance between all wake elements and the wake element of
            %interest:
            rx = vortic(i,1) - vortic(1:it,1);
            rz = vortic(i,2) - vortic(1:it,2);
            r2 = rx.^2 + rz.^2 + 0*cutoff^2;
            
            u11 = rz/(2*pi).*vortic(1:it,3)./r2.*(1 - exp(-r2/cutoff^2));
            w11 = -rx/(2*pi).*vortic(1:it,3)./r2.*(1 - exp(-r2/cutoff^2));
            
            %Get rid of NaN values due to singularity at core:
            u11(isnan(u11)) = 0;
            w11(isnan(w11)) = 0;
            %This needs to be done because in spite of the vortex core
            %model, the code will still get confused at the core's centre.
            %The programme is trying to evaluate 0*0/0, which returns NaN. 
            %This must be removed.
            
            %Sum influence of all wake elements on the element of interest:
            u1 = sum(u11);
            w1 = sum(w11);
            
            u(i) = u(i) + u1;
            w(i) = w(i) + w1;
            
            %uw(i,1) = vortic(i,1) + u(i)*dt;
            %uw(i,2) = vortic(i,2) + w(i)*dt;
            
        end
                
        %Move the wake elements with the local velocities:
        vortic(1:it,1) = vortic(1:it,1) + u(1:it)*dt;
        vortic(1:it,2) = vortic(1:it,2) + w(1:it)*dt;
        
    end
    
    if it == 1
        
        %No circulation at t = 0:
        gamat1 = 0;
        
    end
    
    %Kinematic velocity:
    kin_vel = 0.5*rho*(ut^2 + wt(it)^2);

    %Rate of change of bound circulation with respect to time:
    dgamdt = (gammat - gamat1)/dt;
    
    %Store current value of bound circulation for use in the unsteady
    %Bernoulli equation at the next time step:
    gamat1 = gammat;
    
    AoA_x_pos = 0.75*c*cs + sx(it);
    AoA_z_pos = 0.75*c*sn + sz(it);
    
    %Calculate downwash:
    rx = AoA_x_pos - vortic(1:it,1);
    rz = AoA_z_pos - vortic(1:it,2);
    r2 = rx.^2 + rz.^2 + cutoff^2;
    
    u11 = rz/(2*pi).*vortic(1:it,3)./r2;
    w11 = -rx/(2*pi).*vortic(1:it,3)./r2;
            
    %Sum influence of all wake elements on the element of interest:
    u = sum(u11);
    w = sum(w11);
    
    %Wake-induced downwash:
    ww = u*sn + w*cs;
    
    %Lift force (equation 13.128):
    L(it) = rho*(ut*gammat + dgamdt*c);
    
    %Drag force (equation 13.130):
    D(it) = rho*(-ww*gammat + dgamdt*c*sn);
    %Page 415 - The first term is due to the wake-induced downwash, which
    %in the lumped-vortex element case is evaluated at the panel's
    %three-quarter chord point. The second term is due to the fluid
    %acceleration.
    
    %Lift coefficient:
    cl(it) = L(it)/kin_vel/c;
    
    %Drag coefficient:
    cd(it) = D(it)/kin_vel/c;
    
    %Note the drag force is a lift-induced drag. It is because of the
    %vorticity of the wake.
    
    sx1(it) = sx(it) - ut*dt;
    
end

if dt > 0.001*c/ut

    figure
    plot(vortic(:,1)/c,vortic(:,2)/c,'--')

else

    figure
    plot(vortic(:,1)/c,vortic(:,2)/c,'x')
    
end

%Add some labels:
xlabel('x/c')
ylabel('z/c')

if omega == 0
    
    %set some limits:
    ylim([-0.75 0.75])
    
end

figure
plot(ut*t/c,cl)

%Add some labels:
xlabel('u_\infty t/c')
ylabel('C_l')

figure
plot(ut*t/c,cd)

%Add some labels:
xlabel('u_\infty t/c')
ylabel('C_d')