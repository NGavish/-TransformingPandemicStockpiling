function main
close all;clc;clear all;

gammaI=1/2;gammaH=1/4;pHospitalization=0.12;betaParm=2*gammaI;duration=90;clinicalAttackRate=0.25;N=320e6;
computeNPI(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N,0.72)
set(gcf,'position',[1 1 917 424]);
printGraph('./graphs/NPI_impact')
end

function [Hdynamic,Hflat]=computePeak(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N)
betaParm=2*gammaI;
err=findR0(betaParm,gammaI,gammaH,pHospitalization);

beta=fzero(@(betaParm)findR0(betaParm,gammaI,gammaH,pHospitalization),betaParm);

I0=1e-3;y0=[1-I0;I0;0;0;0];
[t,y]=ode45(@(t,y)SIHR(t,y,beta,gammaI,gammaH,pHospitalization),[0 duration],y0,odeset('RelTol',1e-10));
Z=y(end,5);
[t,y]=ode45(@(t,y)SIHR(t,y,beta,gammaI,gammaH,pHospitalization*clinicalAttackRate/Z),[0 duration],y0,odeset('RelTol',1e-10));

R0=beta/gammaI;

H=y(:,3);

format bank
Hdynamic=max(H)*N;
Hflat=clinicalAttackRate*N*pHospitalization/gammaH/duration;
end

function err=findR0(betaParm,gammaI,gammaH,pHospitalization)

I0=1e-3;
y0=[1-I0;I0;0;0;0];
[t,y]=ode45(@(t,y)SIHR(t,y,betaParm,gammaI,gammaH,pHospitalization),[0 3000],y0,odeset('RelTol',1e-10,'Events', @(t,y) peakEventsFcn(t,y,betaParm,gammaI,gammaH,pHospitalization)));

err=t(end)-21;
end

function [position,isterminal,direction] = peakEventsFcn(t,y,beta,gammaI,gammaH,pHospitalization)

  dydt = SIHR(t,y,beta,gammaI,gammaH,pHospitalization);

  position = dydt(2); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

function dy=SIHR(t,y,beta,gammaI,gammaH,pHospitalization)

S=y(1);I=y(2);H=y(3);R=y(4);ZI=y(5);

dy=[-beta*I*S;
    beta*I*S-gammaI*I;
    pHospitalization*gammaI*I-gammaH*H;
    (1-pHospitalization)*gammaI*I+gammaH*H;
    gammaI*I
    ];
end


function computeNPI(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N,R0ratio)
betaParm=2*gammaI;
err=findR0(betaParm,gammaI,gammaH,pHospitalization);

beta=fzero(@(betaParm)findR0(betaParm,gammaI,gammaH,pHospitalization),betaParm);

I0=1e-3;y0=[1-I0;I0;0;0;0];
[t,y]=ode45(@(t,y)SIHR(t,y,beta,gammaI,gammaH,pHospitalization),[0 duration],y0,odeset('RelTol',1e-10));
Z=y(end,5);
[t,y]=ode45(@(t,y)SIHR(t,y,beta,gammaI,gammaH,pHospitalization*clinicalAttackRate/Z),[0 duration],y0,odeset('RelTol',1e-10));
t=t/7;Z=y(end,5)
R0=beta/gammaI;

R0*[1 R0ratio]
H=y(:,3);

x_limits=[0 245 245 0 0];
    p6 = fill(x_limits, [5e4 5e4 1e5 1e5 5e4], [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');hold on;
    p7 = fill(x_limits, [1e5 1e5 1.5e5 1.5e5 1e5], [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    p8 = yline(1.5e5, 'r--', 'LineWidth', 2);

    
% p6=fill([0 245 245 0 0],5e4+5e4*[0 0 1 1 0],'y','facealpha',0.25);hold on;
% p7=fill([0 245 245 0 0],1e5+5e4*[0 0 1 1 0],'r','facealpha',0.25);hold on;
% p8=fill([0 245 245 0 0],1e5+5e4+5e6*[0 0 1 1 0],'r','facealpha',0.5);hold on;

p3d = fill([t(:,1);flipud(t(:,1))], N*[H*0.7; flipud(H)*0.4],'k','FaceAlpha',0.4);

% plot(t,N*H,'b','LineWidth',1)


hold on;

[t,y]=ode45(@(t,y)SIHR(t,y,R0ratio*beta,gammaI,gammaH,pHospitalization*clinicalAttackRate/Z),[0 duration],y0,odeset('RelTol',1e-10));
H=y(:,3);t=t/7;
p3dR0 = fill([t(:,1);flipud(t(:,1))], N*[H*0.7; flipud(H)*0.4],'b','FaceAlpha',0.4);

axis([0 12 0 1.5e6])

ytickformat('%.0f'); % Reset first
yticks_vals = get(gca, 'YTick');
labels = arrayfun(@(v) sprintf('%.1fM', v/1000000), yticks_vals, 'UniformOutput', false);
% Fix the lower ones to k if you prefer, but uniform M is cleaner for high scale
set(gca, 'YTickLabel', labels);
set(gca, 'FontSize', 11);
xlabel('Weeks into wave', 'FontSize', 12);
ylabel('Daily Hospitalizations', 'FontSize', 12);

   legend([p3d,p3dR0,p6, p7, p8], ...
        {'Hospitalization', 'With NPIs', 'Severe Strain Zone', 'Crisis Zone', 'Severely Exceeding Capacity'}, ...
        'Location', 'northeast');

   title('Impact of NPIs on Hospital Demand')

   y(end,5)*clinicalAttackRate/Z
end