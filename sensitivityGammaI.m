function main
close all;clc;clear all;

gammaI=1/2;gammaH=1/4;pHospitalization=0.12;betaParm=2*gammaI;duration=90;clinicalAttackRate=0.25;N=320e6;midPeriod=21;
computeNPI(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N,0.72,midPeriod)

%% CDC
'CDC'
gammaI=1/2;gammaH=1/4;pHospitalization=0.12;betaParm=2*gammaI;duration=42;clinicalAttackRate=0.25;N=320e6;midPeriod=21;
[Hdynamic,Hflat]=computePeak(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N,midPeriod)

gammaIVec=1./linspace(1.5,10,50);
for ix=1:numel(gammaIVec)
    [Hdynamic(ix),Hflat(ix)]=computePeak(pHospitalization,duration,clinicalAttackRate,gammaIVec(ix),gammaH,N,midPeriod);
end

subplot(2,2,1)
plot(1./gammaIVec,Hdynamic,1./gammaIVec,Hflat,'--',LineWidth=1);xlim([1 10.5]);xlabel('Average infectious period (days)')
ylabel('peak capacity demand')
legend('H^{dynamic}','H^{flat}','Location','best')
printGraph('./graphs/sensitivityGammaI')
%% UK
'UK'
pHospitalization=0.04;duration=15*7;clinicalAttackRate=0.5;N=70e6;midPeriod=52;gammaI=1/5;
[Hdynamic,Hflat]=computePeak(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N,midPeriod)
% 160,000 hospital beds in normal operation, about 13,000 vacant

%% Japan
'Japan'
pHospitalization=0.078;duration=8*7;clinicalAttackRate=0.25;N=124e6;midPeriod=28;
[Hdynamic,Hflat]=computePeak(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N,midPeriod)
% 920,000 beds, about 184,000 vacant
end

function [Hdynamic,Hflat]=computePeak(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N,midPeriod)
betaParm=2*gammaI;
err=findR0(betaParm,gammaI,gammaH,pHospitalization,midPeriod);

beta=fzero(@(betaParm)findR0(betaParm,gammaI,gammaH,pHospitalization,midPeriod),betaParm);

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

function err=findR0(betaParm,gammaI,gammaH,pHospitalization,midPeriod)

I0=1e-3;
y0=[1-I0;I0;0;0;0];
[t,y]=ode45(@(t,y)SIHR(t,y,betaParm,gammaI,gammaH,pHospitalization),[0 3000],y0,odeset('RelTol',1e-10,'Events', @(t,y) peakEventsFcn(t,y,betaParm,gammaI,gammaH,pHospitalization)));

err=t(end)-midPeriod;
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


function computeNPI(pHospitalization,duration,clinicalAttackRate,gammaI,gammaH,N,R0ratio,midPeriod)
betaParm=2*gammaI;
err=findR0(betaParm,gammaI,gammaH,pHospitalization,midPeriod);

beta=fzero(@(betaParm)findR0(betaParm,gammaI,gammaH,pHospitalization,midPeriod),betaParm);

I0=1e-3;y0=[1-I0;I0;0;0;0];
[t,y]=ode45(@(t,y)SIHR(t,y,beta,gammaI,gammaH,pHospitalization),[0 duration],y0,odeset('RelTol',1e-10));
Z=y(end,5);
[t,y]=ode45(@(t,y)SIHR(t,y,beta,gammaI,gammaH,pHospitalization*clinicalAttackRate/Z),[0 duration],y0,odeset('RelTol',1e-10));
t=t/7;Z=y(end,5)
R0=beta/gammaI;

R0*[1 R0ratio]
H=y(:,3);

p6=fill([0 245 245 0 0],5e4+5e4*[0 0 1 1 0],'y','facealpha',0.25);hold on;
p7=fill([0 245 245 0 0],1e5+5e4*[0 0 1 1 0],'r','facealpha',0.25);hold on;
p8=fill([0 245 245 0 0],1e5+5e4+5e6*[0 0 1 1 0],'r','facealpha',0.5);hold on;

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
        {'Hospitalization', 'With NPIs', 'Severe Strain Zone', 'Crisis Zone', 'Severe Overcapacity'}, ...
        'Location', 'northeast');

   title('Impact of NPIs on Hospital Demand')

   y(end,5)*clinicalAttackRate/Z
end