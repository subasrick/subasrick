clc;
gamma_max = 0.93;
c = 1;
f=1e-5;
cmax = 20;
nadd = 20;
Ts=5000;
Kp = 40;
Ki = 1;
Kd = 1;
umim=0;
umax=10;
ymin=0;
ymax=10;
tav=0.01;
ddead=0.001;
ap=0.02;
bp=0.02;
sigmalL=0.2;
utheta=0.2;
ltheta=0.1;
Tstop=5000;
deltay=ymax-ymin;
deltawps=deltay/2;

load ydata.mat
load ydataref.mat
track_error=ydataref.signals.values-ydata.signals.values;
xk=[track_error/deltawps track_error-ymin/deltay];
mean_val=xk(1);
Mi=1;
gsign=1;
mean_square=sqrt(xk).^2;
medata=(((Mi-1)/Mi)*mean_val)+(1/Mi)*xk;
if(c==0)
    c=c+1;
else
    sigamedata=(((Mi-1)/Mi)*mean_val)+(1/Mi)*mean_square;
    lamdata=1/1+((sqrt(xk-medata).^2)+sigamedata-(sqrt(medata).^2));
    lamdanew=lamdata/sum(sum(lamdata));
    
    if(max(lamdata)<max(lamdanew))
        c=c+1;
        fmax=Ki/2;
    else
        fmax=Ki/2;
    end
end
alpha=((umax-umim)/20).*0.1;
deltap=Kp*alpha*gsign.*lamdanew(:,1).*(ydataref.signals.values.*track_error)./ydataref.signals.values;
deltap=sort(deltap/f);
idatap=Ki*alpha*gsign.*lamdanew(:,1).*(ydataref.signals.values.*track_error)./ydataref.signals.values;
idatap=sort(idatap/f);
kdatap=Kd*alpha*gsign.*lamdanew(:,1).*(ydataref.signals.values.*track_error)./ydataref.signals.values;
kdatap=sort(kdatap/f);
open_system('Original')
sim('Original',Tstop)
figure
plot(tout_i(1:500),ydataref.signals.values(1:500),'-r','LineWidth',2)
hold on
plot(tout_i(1:500),cdata.signals.values(1:500),'-b','LineWidth',2)
xlabel('Time (s)')
ylabel('')
legend('Reference','Controlled signal')
title('Reference & Controlled signal')
figure
plot(tout_i(1:500),ydata.signals.values(1:500),'-b','LineWidth',2)
xlabel('Time (s)')
ylabel('')
title('Control signal')
figure
plot(tout_i(1:5000),track_error(1:5000),'-b','LineWidth',2)
xlabel('Time (s)')
ylabel('')
title('Tracking Error')
figure
plot(tout_i(1:5000),deltap(1:5000),'-b','LineWidth',2)
xlabel('Time (s)')
ylabel('P')
title('P')
figure
plot(tout_i(1:5000),idatap(1:5000),'-b','LineWidth',2)
xlabel('Time (s)')
ylabel('I')
title('I')
figure
plot(tout_i(1:5000),kdatap(1:5000),'-b','LineWidth',2)
xlabel('Time (s)')
ylabel('D')
title('D')
