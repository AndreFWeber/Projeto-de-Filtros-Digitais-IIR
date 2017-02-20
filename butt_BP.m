clear all; close all; clc;

syms p s z

Wp1 = 0.3*pi;
Wp2 = 0.4*pi;

Ws1 = 0.2*pi;
Ws2 = 0.6*pi;

lp1 = 2*tan(Wp1/2);
lp2 = 2*tan(Wp2/2);

ls1 = 2*tan(Ws1/2);
ls2 = 2*tan(Ws2/2);

G0 = -3;
H0 = 10^(G0/20);
B = lp2 - lp1;
w0 = sqrt(lp1*lp2);

os = abs((1/B)*  ( (  (w0^2) - (ls1^2) ) /ls1)  );

Ap = 0.1;
As = 20;

n = ceil(log10((10^(As/10)-1)/(10^(Ap/10)-1))/(2*log10(os))) %pagina 200

epsilon = sqrt(10^(Ap/10) - 1);

for k = [1:n]
    pk(k) = (epsilon^-(1/n))*exp(1j*((2*k + n -1)/(2*n))*pi);
end

figure
plot(pk, 'x'); xlim([-2 2]); ylim([-2 2]);
ylabel(['Parte imaginária'],'interpreter','latex')
xlabel(['Parte real'],'interpreter','latex')
legend('Polos e zeros')
grid  on

D = real(poly(pk));
N = abs(prod(pk))*H0;

w = linspace(0,10,1000);

%Hp
Hp(p) = N./poly2sym(D, p);

%Hs
Hs(s) = collect(Hp((1/B)*((s^2 + w0^2)/s)));

%Hz
wz = linspace(0,1,1000)*pi;
Hz(z) = Hs(2*(z-1)/(z+1));

% |H(jw)|^2
Hjw = Hs(1j*w).*Hs(-1j*w);

%Plotando Hp
figure
plot(w, 20*log10(abs(Hp(1j*w)))); hold on;
plot([1 os], [-Ap -As]+G0, '+r');hold off;
legend('protótipo H(p)')
ylabel(['Magnitude em dB'],'interpreter','latex')
xlabel(['Frequencia'],'interpreter','latex')
grid on;

%Plotando Hs
figure
plot(w, 20*log10(abs(Hs(1j*w)))); hold on;
plot([ls1 lp1 ls2 lp2], [-As -Ap -As -Ap]+G0, '+r');hold off;
ylim([-70 0]);
legend('Filtro analógico H(s)')
ylabel(['Magnitude em dB'],'interpreter','latex')
xlabel(['Frequencia'],'interpreter','latex')
grid on;

%Plotando Hz
figure
plot(wz/pi, 20*log10(abs(Hz(exp(1j*wz))))); hold on;
plot([Ws1 Wp1 Ws2 Wp2]/pi, [-As -Ap -As -Ap]+G0, '+r');hold off;
ylim([-70 0]);
legend('Filtro digital H(z)')
ylabel(['Magnitude em dB'],'interpreter','latex')
xlabel(['Frequencia'],'interpreter','latex')
grid on;
