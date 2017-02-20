
clear all; close all; clc;

syms p s z

ws = 0.8*pi;
wp = 0.4*pi;

ls = 2*tan(ws/2);
lp = 2*tan(wp/2);

os = ls/lp;

As = 35;
Ap = 1;

epsilon = sqrt(10^(Ap/10)-1);

n = ceil(acosh(sqrt((10^(As/10) - 1)/(10^(Ap/10)-1)))/acosh(os));

fi2 = (1/n)*asinh(1/epsilon);
for k =[1:n]
    fi1 = ((2*k-1)*pi)/(2*n);
    pk(k) = -sinh(fi2)*sin(fi1) + 1j*cosh(fi2)*cos(fi1);
end

figure
plot(pk, 'x'); xlim([-2 2]); ylim([-2 2]); grid on;
ylabel(['Parte imaginária'],'interpreter','latex')
xlabel(['Parte real'],'interpreter','latex')
legend('Polos e zeros')

D = real(poly(pk));
N = abs(prod(pk));

w = linspace(0,10,1000)/lp;

%Hp
Hp(p) = N./poly2sym(D, p);

%Hs
Hs(s) = Hp(s/lp);

%Hz
wz = linspace(0,1,1000)*pi;
Hz(z) = Hs(2*(z-1)/(z+1));

% |H(jw)|^2
Hjw = Hs(1j*w).*Hs(-1j*w);

%|H(ω)|2 = 1/[1 + (ω / ωc ) 2N)] 

s= w*1j;
p = s/wp;
z= wz*1j

%Plotando Hp
figure
plot(w,20*log10(abs(Hp(s))));hold on;
plot([0 1 os], [0 -Ap -As], '+r');hold off;
legend('protótipo H(p)')
ylabel(['Magnitude em dB'],'interpreter','latex')
xlabel(['Frequencia'],'interpreter','latex')
grid on;

%Plotando Hs
figure
plot(w, 20*log10(abs(Hs(s)))); hold on;
plot([ls lp], [-As -Ap], '+r');hold off;
legend('Filtro analógico H(s)')
ylabel(['Magnitude em dB'],'interpreter','latex')
xlabel(['Frequencia'],'interpreter','latex')
grid on;

%Plotando Hz
figure
plot(wz/pi, 20*log10(abs(Hz(exp(1j*wz))))); hold on;
plot([ws wp]/pi, [-As -Ap], '+r');hold off;
ylim([-50 2]);
legend('Filtro digital H(z)')
ylabel(['Magnitude em dB'],'interpreter','latex')
xlabel(['Frequencia'],'interpreter','latex')
grid on;


