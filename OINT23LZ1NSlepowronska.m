syms x v epsilon eta2 eta4 eta_plus eta_sin eta_exp eta_cos eta_ln eta_mult eta_div eta_sub


n = 1000;
m = 10000;
%przedzial x
x1 = linspace(0,1,n);

%funkcja y zadana przez prowadzacego
y = sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3));

y_diff = diff(y,x);
T = x*y_diff;
T = T/y;
%T obliczony metoda rozniczkowania analitycznego
T = simplify(T);

%T obliczony metoda rachunku epsilonow
Te = @(x) x.*(2.*x.*exp(x.^4).*log(x+3).^2.*(cos(x.^2)+sin(x.^2).*2.*x.^2)+sin(x+3).*log(x+3)+cos(x+3)./(x+3))./(sin(x.^2).*exp(x.^4).*log(x+3).^2-cos(x+3).*log(x+3));

%% 

%zadanie1

%wspolczynnik K1 obliczony metoda rachunku epsilonow
K1e = @(x) (x.^2.*cos(x.^2).*exp(x.^4))./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y1 = sin(v).*exp(x.^4)-(cos(x+3)./log(x+3));
K1 = v.*diff(y1,v)./y1;
%wspolczynnik K1 obliczony metoda rozniczkowania analitycznego
K1 = simplify(subs(K1,v,x.^2))

%wspolczynnik K2 obliczony metoda rachunku epsilonow
K2e = @(x) (sin(x.^2).*x.^4.*exp(x.^4))./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y2 = sin(x^2).*exp(v)-(cos(x+3)./log(x+3));
K2 = v.*diff(y2,v)./y2;
%wspolczynnik K2 obliczony metoda rozniczkowania analitycznego
K2 = simplify(subs(K2,v,x.^4))

%wspolczynnik K3 obliczony metoda rachunku epsilonow
K3e = @(x) ((1./log(x+3)+tan(x+3).*(x+3)).*(cos(x+3)./log(x+3)))./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y3 = sin(x^2).*exp(x.^4)-(cos(v)./log(v));
K3 = v.*diff(y3,v)./y3;
%wspolczynnik K3 obliczony metoda rozniczkowania analitycznego
K3 = simplify(subs(K3,v,x+3))

%wspolczynnik K4 obliczony metoda rachunku epsilonow
K4e = @(x) sin(x.^2).*exp(x.^4)./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y4 = v.*exp(x.^4)-(cos(x+3)./log(x+3));
K4 = v.*diff(y4,v)./y4;
%wspolczynnik K4 obliczony metoda rozniczkowania analitycznego
K4 = simplify(subs(K4,v,sin(x.^2)))

%wspolczynnik K5 obliczony metoda rachunku epsilonow
K5e = @(x) sin(x.^2).*exp(x.^4)./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y5 = sin(x^2).*v-(cos(x+3)./log(x+3));
K5 = v.*diff(y5,v)./y5;
%wspolczynnik K5 obliczony metoda rozniczkowania analitycznego
K5 = simplify(subs(K5,v,exp(x.^4)))

%wspolczynnik K6 obliczony metoda rachunku epsilonow
K6e = @(x) (-cos(x+3)./log(x+3))./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y6 = sin(x^2).*exp(x^4)-(v./log(x+3));
K6 = v.*diff(y6,v)./y6;
%wspolczynnik K6 obliczony metoda rozniczkowania analitycznego
K6 = simplify(subs(K6,v,cos(x+3)))

%wspolczynnik K7 obliczony metoda rachunku epsilonow
K7e = @(x) (cos(x+3)./log(x+3))./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y7 = sin(x^2).*exp(x^4)-(cos(x+3)./v);
K7 = v.*diff(y7,v)./y7;
%wspolczynnik K7 obliczony metoda rozniczkowania analitycznego
K7 = simplify(subs(K7,v,log(x+3)))

%wspolczynnik K8 obliczony metoda rachunku epsilonow
K8e = @(x) sin(x.^2).*exp(x.^4)./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y8 = v-(cos(x+3)./log(x+3));
K8 = v.*diff(y8,v)./y8;
%wspolczynnik K8 obliczony metoda rozniczkowania analitycznego
K8 = simplify(subs(K8,v,sin(x^2).*exp(x^4)))

%wspolczynnik K9 obliczony metoda rachunku epsilonow
K9e = @(x) (-cos(x+3)./log(x+3))./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y9 = sin(x^2).*exp(x^4)-(v);
K9 = v.*diff(y9,v)./y9;
%wspolczynnik K9 obliczony metoda rozniczkowania analitycznego
K9 = simplify(subs(K9,v,cos(x+3)./log(x+3)))

%wspolczynnik K10 obliczony metoda rachunku epsilonow
K10e = @(x) (sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)))./(sin(x.^2).*exp(x.^4)-(cos(x+3)./log(x+3)));
y10 = v;
K10 = v.*diff(y10,v)./y10;
%wspolczynnik K10 obliczony metoda rozniczkowania analitycznego
K10= simplify(subs(K10,v,sin(x.^2).*exp(x.^4) - (cos(x+3)./log(x+3)))) 
%% 

%wykresy

figure(1)
hold off
fplot(K1, [0, 1], 'g');
hold on
plot(x1,K1e(x1),'blue')
legend({'K1(x) - rozniczkowanie analityczne', 'K1(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K1(x)')

figure(2)
hold off
fplot(K2, [0, 1], 'g');
hold on
plot(x1,K2e(x1),'blue')
legend({'K2(x) - rozniczkowanie analityczne', 'K2(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K2(x)')

figure(3)
hold off
fplot(K3, [0, 1], 'g');
hold on
plot(x1,K3e(x1),'blue')
legend({'K3(x) - rozniczkowanie analityczne', 'K3(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K3(x)')

figure(4)
hold off
fplot(K4, [0, 1], 'g');
hold on
plot(x1,K4e(x1),'blue')
legend({'K4(x) - rozniczkowanie analityczne', 'K4(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K4(x)')

figure(5)
hold off
fplot(K5, [0, 1], 'g');
hold on
plot(x1,K5e(x1),'blue')
legend({'K5(x) - rozniczkowanie analityczne', 'K5(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K5(x)')

figure(6)
hold off
fplot(K6, [0, 1], 'g');
hold on
plot(x1,K6e(x1),'blue')
legend({'K6(x) - rozniczkowanie analityczne', 'K6(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K6(x)')

figure(7)
hold off
fplot(K7, [0, 1], 'g');
hold on
plot(x1,K7e(x1),'blue')
legend({'K7(x) - rozniczkowanie analityczne', 'K7(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K7(x)')

figure(8)
hold off
fplot(K8, [0, 1], 'g');
hold on
plot(x1,K8e(x1),'blue')
legend({'K8(x) - rozniczkowanie analityczne', 'K8(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K8(x)')

figure(9)
hold off
fplot(K9, [0, 1], 'g');
hold on
plot(x1,K9e(x1),'blue')
legend({'K9(x) - rozniczkowanie analityczne', 'K9(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K9(x)')

figure(10)
hold off
fplot(K10, [0, 1], 'g');
hold on
plot(x1,K10e(x1),'blue')
legend({'K10(x) - rozniczkowanie analityczne', 'K10(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik K10(x)')

figure(11)
hold off
fplot(T, [0, 1], 'g');
hold on
plot(x1,Te(x1),'blue')
legend({'T(x) - rozniczkowanie analityczne', 'T(x) - rachunek epsilonow'},'Location','northwest','NumColumns',2)
title('wspolczynnik T(x)')

%% 

%zadanie 2

eps = 2*10^(-13);
abs_T = @(x1) abs(Te(x1)); 
abs_K1 = @(x1) abs(K1e(x1)); 
abs_K2 = @(x1) abs(K2e(x1)); 
abs_K3 = @(x1) abs(K3e(x1)); 
abs_K4 = @(x1) abs(K4e(x1)); 
abs_K5 = @(x1) abs(K5e(x1)); 
abs_K6 = @(x1) abs(K6e(x1)); 
abs_K7 = @(x1) abs(K7e(x1)); 
abs_K8 = @(x1) abs(K8e(x1)); 
abs_K9 = @(x1) abs(K9e(x1)); 
abs_K10 = @(x1) abs(K10e(x1)); 

sum_abs = zeros(1,n);

for i = 1:1:n
    sum_abs(i) = abs_T(x1(i))+ abs_K1(x1(i)) + abs_K2(x1(i)) + abs_K3(x1(i)) + abs_K4(x1(i)) + abs_K5(x1(i)) + abs_K6(x1(i)) + abs_K7(x1(i)) + abs_K8(x1(i)) + abs_K9(x1(i)) + abs_K10(x1(i));
end

%blad calkowity
y_sup1 = max(sum_abs)*eps

%%

%zadanie 3

y = @(x) sin(x^2)*exp(x^4)-(cos(x+3)/log(x+3));
y_zaburzony = @(x, epsilon, eta2, eta4, eta_plus, eta_sin, eta_exp, eta_cos, eta_ln, eta_mult, eta_div, eta_sub) ( ( sin(((x*(1+epsilon))^2)*(1+eta2))*(1+eta_sin) * exp(((x*(1+epsilon))^4)*(1+eta4))*(1+eta_exp) )*(1+eta_mult) - ( cos((x*(1+epsilon)+3)*(1+eta_plus))*(1+eta_cos) / log((x*(1+epsilon)+3)*(1+eta_plus))*(1+eta_ln) )*(1+eta_div) )*(1+eta_sub);

x2 = linspace(0,1,m);
tab = (de2bi(0:(2^11-1))-0.5)*2;
matr = @(p,e) eps*tab(p,e);
%wszystkie mozliwe kombinacje wartosci dla wzglednych bledow danych i
%bledow zaokraglen zapisane w macierzy
blad = zeros(2^11, m);

for j = 1:(2^11)
    for k = 1:m
        blad(j,k) = abs(( y_zaburzony(x2(k), matr(j,1), matr(j,2), matr(j,3), matr(j,4), matr(j,5), matr(j,6), matr(j,7), matr(j,8), matr(j,9), matr(j,10), matr(j,11)) - y(x2(k))) / y(x2(k)));
    end
end 

y_sup2 = max(blad(:))
