function [B] = transhod_mat(P, K, body)
%P = pocatecni bod primky (nejdrive souradnice v sjtsk, potom v mistni) →
%matice 2x2
%K = koncovy bod primky (to same)
% body = matice bodu se souradnicemi v mistnim systemu
%B = matice souradnic novych bodu v sjtsk
%% delka
[r,s] = size(P);
n = r;
[r,s] = size(K);
m = r;
if n~=m
    error('matice souradnic P a K musi mit stejny pocet prvku')
end
for i = 1:n
    dY(i) = P(i,1) - K(i,1);
    dX(i) = P(i,2) - K(i,2);
    d(i) = sqrt(dY(i)^2 + dX(i)^2);
end
d = d';
%% meritko

m = d(1)/d(2);
m = round(m*10^9)/10^9

%% uhel stoceni
for i = 1:n
dY(i) = K(i,1) - P(i,1);
dX(i) = K(i,2) - P(i,2);
end
dY = dY';
dX = dX';
beta = (atan2(dY(1),dX(1)))*200/pi; %smernik pro sjtsk
alfa = (atan2(dY(2),dX(2)))*200/pi; %smernik pro mistni soustavu
if beta<0
   beta = beta + 400
elseif alfa <0
    alfa = alfa + 400
end
omega = beta - alfa                


%% posun
tY = P(1,1) - m*(P(2,2)*sin(omega) + P(2,1)*cos(omega))
tX = P(1,2) - m*(P(2,2)*cos(omega) - P(2,1)*sin(omega))

%% vypocet novych bodu
[r,s] = size(body);

for i = 1:r
    X(i) = tX + m*(body(i,2)*cos(omega) - body(i,1)*sin(omega));
    Y(i) = tY + m*(body(i,2)*sin(omega) + body(i,1)*cos(omega));
end
Y = round(Y.*10^3)./10^3;
X = round(X.*10^3)./10^3;
B = [Y', X'];

end