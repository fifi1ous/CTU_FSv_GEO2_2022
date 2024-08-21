
function [Mmean, Mop, V, V2, S, Smean] = zpracuj_zapisnik(M, J,  filename)
% Tato funkce zadjustuje a statisticky zpracuje zápisník vodorovných směrů.
% Vstupní proměnné jsou M - matice obsahujicí data měření a to tak, jak se
% zapisují, dále pak string vektor názvů stanovisek a název textového souboru.
% Pro 8 cílů, je tedy 16 řádků, jeden cíl se zapíše na dva řádky
% a do počtu sloupců, jako je skupin. Výstupní proměnné jsou pak hodnoty
% zprůměrovaných úhlů, hodnoty opravených úhlů, chyby, čterce chyb, S a
% S konečné a to v tomto pořadí! Kromě jiného se taky vytvoří 


[radky, sloupce] = size(M);
Mmean = ones(radky/2, sloupce);
a = 0;
for i = 1:2:radky
    a = a+1;
    for j = 1:sloupce
        if M(i, j) < M(i+1, j)
            Mmean(a, j) = (M(i, j)+M(i+1, j)-200)/2;
        else 
            Mmean(a, j) = (M(i, j)+M(i+1, j)+200)/2;
        end
    end
end

Opravy = Mmean(1, :);
Mop = mod(Mmean-Opravy, 400);

for i = 1:sloupce
    if Mop(radky/2, i) < 10
        Mop(radky/2, i) = Mop(radky/2, i)+400;
    else 
        Mop(radky/2, i) = Mop(radky/2, i);
    end
end

Mfin = mean(Mop, 2);
V = (Mfin-Mop)*10^(3);
V2 = V.^2;

S = ones(radky/2-1,1);
for i = 2:radky/2
    S(i-1) = sqrt(sum(V2(i,:))/(sloupce*(sloupce-1)));
end

S2 = S.^2;
Smean = sqrt(sum(S2)/(radky/2-1));
Mfin = compose('%3.4f', Mfin);
Mop = compose('%3.4f', Mop);
Mmean = compose('%3.4f', Mmean);
M = compose('%3.4f', M);
V = compose('%3.1f', V);
V2 = compose('%3.2f', V2);
S = compose('%3.2f', S);
Smean = compose('%3.2f', Smean);

C = cell(radky, 2*sloupce+2);
C2 = cell(radky/2-1, 2*sloupce+2);
Cvec = 1:2:radky;
for i = 1:length(J)
    C{Cvec(i), 1} = J(i);
    C{Cvec(i)+1, 1} = '';
    C{Cvec(i), 2*sloupce+2} = '';
    C{Cvec(i)+1, 2*sloupce+2} = Mfin(i);
end

for i = 1:sloupce
    for j = 1:radky
        C{j, 2*i} = M(j, i);
    end
end

V = V(2:length(V), :);
V2 = V2(2:length(V2), :);

for i = 1:sloupce
    for j = 1:radky/2
        C{Cvec(j), 2*i+1} = Mmean(j, i);
        C{Cvec(j)+1, 2*i+1} = Mop(j, i);
    end
end

for i = 1:sloupce
    for j = 1:radky/2-1
        C2{j, 2*i} = V(j, i);
        C2{j, 2*i+1} = V2(j, i);
    end
end

for j = 1:radky/2-1
    C2{j, 2*sloupce+2} = S(j);
    C2{j, 1} = J(j);
end


C = string(C);
C = C';
C2 = string(C2);
C2 = C2';
Smean = string(Smean);
fileID = fopen(filename, 'w');
mereni = fprintf(fileID, '%22s %12s %12s %12s %12s %12s %12s %12s\n', C);
Mezery = fprintf(fileID, '\n\n\n');
statistika = fprintf(fileID, '%22s %12s %12s %12s %12s %12s %12s %12s\n', C2);
finalni_S = fprintf(fileID, '%12s', Smean);
fclose(fileID);
end