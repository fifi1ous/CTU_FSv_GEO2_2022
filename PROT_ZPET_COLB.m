function[S]=PROT_ZPET_COLB(p1,fi,p2,p3)
%výpočet metody protínání zpět pomocí collinsonova bodu
%vstup: p1-souřadnice bodu č1-povinná
%       fí-matice úhlů-povinná
%       p2-souřadnice bodu č2 - bod na stejné přímce jako colinsův bod
%       p3-souřadnice bodu č3
%
%       vstupní body mohou být zadány pomocí matice bodů p1[Y;X], nebo jako
%       jednotlivé body p1[Y;X], p2[Y;X] a p3[Y;X]
%výstup:    S-výsledné souřadnice bodu
if length(fi)~=3
    Q1=[5,5;7,0.5];
    Q2=[9,5;7,0.5];
    Q3=[7,3.5;7,0.5];
    figure
    subplot(1,2,1)
    plot(Q1(:,1),Q1(:,2),'Marker','*','MarkerSize',10)
    hold on
    plot(Q2(:,1),Q2(:,2),'Marker','*','MarkerSize',10)
    plot(Q3(:,1),Q3(:,2),'Marker','*','MarkerSize',10)
    xlim([4,10]);ylim([0,6])
    hold off
    Q1=[5,5;7,3];
    Q2=[9,5;7,3];
    Q3=[7,3;7,0.5];
    subplot(1,2,2)
    plot(Q1(:,1),Q1(:,2),'Marker','*','MarkerSize',10)
    hold on
    plot(Q2(:,1),Q2(:,2),'Marker','*','MarkerSize',10)
    plot(Q3(:,1),Q3(:,2),'Marker','*','MarkerSize',10)
    xlim([4,10]);ylim([0,6])
    error('špatně zadaná matice úhlů: špatný počet parametrů')
end

if nargin<3
    if nargin==2
        p3=p1(3,:);
        p2=p1(2,:);
        p1=p1(1,:);
    end
end

alfa=fi(2)-fi(1);
if alfa<0
    alfa=alfa+400;
elseif alfa>400
    alfa=alfa-400;
else
    alfa=alfa;
end
beta=fi(3)-fi(2);
if beta<0
    beta=beta+400;
elseif beta>400
    beta=beta-400;
else
    beta=beta;
end

if alfa>100
    alfa=alfa-200;
end
if beta>100
    beta=beta-200;
end

if alfa<0
    alfa=alfa+400;
elseif alfa>400
    alfa=alfa-400;
else
    alfa=alfa;
end
if beta<0
    beta=beta+400;
elseif beta>400
    beta=beta-400;
else
    beta=beta;
end

a=alfa/200*pi; b=beta/200*pi;

d13=sqrt((p1(1,1)-p3(1,1))^2+(p1(1,2)-p3(1,2))^2);
sm13p=(atan2(p3(1)-p1(1),p3(2)-p1(2)));
if sm13p<0
    sm13p=sm13p+2*pi;
end
sm13=sm13p/pi*200;
sm31p=(atan2(p1(1)-p3(1),p1(2)-p3(2)));
if sm31p<0
    sm31p=sm31p+2*pi;
end
sm31=sm31p/pi*200;

g=pi-a-b; 
if g<0
    g=g+2*pi;
elseif g>pi
    g=g-2*pi;
else
    g=g;
end
gama=g/pi*200;
d1c=d13/sin(g)*sin(a);
d2c=d13/sin(g)*sin(b);

C1=[p1(1)+sin(sm13p-b)*d1c,p1(2)+cos(sm13p-b)*d1c];
C3=[p3(1)+sin(sm31p+a)*d2c,p3(2)+cos(sm31p+a)*d2c];
Croz=abs(C1-C3);
C=[(C1(1)+C3(1))/2,(C1(2)+C3(2))/2];

smc1=atan2(p1(1)-C(1),p1(2)-C(2));
if smc1<0
    smc1=smc1+2*pi;
end
smc1u=smc1/pi*200;
smc2=atan2(p2(1)-C(1),p2(2)-C(2));
if smc2<0
    smc2=smc2+2*pi;
end
smc2u=smc2/pi*200;
smc3=atan2(p3(1)-C(1),p3(2)-C(2));
if smc3<0
    smc3=smc3+2*pi;
end
smc3u=smc3/pi*200;

omega=smc2-smc3;
if omega<0
    omega=omega+2*pi;
elseif omega>(2*pi)
    omega=omega-2*pi;
else
    omega=omega;
end
ro=smc1-smc2;
if ro<0
    ro=ro+2*pi;
elseif ro>(2*pi)
    ro=ro-2*pi;
else
    ro=ro;
end

omp=omega;   
rop=ro;
omega=omega/pi*200;
ro=ro/pi*200;
Tau=200-omega-ro;
w=pi-rop-omp;
d14=d13/sin(w)*sin(rop);
d34=d13/sin(w)*sin(omp);

sm14=sm13+omega;
if sm14<0
    sm14=sm14+400;
elseif sm14>400
    sm14=sm14-400;
else
    sm14=sm14;
end
sm34=sm31-ro;
if sm34<0
    sm34=sm34+400;
elseif sm34>400
    sm34=sm34-400;
else
    sm34=sm34;
end

sm14p=sm14/200*pi; sm34p=sm34/200*pi;

S1=[p1(1)+sin(sm14p)*d14,p1(2)+cos(sm14p)*d14];
S3=[p3(1)+sin(sm34p)*d34,p3(2)+cos(sm34p)*d34];
Sroz=abs(S1-S3);
S=[(S1(1)+S3(1))/2,(S1(2)+S3(2))/2];
%% graf
U1=[p1;S]*(-1);
U2=[p2;S]*(-1);
U3=[p3;S]*(-1);
U4=[C;S]*(-1);
figure
plot(U1(:,1),U1(:,2),'Marker','*','MarkerSize',10,'Color','b')
text((p1(:,1)*(-1)),(p1(:,2)*(-1)),'1')
text((p2(:,1)*(-1)),(p2(:,2)*(-1)),'2')
text((p3(:,1)*(-1)),(p3(:,2)*(-1)),'3')
text((C(:,1)*(-1)),(C(:,2)*(-1)),'C')
text((S(:,1)*(-1)),(S(:,2)*(-1)),'S')
hold on
plot(U2(:,1),U2(:,2),'--','Marker','*','MarkerSize',10,'Color','b')
plot(U3(:,1),U3(:,2),'Marker','*','MarkerSize',10,'Color','b')
plot(U4(:,1),U4(:,2),'Marker','*','MarkerSize',10,'Color','R')
set(gca,'xAxislocation','top');set(gca,'yAxislocation','right');
grid on
hold off
%%
fid=fopen('postup.txt','w');
fprintf(fid,'1 [%13.3f %14.3f]  %10.4f\n',p1,fi(1));
fprintf(fid,'2 [%13.3f %14.3f]  %10.4f\n',p2,fi(2));
fprintf(fid,'3 [%13.3f %14.3f]  %10.4f\n',p3,fi(3));
fprintf(fid,'Alfa v gonech: %10.4f       Beta v gonech: %10.4f    Gama v gonech: %10.4f\n',alfa,beta,gama);
fprintf(fid,'Délka mezi body 1 a 3: %21.6f m\n',d13 );
fprintf(fid,'Délka mezi body 1 a C: %21.6f m     Délka mezi bod 3 a C: %21.6f m\n',d1c,d2c);
fprintf(fid,'Směrník z bodu 1 na bod 3 v gonech: %12.6f  Směrník z bodu 3 na bod 1 v gonech: %12.6f\n',sm13,sm31);
fprintf(fid,'Souřadnice bodu C z bodu 1:[%13.5f %14.5f]\n',C1);
fprintf(fid,'Souřadnice bodu C z bodu 3:[%13.5f %14.5f]\n',C3);
fprintf(fid,'Rozdíl souřadnic bodů C   :[%13.5f %14.5f]\n',Croz);
fprintf(fid,'Souřadnice bodu C         :[%13.5f %14.5f]\n',C);
fprintf(fid,'Směrník z bodu C na 1 v gonech:    %12.6f\n',smc1u);
fprintf(fid,'Směrník z bodu C na 2 v gonech:    %12.6f\n',smc2u);
fprintf(fid,'Směrník z bodu C na 3 v gonech:    %12.6f\n',smc3u);
fprintf(fid,'Fí v gonech: %11.6f   Ró v gonech: %11.6f  Tau v gonech:%11.6f\n',omega,ro,Tau);
fprintf(fid,'Délka mezi body 1 a S: %21.6f m     Délka mezi bod 1 a S: %21.6f m\n',d14,d34);
fprintf(fid,'Souřadnice bodu S z bodu 1:[%13.5f %14.5f]\n',S1);
fprintf(fid,'Souřadnice bodu S z bodu 3:[%13.5f %14.5f]\n',S3);
fprintf(fid,'Rozdíl souřadnic bodů S   :[%13.5f %14.5f]\n',Sroz);
fprintf(fid,'Souřadnice bodu S         :[%13.5f %14.5f]\n',S);
fclose(fid)
end