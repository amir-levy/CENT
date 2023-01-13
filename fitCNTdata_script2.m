%% Data from PRL 116, 154501 (2016), Secchi et al
close all

figure(1); subplot(1,4,1);
figure; 
text(0.8, 0.1,'$3.5$nm','interpreter','latex');

initGuess=[0.028 -10 0.1];
%initGuess=[ 0.0134  -8.7    0.2];
%initGuess=[ 0.0134  -8.7    0.1];


%3.5nm data:
R=35;
L=30000;
clear sigma_A
%initGuess=[ 0.0134  -8.7    0.1];
c=[0.0088    0.0193    0.0481    0.1355    0.3458    0.9632]*1000; 
G=[0.0104    0.0179    0.0259    0.0452    0.0863    0.2148]; 
sigma_A(1,:)=fit1DTransport(c,G,L,R,initGuess,'k'); hold on
hold on; plot(c,G,'Ok','linewidth',2,'markersize',5);

%initGuess=[ 0.0134  -8.7    1];
c=[0.0026    0.0070    0.0188    0.0493    0.1338    0.3373    0.9753]*1000; 
G=[ 0.0201    0.0203    0.0290    0.0414    0.0732    0.1230    0.3142]; 
sigma_A(2,:)=fit1DTransport(c,G,L,R,initGuess,'r'); hold on
hold on; plot(c,G,'Or','linewidth',2,'markersize',5);

initGuess=[ 0.0134  -8.7    0.1];

c=[    0.0026    0.0067    0.0188    0.0487    0.1338    0.3415    0.9395]*1000;
G=[ 0.1776    0.2260    0.2987    0.5701    1.1891    1.4381    2.1303];
sigma_A(3,:)=fit1DTransport(c,G,L,R,initGuess,'g'); hold on
hold on; plot(c,G,'Og','linewidth',2,'markersize',5);

c=[0.0026    0.0188    0.0475    0.1355    0.3458    0.9513]*1000;
G=[ 0.9112    1.7393    2.2986    2.8154    3.3621    4.2237]; 
sigma_A(4,:)=fit1DTransport(c,G,L,R,initGuess,'b');
hold on; plot(c,G,'Ob','linewidth',2,'markersize',5);


xlim([0.5 5000]);
 axis square
 ylabel('$G$[nS]','interpreter','latex','fontsize',15)
xlabel('$c$(mM)','interpreter','latex')

  muFit{1}=10^(-29)*exp(sigma_A(:,2))*L^2/e^2; 
 muFit{1}/(4.8e11)
 sigmaFit=sigma_A(:,1)*0.16*16*1000; 
%%
close all
initGuess=[0.028 -10 0.4];
%initGuess=[0.028 -10 0.1];

R=7.5;
L=100;

c=[3	1458.9
1.5	1045.5
1	766.4
0.6	551.1
0.4	434.7
0.2	333.3
0.1	213.8
0.05	95.6
0.02	55.5
0.01	48.2]

err=[46.0
33.9
48.8
22.1
10.1
64.2
39.1
31.6
26.6
24.1]; 

sigma_A(1,:)=fit1DTransport(c(2:end,1)'*1000,c(2:end,2)'/1000,L,R,initGuess,'b',7,10, 4000);
hold on; errorbar(c(:,1)*1000,c(:,2)/1000,err/1000,'Ob','linewidth',2,'markersize',5);


c(:,2)=[538.3
401.1
241.8
193.5
142.4
90.8
50.0
22.6
14.0
10.8
];
err=[153.2
68.7
72.4
20.2
13.5
6.8
7.8
7.8
9.3
7.1
];

sigma_A(2,:)=fit1DTransport(c(2:end,1)'*1000,c(2:end,2)'/1000,L,R,initGuess,'k',7,10, 3000);
hold on; errorbar(c(:,1)*1000,c(:,2)/1000,err/1000,'Ok','linewidth',2,'markersize',5);
xlim([10 1000]); 
muFit{1}=10^(-29)*(exp(sigma_A(:,2)))*L^2/e;
 muFit{1}/(1e-8)
 sigma_A(:,1)*0.16*16*1000
 axis square
 ylabel('$G$[nS]','interpreter','latex','fontsize',15)
xlabel('$c$(mM)','interpreter','latex')

%%
%10nm data
figure(1); subplot(1,4,2);
text(0.8, 0.1,'$10$nm','interpreter','latex');

R=100;
L=25000;
initGuess=[0.028 -2 1];
%initGuess=[15.4937  -5    0.3];

L=25000;
c=[    0.0088    0.0189    0.0439    0.0915    0.2064    0.4477    0.9522]*1000;
G=[    0.0195    0.0274    0.0522    0.0865    0.2190    0.6152    1.9605];
sigma_B=fit1DTransport(c,G, L, R,initGuess);
hold on; plot(c,G,'Ok','linewidth',2,'markersize',5);
xlim([0.5 5000]);
 axis square
xlabel('$c$(mM)','interpreter','latex')

 muFit{2}=10^(-29)*exp(sigma_B(:,2))*L^2/e^2; 
%%
figure(1); ax=subplot(1,4,3);

%initGuess=[0.028 1];
initGuess=[ 0.0444  -12  1];

R=140;
L=20000/2;
c=[ 0.0010    0.0096    0.1009    0.9592]*1000;
G=[  0.1909    0.3975    0.9389    4.3009];
sigma_C=fit1DTransport(c,G, L, R,initGuess,'k');
hold on; plot(c,G,'Ok','linewidth',2,'markersize',5);
text(0.8, 0.1,'$14$nm','interpreter','latex', 'Units','Normalized' );

c=[ 0.0010    0.0096    0.1009    0.9592]*1000;
G=[ 0.65    0.8745    1.5430    5.1970];
err=[ 0.4581    0.8541    0.6636    1.0651    1.3601    1.7505    5.4919    5.8959];
sigma_C(2,:)=fit1DTransport(c,G, L, R,initGuess,'r');
%hold on; errorbar(c,G,G-err(1:2:end), err(2:2:end)-G,'Or','linewidth',2,'markersize',5);
hold on; plot(c,G,'Or','linewidth',2,'markersize',5);

c=[ 0.0010    0.0096    0.1009    0.9592]*1000;
G=[ 1.2471    1.6434    2.3435    5.5792];
sigma_C(3,:)=fit1DTransport(c,G, L, R,initGuess,'b');
hold on; plot(c,G,'Ob','linewidth',2,'markersize',5);
xlim([0.5 5000]);
xlabel('$c$(mM)','interpreter','latex')
ylabel('$G$[nS]','interpreter','latex','fontsize',15)

axis square
xlabel('$c$(mM)','interpreter','latex')

 muFit{3}=10^(-29)*exp(sigma_C(:,2))*L^2/e^2; 
 mu=10^(-29)*exp(sigma_C(:,2))*L^2/e/(1e-8)
 sigmaFit=sigma_C(:,1)*0.16*16*1000
%%
figure(1); %ax=subplot(1,4,4);
R=350;
L=15000;

initGuess=[0.028 1 1];

c=[0.0010    0.0099    0.0985    0.9559]*1000;
G=[ 0.1222    0.6704    4.8025   42.1554];
sigma_D=fit1DTransport(c,G, L, R,initGuess,'k');
loglog(c,G,'Ok','linewidth',2,'markersize',5); hold on;

c=[0.0010    0.0099    0.0985    0.9559]*1000;
G=[0.5070    1.1144    4.3384   33.1156];
sigma_D(2,:)=fit1DTransport(c,G, L, R,initGuess,'r');
hold on; loglog(c,G,'Or','linewidth',2,'markersize',5);

c=[0.0010    0.0099    0.0985    0.9559]*1000;
G=[0.9447    0.8979    3.4081   24.4134];
sigma_D(3,:)=fit1DTransport(c,G, L, R,initGuess,'b');
loglog(c,G,'Ob','linewidth',2,'markersize',5);

%axis square; 
xlabel('$c$(mM)','interpreter','latex')

 muFit{4}=10^(-29)*exp(sigma_D(:,2))*L^2/e^2; 
 sigmaFit=sigma_D(:,1)*0.16*16*1000; 
 
 
  mu=10^(-29)*exp(sigma_D(:,2))*L^2/e/(1e-8)
 sigmaFit=sigma_D(:,1)*0.16*16*1000
 %% BNNT data (doi:10.1038/nature11876, Siria et al, Nature Letters)
figure(2); subplot(1,6,6);
initGuess=[0.028 1 1];

R=150;
L=8000;

c=[ 0.0001    0.0010    0.0030    0.0102    0.0309    0.1045    0.3259    1.0000]*1000;
g=[1.8443    1.8443    1.8443    1.8443    1.9267    2.4683    8.5188   15.9416];
G=g*(pi*R^2)/L*1e-10;
G=G*1e9; ;
sigmaA=fit1DTransport(c,G, L, R,initGuess,'k');
hold on; plot(c,G,'Ok','linewidth',1,'markersize',5);

 muFit{6}(1)=10^(-29)*exp(sigmaA(1,2))*L^2/e^2; 

R=220;
L=15000;
c=[    0.0005    0.0050    0.0102    0.0516    0.1062    1.0000]*1000;
g=[ 0.7693    0.7806    1.2626    2.3627    3.5019   21.0274]; 
G=g*(pi*R^2)/L*1e-10;
G=G*1e9; 
sigmaA(2,:)=fit1DTransport(c,G, L, R,initGuess,'b');
hold on; plot(c,G,'Ob','linewidth',1,'markersize',5);


xlim([0.5 5000]);
xlabel('$c$(mM)','interpreter','latex')
ylabel('$G$[nS]','interpreter','latex','fontsize',15)
axis square; 
 muFit{6}(2)=10^(-29)*exp(sigmaA(2,2))*L^2/e^2; 
 sigmaFit=sigmaA(:,1)*0.16*16*1000; 

%% fig 2: Data from Tunuguntla et al (Science 2016)
close all
figure(2); 
%figure; 
text(0.8, 0.1,'$14$\AA','interpreter','latex');

initGuess=[0.028 1 0.0949];
%initGuess=[0   1    0.03];
%initGuess=[ 0.0292    -5.5 1];

R=3.4;
L=106;

c=[0.0053    0.0101    0.0205    0.0240    0.0498    0.0995    0.1989    0.5957 1.0017    2.9315]*1000;
G=[4.9434    6.9112   10.6199   16.1791   22.4257   32.4485   38.5319   49.4345  66.2057   69.7082];
sigmaA=fit1DTransport(c,G, L, R,initGuess,'k',112,5000,5);
hold on; plot(c,G,'sk','linewidth',1,'markersize',5,'MarkerFaceColor','k');


c=[0.1018    0.1878    0.4023    0.6027    0.9455    1.4663    2.8979]*1000;
G=[1.2503    3.1084    9.0982   11.7732   20.2288   57.2088   91.7670];
sigmaA(2,:)=fit1DTransport(c,G, L, R,initGuess,'b', 112, 5000,5);
hold on; plot(c,G,'sb','linewidth',1,'markersize',5,'MarkerFaceColor','b');
 ylabel('$G$[pS]','interpreter','latex','fontsize',10)

xlabel('$c$(mM)','interpreter','latex','fontsize',10)
axis square; 

xlim([10 5000]); 
 muFitCNTP=10^(-32)*exp(sigmaA(:,2))*L^2/e^2
sigmaFitCNTP=sigmaA(:,1)*0.16*16*1000
set(gcf,'units','inches','outerposition',[0 0 3.375 3.375]) 

 % gA
 figure(3);  
 initGuess=[0.2 1 ];

 R=2; 
 L=30; 
 c=[0.0857    0.1143    0.1571    0.2143    0.2571    0.4000    0.4286    0.6571 1.3143    3.1000    5.8571]*1000; 
 G=[ 2.8703    4.2629    4.6892    6.4796    7.0195    7.6732    7.8153    9.0657 12.6181   13.0160   15.1190 ]; 
 sigmagA=fit1DTransport(c,G, L, R,initGuess,'k',112*3,5000,60);
hold on; plot(c,G,'sk','linewidth',1,'markersize',5,'MarkerFaceColor','k');
muFitgA=10^(-32)*exp(sigmagA(2))*L^2/e^2
sigmaFit=sigmagA(1)*0.16*16*1000;
axis square
xlabel('$c$(mM)','interpreter','latex','fontsize',10)
 ylabel('$G$[pS]','interpreter','latex','fontsize',10)

set(gcf,'units','inches','outerposition',[0 0 3.375 3.375])

%% Stephan-Maxwell Coupled fluxes
e=1.6e-19;

figure; 
c=[0.0313 0.0939    0.3208    0.5712    1.0016    1.8858    2.8951 ]*1000; 
G=[0.1412 1.2000    4.3235    6.8824    8.1176    8.6471    7.6412]; 


R=2; 
 L=30; 
 initGuess=[10 1]
 %sigmagA=fit1DTransport_SM(c,G, L, R,initGuess,'k',112,4000,1,0.19);
  sigmagA=fit1DTransport_SM(c,G, L, R,initGuess,'k',560,4000,1,0.038);

hold on; plot(c/1000,G,'Ok','linewidth',2,'markersize',5);
muFit_gA_SM{1}=10^(-32)*(sigmagA(:,2))*L^2/e^2; 
sigmaFit(1)=sigmagA(1)*0.16*16*1000;

 p=[ 0.3528478312181129, 1.7171741380249461
0.0978303102928628, 0.5582221253880917
0.6692339815167855, 2.4489787025437098
1.2953102020807232, 2.8685440383463146]
c=p(:,1)'*1000; 
G=p(:,2)'; 
 %sigmagA=fit1DTransport_SM(c,G, L, R,initGuess,'b',560,4000,1, 0.049);
  sigmagA=fit1DTransport_SM(c,G, L, R,initGuess,'r',112,4000,1,0.19);

hold on; plot(c/1000,G,'Or','linewidth',2,'markersize',5);
muFit_gA_SM{2}=10^(-32)*(sigmagA(:,2))*L^2/e^2; 
sigmaFit(2)=sigmagA(1)*0.16*16*1000;


p=[0.0770012891043449, 1.1760989160629656
0.3287506581694718, 2.996851680374747
0.6194709225266445, 4.1259981480472785
1.1833796321513517, 5.134484449044065
2.6356010675962747, 5.059643771447246]; 
c=p(:,1)'*1000; 
G=p(:,2)'; 
 R=2; 
 L=30; 
 initGuess=[10 1]
 sigmagA=fit1DTransport_SM(c,G, L, R,initGuess,'b',112,4000,1,0.19);
 %sigmagA=fit1DTransport_SM(c,G, L, R,initGuess,'r',560,4000,1, 0.049);
hold on; plot(c/1000,G,'Ob','linewidth',2,'markersize',5);
muFit_gA_SM{3}=10^(-32)*(sigmagA(:,2))*L^2/e^2; 
sigmaFit(3)=sigmagA(1)*0.16*16*1000;
 %%
%%
close all
initGuess=[0.028 -10 0.4];
%initGuess=[0.028 -10 0.1];

R=7.5;
L=100;
f=0.95; 

c=[3	1458.9
1.5	1045.5
1	766.4
0.6	551.1
0.4	434.7
0.2	333.3
0.1	213.8
0.05	95.6
0.02	55.5
0.01	48.2]

err=[46.0
33.9
48.8
22.1
10.1
64.2
39.1
31.6
26.6
24.1]; 
%(c,G, L, R,initGuess,'r',112,4000,1,0.19);
sigma_A(1,:)=fit1DTransport_SM(c(3:end,1)'*1000,c(3:end,2)'/1000,L,R,initGuess,'b',7,5000, 10, f);
hold on; errorbar(c(:,1),c(:,2)/1000,err/1000,'Ob','linewidth',2,'markersize',5);


c(:,2)=[538.3
401.1
241.8
193.5
142.4
90.8
50.0
22.6
14.0
10.8
];
err=[153.2
68.7
72.4
20.2
13.5
6.8
7.8
7.8
9.3
7.1
];

sigma_A(2,:)=fit1DTransport_SM(c(2:end,1)'*1000,c(2:end,2)'/1000,L,R,initGuess,'k',7,5000, 10, f);
hold on; errorbar(c(:,1),c(:,2)/1000,err/1000,'Ok','linewidth',2,'markersize',5);


%sigma_A(2,:)=fit1DTransport(c(2:end,1)'*1000,c(2:end,2)'/1000,L,R,initGuess,'k',7,10, 3000);
%hold on; errorbar(c(:,1)*1000,c(:,2)/1000,err/1000,'Ok','linewidth',2,'markersize',5);
%xlim([10 1000]); 
muFit{1}=10^(-29)*exp(sigma_A(:,2))*L^2/e^2; 

muFit{1}=10^(-29)*(sigma_A(:,2))*L^2/e;
 muFit{1}/(1e-8)
 sigma_A(:,1)*0.16*16*1000
 axis square
 ylabel('$G$[nS]','interpreter','latex','fontsize',15)
xlabel('$c$(M)','interpreter','latex')
xlim([0.01 1]); 
figure(100); 
hold on; errorbar([0.02 0.05 0.1 0.6 1] ,[9 5 4 3 2],[4 2 2 2 1],'Ok','linewidth',2,'markersize',5);
 ylabel('cation/anion','interpreter','latex','fontsize',15)
xlabel('$c$(M)','interpreter','latex')


%%
close all;
f=0.076; 
R=3.4;
L=106;
epsPore=280; 
c=[0.0053    0.0101    0.0205    0.0240    0.0498    0.0995    0.1989    0.5957 1.0017    2.9315]*1000;
G=[4.9434    6.9112   10.6199   16.1791   22.4257   32.4485   38.5319   49.4345  66.2057   69.7082];
sigma_A(1,:)=fit1DTransport_SM(c,G/1000,L,R,initGuess,'k',epsPore,5000, 1, f);
%hold on; errorbar(c(:,1),c(:,2)/1000,err/1000,'Ok','linewidth',2,'markersize',5);
hold on; plot(c/1000,G/1000,'sk','linewidth',2,'markersize',5,'MarkerFaceColor','k');


c=[0.1018    0.1878    0.4023    0.6027    0.9455    1.4663    2.8979]*1000;
G=[1.2503    3.1084    9.0982   11.7732   20.2288   57.2088   91.7670];
sigma_A(2,:)=fit1DTransport_SM(c,G/1000,L,R,initGuess,'k',epsPore,5000, 1, f);
hold on; plot(c/1000,G/1000,'sb','linewidth',2,'markersize',5,'MarkerFaceColor','b');
 ylabel('$G$[pS]','interpreter','latex','fontsize',10)

xlabel('$c$(M)','interpreter','latex','fontsize',10)
axis square; 

muFit{1}=10^(-29)*(sigma_A(:,2))*L^2/e;
 muFit{1}/(1e-8)
 sigma_A(:,1)*0.16*16*1000
 
 figure(100); hold on; plot([0.6 1 3], [184 80 10],'Or'); ylim([0 200]); xlim([0.3 3.3]); 
  ylabel('cation/anion','interpreter','latex','fontsize',15)
xlabel('$c$(M)','interpreter','latex')

%%
close all

%figure(1); subplot(1,4,1);
R=35; 
f=1;
L=30000;
clear sigma_A
initGuess=[0.028 -10 0.1];
epsPore=7; 
%initGuess=[ 0.0134  -8.7    0.1];
c=[0.0088    0.0193    0.0481    0.1355    0.3458    0.9632]*1000; 
G=[0.0104    0.0179    0.0259    0.0452    0.0863    0.2148]; 
sigma_A(1,:)=fit1DTransport_SM(c,G,L,R,initGuess,'k',epsPore,5000, 1, f);
hold on; plot(c/1000,G,'Ok','linewidth',2,'markersize',5);

%initGuess=[ 0.0134  -8.7    1];
c=[0.0026    0.0070    0.0188    0.0493    0.1338    0.3373    0.9753]*1000; 
G=[ 0.0201    0.0203    0.0290    0.0414    0.0732    0.1230    0.3142]; 
sigma_A(2,:)=fit1DTransport_SM(c,G,L,R,initGuess,'k',epsPore,5000, 1, f);
hold on; plot(c/1000,G,'Or','linewidth',2,'markersize',5);

initGuess=[ 0.0134  -8.7    0.1];

c=[    0.0026    0.0067    0.0188    0.0487    0.1338    0.3415    0.9395]*1000;
G=[ 0.1776    0.2260    0.2987    0.5701    1.1891    1.4381    2.1303];
sigma_A(3,:)=fit1DTransport_SM(c,G,L,R,initGuess,'k',epsPore,5000, 1, f);
hold on; plot(c/1000,G,'Og','linewidth',2,'markersize',5);

c=[0.0026    0.0188    0.0475    0.1355    0.3458    0.9513]*1000;
G=[ 0.9112    1.7393    2.2986    2.8154    3.3621    4.2237]; 
sigma_A(4,:)=fit1DTransport_SM(c,G,L,R,initGuess,'k',epsPore,5000, 1, f);
hold on; plot(c/1000,G,'Ob','linewidth',2,'markersize',5);




 muFit{1}=10^(-29)*(sigma_A(:,2))*L^2/e;
 muFit{1}/(1e-8)
 sigma_A(:,1)*0.16*16*1000
 sigmaFit=sigma_A(:,1)*0.16*16*1000; 