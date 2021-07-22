%% Code for generating all figures in N. Driscoll et al. Sci Trans Med. 2021.
% Updated 7/22/21
% Nicolette Driscoll, Vitale Lab 

%% Colors
% defines colors used in many plots, though this data is saved along with
% .mat files of each dataset, so this does not need to be run 
colors=[182,37,169;255,168,52;44,117,173;208,246,50;170,170,170;75,75,75]./255; 
cmap=colors;
%% Fig 1A: BODE plots for MXtex elecs of different sizes and Pt 2mm - in SALINE
clear all; load('MXtrode_allsizes_Pt2mm_EISinPBS.mat')
%|Z| plot 
figure()
[hl,hp]=boundedline(f, MX3mmmeanZnew(1:51), MX3mmstdZnew(1:51)', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl.MarkerSize=16; hl.Color=[cmap(1,:)];hp.FaceColor=[cmap(1,:)];
hold on 
[hl2,hp2]=boundedline(f, MX2mmmeanZnew, MX2mmstdZnew', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl2.MarkerSize=16;hl2.Color=[cmap(2,:)];hp2.FaceColor=[cmap(2,:)];
[hl3,hp3]=boundedline(f, MX1mmmeanZnew, MX1mmstdZnew', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl3.MarkerSize=16;hl3.Color=[cmap(3,:)];hp3.FaceColor=[cmap(3,:)];
[hl5,hp5]=boundedline(f, MX500ummeanZnew, MX500umstdZnew', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl5.MarkerSize=16;hl5.Color=[cmap(4,:)];hp5.FaceColor=[cmap(4,:)];
[hl6,hp6]=boundedline(f, Pt_2mm_meanZ, Pt_2mm_stdZ', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl6.MarkerSize=16;hl6.Color=[cmap(5,:)];hp6.FaceColor=[cmap(5,:)];
legend({'MX 3 mm','MX 2 mm','MX 1 mm','MX 500 \mum','Pt 2 mm'},'Fontsize',12);legend boxoff 
ylabel('Impedance, (\Omega)','Fontsize',12); 
xlabel('Frequency (Hz)','Fontsize',12); 
% title('MX textile and Pt size comparisons, |Z|')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,470,350])

% PHASE plot 
figure()
[hl,hp]=boundedline(f, MX3mmmeanpnew(1:51), MX3mmstdpnew(1:51)', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl.MarkerSize=12;hl.Color=[cmap(1,:)];hp.FaceColor=[cmap(1,:)];
hold on 
[hl2,hp2]=boundedline(f, MX2mmmeanpnew, MX2mmstdpnew', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl2.MarkerSize=12;hl2.Color=[cmap(2,:)];hp2.FaceColor=[cmap(2,:)];
[hl3,hp3]=boundedline(f, MX1mmmeanpnew, MX1mmstdpnew', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl3.MarkerSize=12;hl3.Color=[cmap(3,:)];hp3.FaceColor=[cmap(3,:)];
[hl5,hp5]=boundedline(f, MX500ummeanpnew, MX500umstdpnew', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl5.MarkerSize=12;hl5.Color=[cmap(4,:)];hp5.FaceColor=[cmap(4,:)];
[hl6,hp6]=boundedline(f, Pt_2mm_meanp, Pt_2mm_stdp', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl6.MarkerSize=12;hl6.Color=[cmap(5,:)];hp6.FaceColor=[cmap(5,:)];
% legend({'MX 3 mm','MX 2 mm','MX 1 mm','MX 500 \mum','Pt 2.3 mm'},'Fontsize',12,'location','southeast');legend boxoff 
ylabel('Phase, (degrees)','Fontsize',12); 
xlabel('Frequency (Hz)','Fontsize',12); 
% title('MX textile and Pt size comparisons, Phase')
ylim([-90 0])
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[900,400,470,350])
%% Fig 1D: BODE Plots for 3mm pillar vs. non-pillar MXtrodes - ON SKIN
clear all; load('MXtrode_pillar_and_nonpillar_skinEIS.mat');
% |Z| Plot
figure() 
[hl,hp]=boundedline(f, pillmeanZ, pillstdZ', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl.MarkerSize=16; hl.Color=[cmap(1,:)];hp.FaceColor=[cmap(1,:)];
hold on 
[hl2,hp2]=boundedline(f, meanZ, stdZ', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl2.MarkerSize=16; hl2.Color=[cmap(3,:)];hp2.FaceColor=[cmap(3,:)];
ylabel('Impedance, (\Omega)','Fontsize',12); 
xlabel('Frequency (Hz)','Fontsize',12); 
% title('|Z| on skin')
legend({'MX Pillar 3 mm','MX Planar 3 mm'},'Fontsize',12);legend boxoff 
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 16; 
ylim([100 1e5])
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,470,350])

% PHASE plot
figure() 
[hl,hp]=boundedline(f, pillmeanp, pillstdp', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl.MarkerSize=16; hl.Color=[cmap(1,:)];hp.FaceColor=[cmap(1,:)];
hold on 
[hl2,hp2]=boundedline(f, meanp, stdp', '.', 'alpha', 'transparency', 0.08, 'nan', 'remove');
hl2.MarkerSize=16; hl2.Color=[cmap(3,:)];hp2.FaceColor=[cmap(3,:)];
ylabel('Phase, (degrees)','Fontsize',12); 
xlabel('Frequency (Hz)','Fontsize',12); 
% title('Phase on skin')
legend({'Pillar','Planar'},'Fontsize',12,'location','southeast');legend boxoff 
set(gca, 'XScale', 'log')
ylim([-90 0])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[900,400,470,350])

%% Fig 1B: CV plot - MXtrode & Pt clinical 2 mm - intersection window: -0.6 to +0.6 V 
clear all; load('MXtrode_allsizes_Pt2mm_CVsinPBS.mat') 
d=0.2; %electrode diameter is 2 mm = 0.2 cm [cm]
gsa=pi*((d/2)^2); %geometric surface area [cm^2]
v=Intwin_v; % voltage vector corresponding to intersection window: -0.6 to +0.6 V 

figure() % MXene 2mm and Pt 2mm overlay 
plot(v,(MX2_Intwin_curm./gsa)*1000,'.','MarkerFaceColor',[cmap(1,:)],'MarkerEdgeColor',[cmap(1,:)]) % mean current for 3, 3mm MXene elecs
hold on 
plot(v,(Pt_2mm_Intwin_curm./gsa)*1000,'.','MarkerFaceColor',[cmap(5,:)],'MarkerEdgeColor',[cmap(5,:)]) % mean current for 3, 2.3mm Pt clinical elecs
ylabel('Current Density (mA/cm^2)','Fontsize',12); 
xlabel('Voltage (V)','Fontsize',12); 
% title('CV in PBS')
legend({'MX 2 mm','Pt 2 mm'},'location','northwest','Fontsize',12);legend boxoff 
xlim([-0.7 0.7])
ylim([-11 13])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,430,350])

figure() % Pt only plot 
plot(v,(Pt_2mm_Intwin_cur(:,1)./gsa)*1000,'.','MarkerFaceColor',[cmap(5,:)],'MarkerEdgeColor',[cmap(5,:)]) % mean current for 3, 2.3mm Pt clinical elecs
ylabel('Current Density (mA/cm^2)','Fontsize',12); 
xlabel('Voltage (V)','Fontsize',12); 
% title('CV in PBS')
xlim([-0.7 0.7])
ylim([-0.5 0.2])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[900,400,430,350])

%% Supp fig 7A,B - probing CV voltage limits for MXene  
clear all; load('MXtrode_CV_voltagelimittest.mat') 
d=0.3; %electrode diameter is 3 mm = 0.3 cm [cm]
gsa=pi*((d/2)^2); %geometric surface area [cm^2]

figure()%negative limit
% plot(neglim7(:,1),(neglim7(:,2)./gsa)*1000,'.','MarkerFaceColor',[cmap(6,:)],'MarkerEdgeColor',[cmap(6,:)]) % mean current for 3, 3mm MXene elecs
% hold on 
plot(neglim6(:,1),(neglim6(:,2)./gsa)*1000,'-','Color',[cmap(1,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
hold on 
plot(neglim5(:,1),(neglim5(:,2)./gsa)*1000,'-','Color',[cmap(2,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
plot(neglim4(:,1),(neglim4(:,2)./gsa)*1000,'-','Color',[cmap(3,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
plot(neglim3(:,1),(neglim3(:,2)./gsa)*1000,'-','Color',[cmap(4,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
plot(neglim2(:,1),(neglim2(:,2)./gsa)*1000,'-','Color',[cmap(5,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
plot(neglim1(:,1),(neglim1(:,2)./gsa)*1000,'-','Color',[cmap(6,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
ylabel('Current Density (mA/cm^2)','Fontsize',12); 
xlabel('Voltage (V)','Fontsize',12); 
% title('CV in PBS')
xlim([-2 0.1])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,450,350])

figure()%positive limit
plot(poslim6(:,1),(poslim6(:,2)./gsa)*1000,'-','Color',[cmap(1,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
hold on 
plot(poslim5(:,1),(poslim5(:,2)./gsa)*1000,'-','Color',[cmap(2,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
plot(poslim4(:,1),(poslim4(:,2)./gsa)*1000,'-','Color',[cmap(3,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
plot(poslim3(:,1),(poslim3(:,2)./gsa)*1000,'-','Color',[cmap(4,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
plot(poslim2(:,1),(poslim2(:,2)./gsa)*1000,'-','Color',[cmap(5,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
plot(poslim1(:,1),(poslim1(:,2)./gsa)*1000,'-','Color',[cmap(6,:)],'Linewidth',2.5) % mean current for 3, 3mm MXene elecs
ylabel('Current Density (mA/cm^2)','Fontsize',12); 
xlabel('Voltage (V)','Fontsize',12); 
% title('CV in PBS')
xlim([-0.6 0.9])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,450,350])

%% Supp fig S7D,E - CVs for each MXtrode size overlaid - in MXene and MX-Pt intersection windows 
clear all; load('MXtrode_allsizes_Pt2mm_CVsinPBS.mat')
d3=0.3; %electrode diameter is 3 mm = 0.3 cm [cm]
gsa3=pi*((d3/2)^2); %geometric surface area [cm^2]
d2=0.2; %electrode diameter is 2 mm = 0.2 cm [cm]
gsa2=pi*((d2/2)^2); %geometric surface area [cm^2]
d1=0.1; %electrode diameter is 1 mm = 0.1 cm [cm]
gsa1=pi*((d1/2)^2); %geometric surface area [cm^2]
d5=0.05; %electrode diameter is 500 um = 0.05 cm [cm]
gsa5=pi*((d5/2)^2); %geometric surface area [cm^2]
v=Intwin_v; % voltage corresponding to intersection window: -0.6 - +0.6 V 
v2=MXwin_v; % voltage corresponding to MXene window: -1.7 - +0.6 V 

% CV, Current Density - all electrodes, Intersection window (Supp fig S7E) 
figure() 
hold on
plot(v,(MX3_Intwin_curm./gsa3)*1000,'-','Color',[cmap(1,:)],'Linewidth',2.5)
plot(v,(MX2_Intwin_curm./gsa2)*1000,'-','Color',[cmap(2,:)],'Linewidth',2.5)
plot(v,(MX1_Intwin_curm(:,1)./gsa1)*1000,'-','Color',[cmap(3,:)],'Linewidth',2.5)
plot(v,(MX500_Intwin_curm(:,1)./gsa5)*1000,'-','Color',[cmap(4,:)],'Linewidth',2.5)
plot(v,(Pt_2mm_Intwin_curm./gsa2)*1000,'-','Color',[cmap(5,:)],'Linewidth',2.5)
ylabel('Current Density (mA/cm^2)','Fontsize',12); 
xlabel('Voltage (V)','Fontsize',12); 
% title('Overlaid CVs')
legend({'MX 3 mm','MX 2 mm','MX 1 mm','MX 500 \mum','Pt 2 mm'},'location','northwest','Fontsize',12);legend boxoff 
xlim([-0.7 0.7])
ylim([-25 45])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[900,400,430,350])

% CV, Current Density - MXtrodes only, MXene window (Supp fig S7D)
figure() 
hold on
plot(v2,(MX3_MXwin_curm./gsa3)*1000,'-','Color',[cmap(1,:)],'Linewidth',2.5)
plot(v2,(MX2_MXwin_curm./gsa2)*1000,'-','Color',[cmap(2,:)],'Linewidth',2.5)
plot(v2,(MX1_MXwin_curm(:,1)./gsa1)*1000,'-','Color',[cmap(3,:)],'Linewidth',2.5)
plot(v2,(MX500_MXwin_curm(:,1)./gsa5)*1000,'-','Color',[cmap(4,:)],'Linewidth',2.5)
ylabel('Current Density (mA/cm^2)','Fontsize',12); 
xlabel('Voltage (V)','Fontsize',12); 
% title('Overlaid CVs')
legend({'MX 3 mm','MX 2 mm','MX 1 mm','MX 500 \mum'},'location','northwest','Fontsize',12);legend boxoff 
xlim([-1.8 0.7])
ylim([-50 80])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[900,400,430,350])
%% Caclulating CSCs - time integral of cathodal current density (area under time-current density curve)
clear all; load('MXtrode_allsizes_Pt2mm_CVsinPBS.mat')
% change the following parameters to define which electrode and voltage
% window to calculate CSC for: 
d=0.2; %input electrode diameter in [cm] 
cur=Pt_2mm_Intwin_cur; % choose CVs corresponding to desired electrode and voltage window: 
% MX1 = MXtrode 1 mm, MX2= MXtrode 2 mm, etc.
% Intwin = -0.6 to +06 V, MXwin= -1.7 to 0.6 V, Ptwin = -0.6 to +0.8 V 
t=Intwin_t; % choose the time vector corresponding to the data of interest 
% Ptwin_t = time vector for -0.6 to +0.8 V, MXwin_t = time vector for -1.7
% to +0.6 V, Inttwin_t = time vector for -0.6 to +0.6 V 

gsa=pi*((d/2)^2); %geometric surface area [cm^2]

CSC=[];
for num=1:size(cur,2) % number of electrodes in the sample 
inds=find(cur(:,num)<0); %look only at the cathodal current 
CSC(num)=abs(trapz(t(inds),((cur(inds,num)).*1000)))/gsa; % units are [mC cm^-2]
end
mean(CSC)
std(CSC)

%% Fig 1C: Voltage Transient plots for 2 mm MXtrode and 2 mm Pt electrodes 
clear all; load('MXtrode_allsizes_Pt2mm_voltagetransient_data.mat')
start=8624; %start index for snippet
stop=9524; %stop index for snippet

figure() % MXene 2mm figure 
hold on
for i=1:size(MX2mm_Vnew,2)
plot(t(start:stop).*1000, MX2mm_Vnew(start:stop,i), '-','Color',cmap(i,:),'linewidth',1.5);
end
xlabel('Time (ms)','Fontsize',12); 
ylabel('Voltage (V)','Fontsize',12); 
title('MX 2mm')
legend({'5 mA','4 mA','3 mA','2 mA','1 mA'},'Fontsize',12,'location','southeast');legend boxoff 
xlim([29 31])
ylim([-3 3])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,455,350])

figure() % Pt 2mm figure 
hold on
for i=1:size(Pt_Vnew,2)
plot(t(start:stop).*1000, Pt_Vnew(start:stop,i), '-','Color',cmap(i,:),'linewidth',1.5);
end
xlabel('Time (ms)','Fontsize',12); 
ylabel('Voltage (V)','Fontsize',12); 
title('Pt 2mm')
legend({'5 mA','4 mA','3 mA','2 mA','1 mA'},'Fontsize',12,'location','southeast');legend boxoff 
xlim([29 31])
ylim([-3 3])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,455,350])

%% Supp fig S7 G-J: Voltage Transient plots for all MXtrode sizes 
clear all; load('MXtrode_allsizes_Pt2mm_voltagetransient_data.mat') 
start=8624; %start index for snippet
stop=9524; %stop index for snippet
figure() % MXene 3mm figure 
hold on
for i=1:size(MX3mm_V,2)
plot(t(start:stop).*1000, MX3mm_Vnew(start:stop,i), '-','Color',cmap(i,:),'linewidth',1.5);
end
xlabel('Time (ms)','Fontsize',12); 
ylabel('Voltage (V)','Fontsize',12); 
% title('MX 3mm')
text(29.2,3,'MX 3mm','Fontsize',12,'Fontweight','bold')
legend({'5 mA','4 mA','3 mA','2 mA','1 mA'},'Fontsize',12,'location','southeast');legend boxoff 
xlim([29 31])
ylim([-4 4])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,410,350])

figure() % MXene 2mm figure 
hold on
for i=1:size(MX2mm_Vnew,2)
plot(t(start:stop).*1000, MX2mm_Vnew(start:stop,i), '-','Color',cmap(i,:),'linewidth',1.5);
end
xlabel('Time (ms)','Fontsize',12); 
ylabel('Voltage (V)','Fontsize',12); 
% title('MX 2mm')
text(29.2,3,'MX 2mm','Fontsize',12,'Fontweight','bold')
legend({'5 mA','4 mA','3 mA','2 mA','1 mA'},'Fontsize',12,'location','southeast');legend boxoff 
xlim([29 31])
ylim([-4 4])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[600,400,410,350])

figure() % MXene 1mm figure 
hold on
for i=1:size(MX1mm_V,2)
plot(t(start:stop).*1000, MX1mm_V(start:stop,i), '-','Color',cmap(i,:),'linewidth',1.5);
end
xlabel('Time (ms)','Fontsize',12); 
ylabel('Voltage (V)','Fontsize',12); 
% title('MX 1mm')
text(29.2,3,'MX 1mm','Fontsize',12,'Fontweight','bold')
legend({'5 mA','4 mA','3 mA','2 mA','1 mA'},'Fontsize',12,'location','southeast');legend boxoff 
xlim([29 31])
ylim([-4 4])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[800,400,410,350])

figure() % MXene 500um figure 
hold on
for i=1:size(MX500um_V,2)
plot(t(start:stop).*1000, MX500um_V(start:stop,i), '-','Color',cmap(i,:),'linewidth',1.5);
end
xlabel('Time (ms)','Fontsize',12); 
ylabel('Voltage (V)','Fontsize',12); 
% title('MX 500um')
text(29.2,3,'MX 500um','Fontsize',12,'Fontweight','bold')
legend({'2 mA','1.5 mA','1.2 mA','1 mA','800 uA','600 uA'},'Fontsize',12,'location','southeast');legend boxoff 
xlim([29 31])
ylim([-4 4])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[1000,400,410,350])

%% Supp fig S7F,K - CSC and CIC scaling relations
clear all; load('CSC_CIC_scalingdata.mat')
[f,gof]=fit(diam,CSC(:,1),'b*x^m'); %CSC fitting
[f2,gof2]=fit(diam,CIC(:,1),'b*x^m'); %CIC fitting

figure() %CSC scaling 
hold on 
errorbar(diam,CSC(:,1),CSC(:,2),'.','Color',[cmap(5,:)],'Linewidth',1.1,'MarkerSize',22,'MarkerEdgeColor',[cmap(1,:)],'MarkerFaceColor',[cmap(1,:)]) 
plot(f,'k') % fit for CSC in MXene window
legend off
ylabel('CSC (mC cm^{-2})','Fontsize',12); 
xlabel('Diameter (mm)','Fontsize',12); 
xlim([0 4])
ylim([0 1700])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[900,400,430,350])

figure() %CIC scaling 
hold on 
errorbar(diam,CIC(:,1),CIC(:,2),'.','Color',[cmap(5,:)],'Linewidth',1.1,'MarkerSize',22,'MarkerEdgeColor',[cmap(1,:)],'MarkerFaceColor',[cmap(1,:)]) 
plot(f2,'k') 
legend off
ylabel('CIC (mC cm^{-2})','Fontsize',12); 
xlabel('Diameter (mm)','Fontsize',12); 
xlim([0 4])
ylim([0 3])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[900,400,400,350])

%% Supp fig S6 - DC conductivity of textiles with different inks 
clear all; load('MX_rGO_PEDOT_DCconductivity.mat')
figure() % MXene 
plot(length,MX_crwmean,'.','Color',cmap(1,:),'markersize',18)
ylabel('Resistance (\Omega)'); 
xlabel('Length (cm)')
% title('MXene') 
hold on
    [c,s] = polyfit(length,MX_crwmean,1);
    Rsq=1-(s.normr/norm(MX_crwmean - mean(MX_crwmean)))^2;
    % Display evaluated equation y = m*x + b
    txt1=['y = ' num2str(c(1)) ' x + ' num2str(c(2))'']; 
    txt2=['R^{2} = ' num2str(Rsq)];
    text(2, 80, txt1,'Fontsize',12);
    text(2, 72, txt2,'Fontsize',12);
    % Evaluate fit equation using polyval
    y_est1 = polyval(c,length);
    % Add trend line to plot
    plot(length,y_est1,'k--','LineWidth',1)
xlim([0 20])
ylim([10 85])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[1000,400,400,350])

figure() % PEDOT:PSS 
plot(length,pedotmean,'.','Color',cmap(2,:),'markersize',18)
ylabel('Resistance (k\Omega)'); 
xlabel('Length (cm)')
% title('PEDOT') 
hold on
    [c,s] = polyfit(length,pedotmean,1);
    Rsq=1-(s.normr/norm(pedotmean - mean(pedotmean)))^2;
    % Display evaluated equation y = m*x + b
    txt1=['y = ' num2str(c(1)) ' x + ' num2str(c(2))'']; 
    txt2=['R^{2} = ' num2str(Rsq)];
    text(2, 28, txt1,'Fontsize',12);
    text(2, 25, txt2,'Fontsize',12);
    % Evaluate fit equation using polyval
    y_est2 = polyval(c,length);
    % Add trend line to plot
    plot(length,y_est2,'k--','LineWidth',1)
xlim([0 20])
% ylim([10 85])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[1000,400,400,350])

figure() % rGO
plot(length(1:4),rGOmean(1:4),'.','Color',cmap(3,:),'markersize',18)
ylabel('Resistance (M\Omega)'); 
xlabel('Length (cm)')
% title('rGO') 
hold on
    [c,s] = polyfit(length(1:4),rGOmean(1:4),1);
    Rsq=1-(s.normr/norm(rGOmean(1:4) - mean(rGOmean(1:4))))^2;
    % Display evaluated equation y = m*x + b
    txt1=['y = ' num2str(c(1)) ' x + ' num2str(c(2))'']; 
    txt2=['R^{2} = ' num2str(Rsq)];
    text(2, 28, txt1,'Fontsize',12);
    text(2, 25, txt2,'Fontsize',12);
    % Evaluate fit equation using polyval
    y_est3 = polyval(c,length(1:4));
    % Add trend line to plot
    plot(length(1:4),y_est3,'k--','LineWidth',1)
xlim([0 20])
ylim([0 30])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[1000,400,400,350])

figure() %combined plot to show magnitude differences 
hold on 
plot(length,MX_crwmean,'.','Color',cmap(1,:),'markersize',18)
plot(length,pedotmean.*1e3,'.','Color',cmap(2,:),'markersize',18)
plot(length(1:4),rGOmean(1:4).*1e6,'.','Color',cmap(3,:),'markersize',18)
plot(length,y_est1,'k--','LineWidth',1)
plot(length,y_est2.*1e3,'k--','LineWidth',1)
plot(length(1:4),y_est3.*1e6,'k--','LineWidth',1)
ylabel('Resistance (\Omega)'); 
xlabel('Length (cm)')
legend({'MXene','PEDOT:PSS','rGO'},'Fontsize',12);legend boxoff 
xlim([0 20])
ylim([1 30e6])
ax = gca;
ax.FontSize = 16; 
box on
set(gca, 'YScale', 'log')
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[1000,400,400,350])


%% EEG plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 3C: Impedance magnitude map
clear all; load('EEG_initial_occipital.mat')
pos = [2,0.900000000000000,1,1;1,2,1,1;1,3.50000000000000,1,1;2,4.60000000000000,1,1;3.50000000000000,4.60000000000000,1,1;4.50000000000000,3.50000000000000,1,1;4.50000000000000,2,1,1;3.50000000000000,0.900000000000000,1,1;2.75000000000000,2.75000000000000,1,1];
for i=1:length(amplifier_channels)
    Z(i,:)=amplifier_channels(i).electrode_impedance_magnitude./1000;
end
order=[5,6,7,8,1,2,3,4,9];
figure()
colormap(flip(cbrewer('seq', 'YlGnBu', 100)))
c=[Z;NaN]';
c2=c(order);
scatter(pos(:,1)',pos(:,2)',3000,c2,'filled')
h= colorbar;     
ylabel(h,'1 kHz |Z| (k\Omega)','FontSize',10,'FontWeight','bold')
caxis([1 5])
ylim([0 5.5]);xlim([0 5.5]);
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Electrode Impedances')
%% Fig 3D,E: Eyes open occipital EEG plots
clear all; load('MXEEG_EOdata.mat')
% 60 Hz notch filter - apply to raw data 
dat_filt=[];
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
               'DesignMethod','butter','SampleRate',fs);
for ch=1:size(amplifier_data,1)
    dat_filt(ch,:)= filtfilt(d,amplifier_data(ch,:));
end
% bandpass filter - apply to notch filtered data 
dat_filt2=[];
[b,a]=butter(2,[0.1,100]/(fs/2),'bandpass');
for ch=1:size(amplifier_data,1)
    dat_filt2(ch,:)=filtfilt(b,a,dat_filt(ch,:));
end

% Fig 3D, left: plot of filtered data from all channels - Eyes Open 
% NOTE: Ag/AgCl cup electrode is ch9
figure()
for i=1:size(amplifier_data,1)
plot(t_amplifier,dat_filt2(i,:)-(100*(i-1)),'k')
hold on 
end
title('Eyes Open') 
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel('Signal (uV)','FontSize',12,'FontWeight','bold')
xlim([51.4 56.4])

% Fig 3E, top: Spectrogram, eyes open
MX=2; % choose which MXene channel to plot spectrogram for
cup=9; % gelled Ag/AgCl cup electrode is Ch9
frq = [1 58]; % frequency range (Hz)
x=amplifier_data(MX,:);
figure() 
colormap(flip(cbrewer('div', 'RdYlBu', 100)))
spectrogram(x,1*fix(fs),0.9*fix(fs),logspace(log10(frq(1)),log10(frq(2)),250),fs,'yaxis');
caxis([-20 20])
ax = gca;
ax.YScale = 'log';
% title('Eyes Open - MX')
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,800,250])

%% Fig 3D,E: Eyes closed occipital EEG plots
clear all; load('MXEEG_ECdata.mat')
% 60 Hz notch filter - apply to raw data 
dat_filt=[];
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
               'DesignMethod','butter','SampleRate',fs);
for ch=1:size(amplifier_data,1)
    dat_filt(ch,:)= filtfilt(d,amplifier_data(ch,:));
end
% bandpass filter - apply to notch filtered data 
dat_filt2=[];
[b,a]=butter(2,[0.1,100]/(fs/2),'bandpass');
for ch=1:size(amplifier_data,1)
    dat_filt2(ch,:)=filtfilt(b,a,dat_filt(ch,:));
end

% Fig 3D, right: plot of filtered data from all channels - Eyes Closed 
% NOTE: Ag/AgCl cup electrode is ch9
figure()
for i=1:size(amplifier_data,1)
plot(t_amplifier,dat_filt2(i,:)-(100*(i-1)),'k')
hold on 
end
title('Eyes Closed') 
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel('Signal (uV)','FontSize',12,'FontWeight','bold')
xlim([56 61])

% Fig 3E, bottom: Spectrogram, eyes closed
MX=2; % choose which MXene channel to plot spectrogram for
cup=9; % gelled Ag/AgCl cup electrode is Ch9
frq = [1 58]; % frequency range (Hz)
x=amplifier_data(MX,:);
figure() 
colormap(flip(cbrewer('div', 'RdYlBu', 100)))
spectrogram(x,1*fix(fs),0.9*fix(fs),logspace(log10(frq(1)),log10(frq(2)),250),fs,'yaxis');
caxis([-20 20])
ax = gca;
ax.YScale = 'log';
% title('Eyes Open - MX')
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,800,250])

%% Alpha Bandpower Calculations - for Supp fig S8 and Movie S1
% load('MXEEG_ECdata.mat') % load this for data from eyes closed recording
% load('MXEEG_EOdata.mat') % load this for data from eyes open recording
alphap=[]; 
band=[8 12];
winLen=1; %units are seconds
winDisp=0.5; %units are seconds
for ch=1:size(amplifier_data,1)
    alphap(ch,:)=MovingWin_bandpower(amplifier_data(ch,:),fs,band,winLen,winDisp);
end
% clean up the artifacts by thresholding the bandpower 
lim=300; % set the upper limit for the bandpower
alphap_new=alphap; % threshold the artifacts 
alphap_new(alphap_new>lim)=lim; % threshold the artifacts 
alphaEC=alphap_new;
meanalphapEC=mean(alphaEC,2);
% alphaEC=alphap_new; %save this as new variable for eyes closed state 
% alphaEO=alphap_new; %save this as new variable for eyes open state 
% meanalphapEC=mean(alphaEC,2); %mean alpha, eyes closed
% meanalphapEO=mean(alphaEO,2); %mean alpha, eyes open

%% Supp fig S8A,B: alpha bandpower maps from occipital EEG recordings
clear all; load('EEG_ECvsEOalphabandpower_values.mat')
pos = [2,0.900000000000000,1,1;1,2,1,1;1,3.50000000000000,1,1;2,4.60000000000000,1,1;3.50000000000000,4.60000000000000,1,1;4.50000000000000,3.50000000000000,1,1;4.50000000000000,2,1,1;3.50000000000000,0.900000000000000,1,1;2.75000000000000,2.75000000000000,1,1];
order=[5,6,7,8,1,2,3,4,9]; % map channels onto order of dots (in pos) 
lim=50; % set upper limit of colorbar scale 
dotsize=[1500,1500,1500,1500,1500,1500,1500,1500,5000];
    
% Supp Fig S8A: alpha power, eyes open
figure() 
colormap(flip(cbrewer('div', 'RdYlBu', 100)))
c3=[meanalphapEO]';
c4=c3(order);
hold on
scatter(pos(:,1)',pos(:,2)',dotsize,c4,'filled')
h=colorbar; caxis([15 lim])
ylabel(h,'Bandpower (dB)','FontSize',10,'FontWeight','bold')
ylim([0 5.5]);xlim([0 5.5]);
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Eyes Open, mean alpha power','FontSize',12,'FontWeight','bold')
ax = gca;
ax.FontSize = 16; 
box off
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[200,200,500,400])

% Supp Fig S8B: alpha power, eyes closed
figure() 
colormap(flip(cbrewer('div', 'RdYlBu', 100)))
c3=[meanalphapEC]';
c4=c3(order);
hold on
scatter(pos(:,1)',pos(:,2)',dotsize,c4,'filled')
h=colorbar; caxis([15 lim])
ylabel(h,'Bandpower (dB)','FontSize',10,'FontWeight','bold')
ylim([0 5.5]);xlim([0 5.5]);
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Eyes Open, mean alpha power','FontSize',12,'FontWeight','bold')
ax = gca;
ax.FontSize = 16; 
box off
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[200,200,500,400])

%% Supp fig S8 bottom panels: alpha bandpower bar plots
clear all; load('EEG_ECvsEOalphabandpower_values.mat')
EOalpha=mean(alphaEO,2);
ECalpha=mean(alphaEC,2);

ECerr=std(alphaEC,0,2);
ECerr=ECerr./sqrt(length(alphaEC));
EOerr=std(alphaEO,0,2);
EOerr=EOerr./sqrt(length(alphaEO));

figure() 
bar(EOalpha)
hold on
errorbar(EOalpha,EOerr,'.');
ylabel('Mean \alpha bandpower')
xlabel('Electrode')
title('Eyes Open')
ylim([15 55])
set(gcf,'position',[200,200,350,250])

figure() 
bar(ECalpha)
hold on
errorbar(ECalpha,ECerr,'.');
ylabel('Mean \alpha bandpower')
xlabel('Electrode')
title('Eyes Closed')
ylim([15 55])
set(gcf,'position',[200,200,350,250])

%% Statistics on the bar plots in fig. S8
% one-way ANOVA
[p,tbl] = anova1(alphaEC');
[p2,tbl2] = anova1(alphaEO');

%% EMG Plots - Anterior Policis Brevis (APB) - base of thumb 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; load('MXtrode_EMG_thumb_stim.mat') % load this data for the next few sections
%% Fig 4A: plot average evoked response over grid with detected peaks - APB
grid=[13,11,20,18,16,14,12,19,17,15,1,3,5,7,9,2,4,6,8,10]; % channel mapping 
figure() % all channels in grid plot
hold on
for i=1:length(grid); 
    ax(i)=subplot(4,5,i);
    hold on 
    plot((0:1/fs:(400+on)/fs),-squeeze(mean(segs(grid(i),:,:),2)),'k','linewidth',1)
    plot(peaks_ind2(grid(i))/fs,vals(grid(i))+500,'.r') % plotting detected peak 
end
linkaxes([ax],'xy')
set(gcf,'position',[200,200,600,400])

%% Fig 4B: latency map, interpolated, APB 
[X,Y] = meshgrid(1:size(latency_grid,2), 1:size(latency_grid,1));
% define the resolution by changing the delta from 0.01 to whatever (0.001 = more pixelation, 0.1 = less pixelation)
[X2,Y2] = meshgrid(1:0.01:size(latency_grid,2), 1:0.01:size(latency_grid,1));
latency_interp = interp2(X, Y, latency_grid, X2, Y2); 

for i=1:size(X,2)
    test(i)=find(X2(1,:)==X(1,i));
end
for i=1:size(Y,1)
    test2(i)=find(Y2(:,1)==Y(i,1));
end
elec_gridx=repelem(test(1:5),4);
elec_gridy=[test2(1:4),test2(1:4),test2(1:4),test2(1:4),test2(1:4)];

figure()
imagesc(latency_interp)
c=colorbar;
c.Label.String='Latency (ms)';
c.FontSize=14;
set(gca,'xtick',[])
set(gca,'ytick',[])
% title('Peak Latency Map')
hold on 
plot(elec_gridx,elec_gridy,'.k','markersize',16) % add electrode locations 
colormap(flip(cbrewer('div', 'RdYlBu', 100)))

%% Supp Fig S11 A: Impedance map - APB 
Z=[amplifier_channels.electrode_impedance_magnitude];
disp('Average |Z|1kHz (in kOhm):')
[mean(Z)/1000 std(Z)/1000]
Zmap=Z(grid);
Zmap2=reshape(Zmap,5,4)';

figure() 
imagesc(Zmap2./1000)
c=colorbar;
c.Label.String='1 kHz |Z| (k\Omega)';
c.FontSize=14;
set(gca,'xtick',[])
set(gca,'ytick',[])
% title('Impedance Map')
colormap(flip(cbrewer('seq', 'YlGnBu', 100)))
caxis([0 150]);

%% EMG Plots - Biceps 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; load('MXtrode_EMG_biceps_stim.mat') % load this data for the next few sections 
%% Fig 4C: plot average evoked response over grid with detected peaks - biceps 
grid=[19,17,15,13,11,9,7,5,3,1,20,18,16,14,12,10,8,6,4,2,40,38,36,34,32,30,28,26,24,22,...
    39,37,35,33,31,29,27,25,23,21]; % channel mapping 
figure() % all channels in grid plot
hold on
for i=1:length(grid); 
    ax(i)=subplot(4,10,i);
    hold on 
    plot((0:1/fs:(400+on)/fs),squeeze(mean(segs(grid(i),:,:),2)),'k','linewidth',1)
%     plot((onset_ind(grid(i))+100)/fs,squeeze(mean(segs(grid(i),:,100+onset_ind(grid(i))),2)),'r+') % plotting detected peak 
    plot((peaks_ind(grid(i))+100)/fs,squeeze(mean(segs(grid(i),:,100+peaks_ind(grid(i))),2))+500,'.r') % plotting detected peak +500,'.r') % plotting detected peak 
    plot(on/fs,4000,'.b')
end
linkaxes([ax],'xy')
set(gcf,'position',[100,100,1200,400])

%% Fig 4D: Latency map - peak response  - interpolated - biceps 
[X,Y] = meshgrid(1:size(latency_grid2,2), 1:size(latency_grid2,1));
% define the resolution by changing the delta from 0.01 to whatever (0.001 = more pixelation, 0.1 = less pixelation)
[X2,Y2] = meshgrid(1:0.01:size(latency_grid2,2), 1:0.01:size(latency_grid2,1));
latency_interp = interp2(X, Y, latency_grid2, X2, Y2); 

for i=1:size(X,2)
    test(i)=find(X2(1,:)==X(1,i));
end
for i=1:size(Y,1)
    test2(i)=find(Y2(:,1)==Y(i,1));
end
elec_gridx=repelem(test(1:10),4);
elec_gridy=[test2(1:4),test2(1:4),test2(1:4),test2(1:4),test2(1:4),test2(1:4),test2(1:4),test2(1:4),test2(1:4),test2(1:4)];

figure()
imagesc(latency_interp)
c=colorbar;
c.Label.String='Latency (ms)';
c.FontSize=14;
set(gca,'xtick',[])
set(gca,'ytick',[])
% title('Peak Latency Map')
hold on 
plot(elec_gridx,elec_gridy,'.k','markersize',16) % add electrode locations 
colormap(flip(cbrewer('div', 'RdYlBu', 100)))
set(gcf,'position',[100,100,800,300])

%% Supp fig S11 B: Impedance map - Biceps
Z=[amplifier_channels.electrode_impedance_magnitude];
grid=[19,17,15,13,11,9,7,5,3,1,20,18,16,14,12,10,8,6,4,2,40,38,36,34,32,30,28,26,24,22,...
    39,37,35,33,31,29,27,25,23,21];
disp('Average |Z|1kHz (in kOhm):')
[mean(Z)/1000 std(Z)/1000] 
Zmap=Z(grid);
Zmap2=reshape(Zmap,10,4)';

figure() 
imagesc(Zmap2./1000)
c=colorbar;
c.Label.String='1 kHz |Z| (k\Omega)';
c.FontSize=14;
set(gca,'xtick',[])
set(gca,'ytick',[])
% title('Impedance Map')
colormap(flip(cbrewer('seq', 'YlGnBu', 100)))
caxis([0 100]);

%% Supp fig S11 D: Bipolar subtraction - Biceps
clear all; load('MXtrode_EMGbiceps_resistedflexion.mat')
xlims=[22.80 22.9];%[22.54 22.62];  %[21.57 21.65];
figure() % top row of array (4th column, rotated)
hold on 
for i=1:9
    plot(t,(bip_dat(i,:)./1e3)-(2*(i-1)),'k')
end
xlabel('Time (s)')
ylabel('Signal (mV)')
xlim([xlims])
title('Row 1')
set(gcf,'position',[200,50,230,700])

figure() % second row of array (3rd column, rotated)
hold on 
for i=10:18
    plot(t,(bip_dat(i,:)./1e3)-(2*(i-1-9)),'k')
end
xlabel('Time (s)')
ylabel('Signal (mV)')
xlim([xlims])
title('Row 2')
set(gcf,'position',[200,50,230,700])

figure() % third row of array (2nd column, rotated)
hold on 
for i=19:27
    plot(t,(bip_dat(i,:)./1e3)-(2*(i-19)),'k')
end
xlabel('Time (s)')
ylabel('Signal (mV)')
xlim([xlims])
title('Row 3')
set(gcf,'position',[200,50,230,700])

figure() % fourth row of array (1st column, rotated)
hold on 
for i=28:36
    plot(t,(bip_dat(i,:)./1e3)-(2*(i-28)),'k')
end
xlabel('Time (s)')
ylabel('Signal (mV)')
xlim([xlims])
title('Row 4')
set(gcf,'position',[200,50,230,700])

%% Cardiac ECG recording 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; load('MXtrode_ECGrecording.mat') % load this data for the next few sections 
%% Fig 5B: ECG data plots - MXene and gelled Natus
% plot MXene ECG data (Fig 5B) 
figure() 
plot(t_MX-t_MX(1)-0.6,filt_MX2(1,:)./1e3,'k') 
title('MXene','fontsize',12)
xlabel('Time (s)','Fontsize',12) 
ylabel('Signal (mV)','Fontsize',12)
ylim([-0.8 3.4])
xlim([0 10.5])
ax.FontSize = 16; 
set(gcf,'position',[200,200,600,200])

% plot Natus ECG data (Fig 5B) 
figure() 
plot(t_Natus-t_Natus(1)-50,filt_Natus2(1,:)./1e3,'b')
title('Natus','fontsize',12)
xlabel('Time (s)','fontsize',12) 
ylabel('Signal (mV)','fontsize',12)
ylim([-0.8 3.4])
xlim([0 10.5])
ax.FontSize = 16; 
set(gcf,'position',[200,200,600,200])

%% Fig 5C: overlaid average ECG waveforms - MXene and Natus
figure() % overlaid waveforms
hold on 
plot((0:1/fs:(start+stop)/fs),(mean(segs_MX,1))./1000,'k','linewidth',1)
plot((0:1/fs:(start+stop)/fs),(mean(segs_N,1))./1000,'b','linewidth',1)
ylim([-1 3])
legend({'MXene','Natus'},'fontsize',12)
xlabel('Time (s)','fontsize',12)
ylabel('Signal (mV)','fontsize',12)
ax.FontSize = 16; 
set(gcf,'position',[200,200,350,330])
box on 
%% SNR calculations of the ECG waveforms 
% R peak amplitudes: 
MX_Rpeak=max(mean(segs_MX)) % mean R peak amp for MXene 
a=find(mean(segs_MX)==max(mean(segs_MX))); % find index where the peak occurs 
std(segs_MX(:,a)) % standard deviation of R peak amp for MXene 

NR_peak=max(mean(segs_N)) % mean R peak amp for Natus 
a=find(mean(segs_N)==max(mean(segs_N))); % find index where the peak occurs 
std(segs_N(:,a)) % standard deviation of R peak amp for Natus 

% baseline noise values - choose time between two heartbeats and compute
% the RMS of the signal there
t1=[73.5,74.3]; t2=[74.9, 75.7]; t3=[76.34,77.09]; % a few time segments of baseline - MXene 
inds=[min(find((t_MX-t1(1))>0)),min(find((t_MX-t1(2))>0));...
    min(find((t_MX-t2(1))>0)),min(find((t_MX-t2(2))>0));...
    min(find((t_MX-t3(1))>0)),min(find((t_MX-t3(2))>0))];

for i=1:3
    rmsMX(i)=rms(amplifier_data_MX(:,[inds(i,1):inds(i,2)]));
end
mean_rmsMX=mean(rmsMX);
SNR_MX=MX_Rpeak/mean_rmsMX

%Natus
t1=[28.08,28.9]; t2=[29.5, 30.25]; t3=[30.84,31.53]; % a few time segments of baseline - MXene 
inds=[min(find((t_Natus-t1(1))>0)),min(find((t_Natus-t1(2))>0));...
    min(find((t_Natus-t2(1))>0)),min(find((t_Natus-t2(2))>0));...
    min(find((t_Natus-t3(1))>0)),min(find((t_Natus-t3(2))>0))];

for i=1:3
    rmsN(i)=rms(amplifier_data_Natus(:,[inds(i,1):inds(i,2)]));
end
mean_rmsN=mean(rmsN);
SNR_N=NR_peak/mean_rmsN

%% ECoG Recordings - pig 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; load('MXtrode_ECoG_pig_recording.mat') % load this data for the next few sections 
%% Fig 6B: ECoG segment recorded in pig
% plot segment of raw data, vertically staggered
figure()
hold on 
for i=1:size(Grid,1)
    plot((tvector-tvector(1))./1e6,(Grid(i,:)-7000*(i-1))./1000,'k')
end
ylabel('Signal (mV)')
xlabel('Time (s)') 
title('Raw ECoG Signal')
xlim([89 94])

%% Fig 6C: PSD plot from pig ECoG recordings 
t1=13;
t2=123;
i1=(t1*fs);
i2=(t2*fs);
frq = [0.3 600]; % frequency range (Hz)
figure('Name','freq','NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
for ich = 1
    dat = Grid(ich,i1:i2);
    dat = dat-median(dat);
    [P,f] = pwelch(dat,5*fix(fs),1*fix(fs),logspace(log10(frq(1)),log10(frq(2)),500),fs); 
    ind = find(f>=frq(1) & f<=frq(2));
    f = f(ind); P = P(ind);
    plot(f,smooth(P,3),'Color','k','Tag',num2str(ich),'linewidth',1.2); hold on;
end
hold on 
set(gca,'Box','off','XScale','log','YScale','log','XLim',frq); axis tight;
xlabel('Frequency (Hz)','fontweight','bold','fontsize', 12); 
ylabel('Power (\muV^2/Hz)','fontweight','bold','fontsize', 12); axis square;

%% Fig 6D: plotting ECoG data segment spatially across array 
t1=88;
t2=93;
i1=(t1*fs);
i2=(t2*fs);
mapping=[2,3,6,5,1,4]; %maps channel numbers to subplot locations 

figure() 
hold on 
ax1=subplot(3,2,1);
plot((tvector(i1:i2)-tvector(1))./1e6,(dat_filt2(mapping(1),i1:i2))./1e3,'k');
ax2=subplot(3,2,2);
plot((tvector(i1:i2)-tvector(1))./1e6,(dat_filt2(mapping(2),i1:i2))./1e3,'k');
ax3=subplot(3,2,3);
plot((tvector(i1:i2)-tvector(1))./1e6,(dat_filt2(mapping(3),i1:i2))./1e3,'k');
ax4=subplot(3,2,4);
plot((tvector(i1:i2)-tvector(1))./1e6,(dat_filt2(mapping(4),i1:i2))./1e3,'k');
ax5=subplot(3,2,5);
plot((tvector(i1:i2)-tvector(1))./1e6,(dat_filt2(mapping(5),i1:i2))./1e3,'k');
ax6=subplot(3,2,6);
plot((tvector(i1:i2)-tvector(1))./1e6,(dat_filt2(mapping(6),i1:i2))./1e3,'k');
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'xy')
ylim([-6 6])
xlim([88.9 92.95])

%% Fig 6E: Voltage mapping spatially across array 
%choose time point to plot: 
t=89.22231; % down state 1
% t=89.79799; % up state 2
% t=90.82902; % up state 3
% t=92.65227; % down state 4
ind=min(find((tvector-tvector(1))>(t*1e6)));
for i=1:size(dat_filt2,1)
    voltage(i)=dat_filt2(i,ind);
end
voltage2=voltage(mapping);
voltage2=reshape(voltage2,2,3)';
voltage2= voltage2./1000;
% normalize values 
voltage3=(voltage2-min(voltage2(:)))./(max(voltage2(:))-min(voltage2(:)));

% plotting normalized voltages: 
[X,Y] = meshgrid(1:size(voltage3,2), 1:size(voltage3,1));
% define the resolution by changing the delta from 0.01 to whatever (0.001 = more pixelation, 0.1 = less pixelation)
[X2,Y2] = meshgrid(1:0.15:size(voltage3,2), 1:0.15:size(voltage3,1));
interp = interp2(X, Y, voltage3, X2, Y2); 

figure()
imagesc(interp)
c=colorbar;
c.Label.String='Amplitude (mV)';
set(gca,'xtick',[])
set(gca,'ytick',[])
str = sprintf('Voltage Amplitude Map, t=%f', t);
title(str)
% colormap jet 
colormap(flip(cbrewer('div', 'RdYlBu', 100)))
set(gcf,'position',[200,200,230,400])

%% Cortical Whisker stimulation - rat 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; load('MXtrode_whisker_stim.mat') % load this data for the next few sections 
%% Fig 6G: stimulation pulse trains and response, 1.4 mA stimulation pulses 
figure() 
plot(t, stim_data(3,:)./1e3,'k') 
xlabel('time (s)') 
ylabel('Current (mA)') 
title('Stim Train') 
ylim([-1.7 1.7])
xlim([2 11])
set(gcf,'position',[200,50,500,150])

% plot whisker deflection filtered data 
test=[datfilt3.*(28/5)]; % this is a conversion factor from the optical micrometer - supposedly converts to [mm] units
test(test<20) = 20;
figure() 
plot(t,test,'k') 
xlabel('time (s)') 
ylabel('Whisker deflection (mm)') 
title('Whisker deflection - 1.4 mA') 
ylim([19 25])
xlim([2 11])
set(gcf,'position',[200,50,500,150])
%% Fig 6H - mean whisker deflection
figure() 
hold on 
plot(tseg,smooth(mean1_4,200)*(30),'color',cmap(1,:),'linewidth',1.5)
plot(tseg,smooth(mean1_2,200)*(30),'color',cmap(3,:),'linewidth',1.5) 
plot(tseg,smooth(-mean1_1,200)*(30),'color',cmap(2,:),'linewidth',1.5)
plot(tseg,smooth(mean1_0,200)*(30),'color',cmap(4,:),'linewidth',1.5)
box on 
legend('1.4mA','1.2mA','1.1mA','1.0mA','location','northwest')
legend boxoff 
xlabel('Time (s)','Fontsize',12)
ylabel('Whisker deflection (A.U.)','Fontsize',12)
xlim([0.08 0.21])
set(gcf,'position',[200,200,380,280])

%% Fig 6I: average magnitude of first whisker deflection for each of the stimulation amplitudes 
a0=0;
a1=max(findpeaks(smooth(-mean1_1(1:6000),200).*30,'MinPeakDistance',1000));
a2=max(findpeaks(smooth(mean1_2(1:6000),200).*30,'MinPeakDistance',1000));
a3=max(findpeaks(smooth(-mean1_3(1:6000),200).*30,'MinPeakDistance',1000));
a4=max(findpeaks(smooth(mean1_4(1:6000),200).*30,'MinPeakDistance',1000));
a5=max(findpeaks(smooth(mean1_5(1:6000),200).*30,'MinPeakDistance',1000));

stimamps=[1,1.1,1.2,1.4];
figure()
plot(stimamps,[a0,a1,a2,a4],'.r','Markersize',25)
xlim([0.9 1.5])
% ylim([0 0.4])
xlabel('Stimulation Amplitude (mA)','Fontsize',12) 
ylabel('Whisker Response (A.U.)','Fontsize',12)
set(gcf,'position',[200,200,380,280])

%% fig S14: Human EOG eye tracking 
clear all; load('MXtrode_EOG.mat')
% fig S14B
figure()
plot(t_updown-3,eog_updown./1e3,'k')
xlabel('Time (s)')
ylabel('Signal (mV)')
title('EOG up down 60 Hz notch + 0.1-20 Hz bp filt') 
xlim([0 17])
ylim([-0.65 0.5])
set(gcf,'position',[200,200,800,300])
% fig S14D
figure() 
plot(t_leftright-4,eog_leftright./1e3,'k')
xlabel('Time (s)')
ylabel('Signal (mV)')
title('EOG left right 60 Hz notch + 0.1-20 Hz bp filt') 
xlim([0 17])
set(gcf,'position',[200,200,800,300])
ylim([-0.85 0.85])

%% fig S16D: Magnetic Susceptibility 
clear all; load('mag_susceptibility_Ti3C2_data2.mat') 
% fig S16D
% field [Oe]: 1 Oe = 1000/4pi A/m
% moment [emu]: 1 emu/cm^3 = 1000 A/m 
% vol: sample volume [cm^3] (sample mass=4.82 mg, Ti3C2 density=3.7 g/cm^3)
vol = (4.82/1000)/3.7;
inds=[10:955];
% X = M/H = [(moment/vol)/1000] / [(field*4pi/1000)] - dimensionless
chi=((moment(inds)./vol)./1000)./((field(inds).*4.*pi)./1000);

figure() 
% plot(field(inds).*1e-4,moment(inds),'.k') 
plot(field(inds).*1e-4,(moment(inds)./vol),'.k') % field in Tesla, moment in emu/cm^3
xlabel('Magnetic Field (Tesla)') 
ylabel('Magnetization (emu/cm^3)')
xlim([0 9])

figure()
plot(((field(inds).*4.*pi)./1000),((moment(inds)./vol)./1000),'.k')% field and moment in A/m
xlabel('Magnetic Field (A/m)') 
ylabel('Magnetization (A/m)')
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
xlim([0 1150])

susceptibility=mean(chi)
var=std(chi)

%% fig S1B: Generalizability to different textile substrates - conductiviy of strips
clear all; load('textile_generalizability.mat') 
types=categorical({'Technicloth','TexVantage','Cotton','Polyester'}); 
% bulk conductivity calculations
% rho=R.*(xArea/L); %resistivity - Ohm*m
% stdrho=stdR.*(xArea/L); %stdev resistivity - Ohm*m 
% dev=stdrho./rho; %error as a fraction of the resistivity 
% sigma=1./rho; %conductivity - S/m
% ssigma=dev.*sigma; %stdev conductivity - S/m --- multiply fractional error times conductivity 

figure() 
bar(types,sigma([1,2,3,5])./100)
hold on 
errorbar(types,sigma([1,2,3,5])./100,ssigma([1,2,3,5])./100,'.','Color','k','Linewidth',1.1);
ylabel('Conductivity, \sigma (S/cm)');
ax = gca;
ax.FontSize = 14; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[200,200,410,380])

%% fig S1F,G: Percolation threshold - resistance vs. MXene concentration
clear all; load('percolation_thresh.mat')
% meanR=mean(R,2);
% stdR=std(R')';
% Conc=categorical({'26','20','10','5','2','1','0.5'});
% Conc = reordercats(Conc,{'26','20','10','5','2','1','0.5'});

% bar plot - resistance - fig S1F
figure() 
bar(Conc,meanR)
hold on 
errorbar(Conc,meanR,stdR,'.','Color','k','Linewidth',1.1);
xlabel('Concentration (mg/mL)'); 
ylabel('Resistance (\Omega)');
set(gca,'FontSize',14); 
% ylim([4830-(650/2) 4830+(650/2)])
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 14; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[200,200,390,340])

% bar plot - bulk conductivity - fig S1G
A=8.55e-7; %cross-sectional area: 285um x 3mm - converted to m^2
L=0.05; %length of strip: 5cm = 0.05 m 
rho=meanR.*(A/L); %resistivity - Ohm*m
stdrho=stdR.*(A/L); %stdev resistivity - Ohm*m 
dev=stdrho./rho; %error as a % of original value for resistivity (fraction)
sigma=1./rho; %conductivity - S/m
ssigma=dev.*sigma; %stdev conductivity - S/m --- multiply fractional error times conductivity 

figure() 
bar(Conc,sigma./100)
hold on 
errorbar(Conc,sigma./100,ssigma./100,'.','Color','k','Linewidth',1.1);
xlabel('Concentration (mg/mL)'); 
ylabel('Conductivity, \sigma (S/cm)');
set(gca,'FontSize',14); 
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 14; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[200,200,390,340])

%% fig S1D: Durability of materials - bending tests
clear all; load('bending1000cycletest.mat') % located in \MXtrodes projects\Instron mechanical tests

% Plotting materials together - bending cycle # axis 
figure() 
plot(cycles,smooth((V_technicloth2./I_technicloth2),10)-(V_technicloth2(1)./I_technicloth2(1)),'k') 
hold on 
plot(cycles2,smooth((V_texvantage./I_texvantage),10)-(V_texvantage(1)./I_texvantage(1)),'b') 
plot(cycles2,smooth((V_cotton2./I_cotton2),10)-(V_cotton2(1)./I_cotton2(1)),'r') 
plot(cycles2,smooth((V_polyester1./I_polyester1),10)-(V_polyester1(1)./I_polyester1(1)),'g') 
xlabel('Cycle #','Fontsize',12) 
ylabel('R - R_0 (\Omega)','Fontsize',12) 
ax = gca; ax.FontSize = 14; 
xlim([0 1000]); set(gcf,'position',[200,200,700,250])
ylim([-5 10])
legend('TC','TV','Cotton','Polyester','Cellulose') %TC = Technicloth, TV = TexVantage
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];

%% Individual plots for different materials - resistance change during bending cycles 

% Technicloth 2 - 1000 cycles 
figure() 
plot(time,(V_technicloth2./I_technicloth2),'k') 
title('Technicloth 2, 1000 cycles')
xlabel('Time (s)','Fontsize',12) 
ylabel('Resistance, \Omega','Fontsize',12) 
ax = gca; ax.FontSize = 14; 
xlim([0 time(end)]); set(gcf,'position',[200,200,900,250])
% ylim([39.8 41.5])

% Cotton 2 - 1000 cycles 
figure() 
plot(time2,(V_cotton2./I_cotton2),'k') 
title('Cotton 2, 1000 cycles')
xlabel('Time (s)','Fontsize',12) 
ylabel('Resistance, \Omega','Fontsize',12) 
ax = gca; ax.FontSize = 14; 
xlim([0 time(end)]); set(gcf,'position',[200,200,900,250])
ylim([39.5 43.5])

% Texvantage - 1000 cycles 
figure() 
plot(time2,(V_texvantage./I_texvantage),'k') 
title('Texvantage, 1000 cycles')
xlabel('Time (s)','Fontsize',12) 
ylabel('Resistance, \Omega','Fontsize',12) 
ax = gca; ax.FontSize = 14; 
xlim([0 time(end)]); set(gcf,'position',[200,200,900,250])
% ylim([39.8 41.5])

% Polyester 1 - 1000 cycles 
figure() 
% plot(time2,(V_polyester1./I_polyester1),'k') 
plot(time2,smooth((V_polyester1./I_polyester1),8),'k') 
title('Polyester 1, 1000 cycles')
xlabel('Time (s)','Fontsize',12) 
ylabel('Resistance, \Omega','Fontsize',12) 
ax = gca; ax.FontSize = 14; 
xlim([0 time(end)]); set(gcf,'position',[200,200,900,250])
% ylim([39.8 41.5])

%% fig S1E: Durability of materials - bending tests - range and %change calcs

% Technicloth samples - 1000 cycles 
R1=[];
R1=V_technicloth2./I_technicloth2;
range1=max(R1)-min(R1)
variance1=sum((R1-mean(R1).^2))/length(R1);
change1=R1(end)-R1(1);
change2=(change1/R1(1)*100)

% Cotton samples - 1000 cycles 
R1=[];
R1=V_cotton2./I_cotton2;
range1=max(R1)-min(R1)
variance1=sum((R1-mean(R1).^2))/length(R1);
change1=R1(end)-R1(1);
change2=abs(change1/R1(1)*100)

% Texvantage - 1000 cycles 
R1=[];
R1=V_texvantage./I_texvantage;
range1=max(R1)-min(R1)
variance1=sum((R1-mean(R1).^2))/length(R1);
change1=R1(end)-R1(1);
change2=(change1/R1(1)*100)

% Polyester samples - 1000 cycles  --- use 1st sample only, second is a bad
% measurement 
R1=[];
R1=V_polyester1(1:3401)./I_polyester1(1:3401);
range1=max(R1)-min(R1)
variance1=sum((R1-mean(R1).^2))/length(R1);
change1=R1(end)-R1(1);
change2=(change1/R1(1)*100)

%% fig S4: Instron - pillar strength testing
clear all; load('Instron_pillars.mat') % location:  C:\Users\nicki\Dropbox\MXtrodes projects\Instron mechanical tests\tesion to failure - pillars 
% area=pi*(1.5^2); % area in mm^2

figure()  % stress-strain plot - lateral shear force - fig S4A
plot(ext_lateral,F_lateral./area,'LineWidth',1.5)
xlabel('Strain,\epsilon (mm)','Fontsize',12) 
ylabel('Stress, \sigma (MPa)','Fontsize',12) 
ax = gca; ax.FontSize = 14; 
set(gcf,'position',[200,200,400,300])

figure()  % stress-strain plot - vertical pull force - fig S4B
plot(ext_vertical,F_vertical./area,'LineWidth',1.5)
xlabel('Strain,\epsilon (mm)','Fontsize',12) 
ylabel('Stress, \sigma (MPa)','Fontsize',12) 
ax = gca; ax.FontSize = 14; 
set(gcf,'position',[200,200,400,300])
xlim([0 0.5])

%% Motion artifacts testing - fig S12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; load('motionartifact_EMG.mat')
%% fig S12B, left: plot 60 Hz notch filtered data w/ RMS envelope - STACKED EMG PLOT
figure() 
hold on 
plot(t_EMG,(datfiltEMG(1,:))./1e3,'b')
% plot(t_EMG,(RMSallEMG(1,:)./1e3),'r')
for i=2:size(amplifier_dataEMG,1)
    plot(t_EMG,(datfiltEMG(i,:)-700*(i-1))./1e3,'k')
%     plot(t_EMG,(RMSallEMG(i,:)-700*(i-1))./1e3,'r')
end
title('Pinch grip, raw data') 
ylabel('Signal (mV)','Fontsize',12)
xlabel('Time (s)','Fontsize',12)
ax = gca; ax.FontSize = 14; set(gcf,'position',[100,100,400,650])
xlim([20 25]); ylim([-4 0.4]) %limits for EMG contraction plot
%% plot load cell data corresponding to this pinch 
figure() 
plot(t_EMG,loadcell_EMG./1e3,'m')
xlabel('Time (s)','Fontsize',12)
ax = gca; ax.FontSize = 14; 
xlim([20 25]); 
%% fig S12B, right: plot BANDPASS filtered data w/ envelope - STACKED EMG PLOT
figure() 
hold on 
plot(t_EMG,(datfilt2EMG(1,:))./1e3,'b')
% plot(t_EMG,(RMSallEMG(1,:)./1e3),'r')
for i=2:size(amplifier_dataEMG,1)
    plot(t_EMG,(datfilt2EMG(i,:)-700*(i-1))./1e3,'k')
%     plot(t_EMG,(RMSallEMG(i,:)-700*(i-1))./1e3,'r')
end
% title('Pinch grip, filtered data') 
ylabel('Signal (mV)','Fontsize',12)
xlabel('Time (s)','Fontsize',12)
ax = gca; ax.FontSize = 14; set(gcf,'position',[100,100,400,650])
xlim([20 25]); ylim([-4 0.4]) %limits for EMG contraction plot

%% fig S12C, left: plot 60 Hz notch filtered data w/ envelope - STACKED MOTION ARTIFACT PLOT
figure() 
hold on 
plot(t_MA,(datfiltMA(1,:))./1e3,'b')
% plot(t_MA,(RMSallMA(1,:)./1e3),'r')
for i=[2,3,4,5,6] 
    plot(t_MA,(datfiltMA(i,:)-700*(i-1))./1e3,'k')
%     plot(t_MA,(RMSallMA(i,:)-700*(i-1))./1e3,'r')
end
title('Motion only, raw data') 
ylabel('Signal (mV)','Fontsize',12)
xlabel('Time (s)','Fontsize',12)
ax = gca; ax.FontSize = 14; set(gcf,'position',[100,100,400,650])
xlim([18.5 23.5]); ylim([-4 0.4]) %limits for MA plot

%% fig S12C, right: plot BANDPASS filtered data w/ envelope - STACKED MOTION ARTIFACT PLOT
figure() 
hold on 
plot(t_MA,(datfilt2MA(1,:))./1e3,'b')
% plot(t_MA,(RMSallMA(1,:)./1e3),'r')
for i=[2,3,4,5,6]  %note: there's a few bad channel in this recording that we'll discard
    plot(t_MA,(datfilt2MA(i,:)-700*(i-1))./1e3,'k')
%     plot(t_MA,(RMSallMA(i,:)-700*(i-1))./1e3,'r')
end
title('Motion only, filtered data') 
ylabel('Signal (mV)','Fontsize',12)
xlabel('Time (s)','Fontsize',12)
ax = gca; ax.FontSize = 14; set(gcf,'position',[100,100,400,650])
xlim([18.5 23.5]); ylim([-4 0.4]) %limits for MA plot

%% fig S2A: XRD analysis
clear all; load('XRD_of_MXtrodes.mat')
figure() 
plot(Theta,MXtrode,'k')
hold on 
plot(Theta,(Textile-(1.5e5))*0.5,'b')
xlim([2 60])
ylim([-1.2e5 3.7e5])
ylabel('Intensity (a.u.)')
xlabel('2\Theta (deg.)')
legend('MXtrode','pristine textile')
ax = gca; ax.FontSize = 14; set(gcf,'position',[100,100,400,300])

%% fig S2B: Raman spectroscopy
clear all; load('Raman_of_MXtrodes.mat') 
figure() 
plot(Lambda, smooth(MXtrode.*1.4-200,3), 'k') 
hold on 
plot(Lambda, smooth(Textile-1400,3), 'b') 
ylabel('Intensity (a.u.)')
xlabel('Raman Shift (cm^-1)')
% legend('MXtrode','pristine textile')
ax = gca; ax.FontSize = 14; set(gcf,'position',[100,100,400,300])
xlim([min(Lambda) max(Lambda)])
ylim([-900 1600])

%% fig S7C: CV cycling stability in MXene window
clear all; load('CVstability_MX500um_50cycles.mat')
% input=readtable('CV_cycling_stability_MX500um.xlsx','Sheet','MX 500 um', 'Range', 'A5:L4605');
% mat=table2array(input);
% V=mat(:,1);
% cycles=[2,10,20,30,40,50]; % choose a subset of the CV cycles to plot 
% I=mat(:,2:2:6);
d=0.05; %electrode diameter is 500 um = 0.05 cm [cm]
gsa=pi*((d/2)^2); %geometric surface area [cm^2]

figure() 
plot(V,(I(:,1)./gsa)*1000,'-','Color',[cmap(1,:)],'Linewidth',2.5)
hold on 
plot(V,(I(:,2)./gsa)*1000,'-','Color',[cmap(2,:)],'Linewidth',2.5)
plot(V,(I(:,3)./gsa)*1000,'-','Color',[cmap(3,:)],'Linewidth',2.5) 
plot(V,(I(:,4)./gsa)*1000,'-','Color',[cmap(4,:)],'Linewidth',2.5) 
plot(V,(I(:,5)./gsa)*1000,'-','Color',[cmap(5,:)],'Linewidth',2.5) 
plot(V,(I(:,6)./gsa)*1000,'-','Color',[cmap(6,:)],'Linewidth',2.5) 
ylabel('Current Density (mA/cm^2)','Fontsize',12); 
xlabel('Voltage (V)','Fontsize',12); 
% title('CV in PBS')
xlim([-1.8 0.7])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[100,400,450,350])
ylim([-75 100])

%% fig S15: Stimulation stability
% 500 um MXtrodes, EIS pre and post 1000 cycles of 1.2 mA stimulation 
clear all; load('MX_stimstability.mat') 
% input=readtable('Stim_stability.xlsx','Sheet','CVs', 'Range', 'A5:D4605');
% mat=table2array(input);
% V=mat(:,1);
% I_post=mat(:,2);
% I_pre=mat(:,4);

%|Z| plot 
figure()
plot(f,z_ch2_pre,'.','MarkerSize',16,'Color',cmap(1,:))
hold on 
plot(f,z_ch2_post,'.','MarkerSize',16,'Color',cmap(3,:))
legend({'Pre Stim','Post stim'},'Fontsize',12);legend boxoff 
ylabel('Impedance, (\Omega)','Fontsize',12); 
xlabel('Frequency (Hz)','Fontsize',12); 
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
ax = gca; ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[400,400,470,350])
ylim([1e2 1e4])

% PHASE plot 
figure()
plot(f,p_ch2_pre,'.','MarkerSize',16,'Color',cmap(1,:))
hold on 
plot(f,p_ch2_post,'.','MarkerSize',16,'Color',cmap(3,:))
ylabel('Phase, (degrees)','Fontsize',12); 
xlabel('Frequency (Hz)','Fontsize',12); 
ylim([-30 0])
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1; ax.TickLength=[0.02 0.025];
set(gcf,'position',[900,400,470,350])

% calculating |Z| % change at 1 kHz
abs((z_ch2_pre(21)-z_ch2_post(21))/z_ch2_pre(21))

% change in CV pre and post stim 
d=0.05; %electrode diameter is 500 um = 0.05 cm [cm]
gsa=pi*((d/2)^2); %geometric surface area [cm^2]
% figure() 
% plot(V,(I_pre./gsa)*1000,'-','Color',[cmap(1,:)],'Linewidth',2.5) 
% hold on 
% plot(V,(I_post./gsa)*1000,'-','Color',[cmap(3,:)],'Linewidth',2.5)
% ylabel('Current Density (mA/cm^2)','Fontsize',12); 
% xlabel('Voltage (V)','Fontsize',12); 
% legend('pre','post')
% % title('CV in PBS')
% xlim([-1.8 0.7])
% ylim([-75 100])
% ax = gca;
% ax.FontSize = 16; 
% box on
% ax.LineWidth=1.1; ax.TickLength=[0.02 0.025];
% set(gcf,'position',[900,400,470,350])

% change in CSC
% CSCpre=[];CSCpost=[];
% inds=find(I_pre(:,num)<0); %look only at the cathodal current 
% CSCpre(num)=abs(trapz(t(inds),((I_pre(inds,num)).*1000)))/gsa % units are [mC cm^-2]
% inds=find(I_post(:,num)<0); %look only at the cathodal current 
% CSCpost(num)=abs(trapz(t(inds),((I_post(inds,num)).*1000)))/gsa % units are [mC cm^-2]
% abs((CSCpre-CSCpost)/CSCpre)

%% fig S16C: MRI human scan - Neoptix temperature data plot
clear all; load('Neoptix_MRI_temp_data.mat')
figure()
plot(t(1:7741)./60,temp(1:7741),'k') 
xlabel('Time (min)') 
ylabel(['Temperature (' char(176) 'C)']) 
ylim([18.1 18.6])
xlim([0 t(7741)/60])
ax = gca;
ax.FontSize = 16; 
box on
ax.LineWidth=1.1;
ax.TickLength=[0.02 0.025];
set(gcf,'position',[100,400,800,250])
