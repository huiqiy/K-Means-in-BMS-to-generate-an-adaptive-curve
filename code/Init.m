clc; clear;clear all;
%% Initialisierung der Parameter für das Zellmodell der Parallelschaltung
Data=load('LGHE_Projekt_Parameter.mat');
data_FUDS=load('FUDS_LGHE4-2.mat');
data_DST=load('DST_LGHE4-2.mat');
%% An dieser Stelle sollen die Parameter neu berechnet werden (R0,R1,Tau1)
%n_par = 4;
C_cell = Data.LGHE4.QOCV.C;
CN = C_cell;
eta2 = 0.985; 
tsample=1000;
R0=[Data.LGHE4.Parameter.R0];
R1=[Data.LGHE4.Parameter.R1];
R2=[Data.LGHE4.Parameter.R2];
Tau1=[Data.LGHE4.Parameter.Tau1];
Tau2=[Data.LGHE4.Parameter.Tau2];

SOCData=Data.LGHE4.QOCV.Ah/Data.LGHE4.QOCV.C;
SOC_Param=Data.LGHE4.Parameter.SOC;
OCVData=Data.LGHE4.QOCV.Voltage;

% Current=[data_FUDS.data.time,data_FUDS.data.current];
% Voltage=[data_FUDS.data.time,data_FUDS.data.voltage];
Temperature=[data_DST.data.time,data_DST.data.temperature];
%Current=[data_DST.data.time,data_DST.data.current];
%Voltage=[data_DST.data.time,data_DST.data.voltage];
%Temperature=[data_DST.data.time,data_DST.data.temperature];
Current=[data_DST.data.time(3:end)-data_DST.data.time(3),data_DST.data.current(3:end)];
%Voltage=[data_FUDS.data.time(2:end)-data_FUDS.data.time(2),data_FUDS.data.voltage(2:end)];


%%
% plot(data_DST.data.time(2:end)-data_DST.data.time(2),data_DST.data.current(2:end));

