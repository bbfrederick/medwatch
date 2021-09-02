clear all;
close all;
clc;

% a=importdata('result0.txt');

fid = fopen('result0.txt');  %open file in binary mode
bytes = fread(fid, [1 Inf], 'uint8');        %read the whole lot as bytes
fclose(fid);

spo2=[];
for i=1:22:length(bytes),  
    spo2=[spo2 bytes(i+11)];  
end
plot(spo2);