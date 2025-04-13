%%% This code uses 'Batch_smFISH_v2_10.m' file as a function to run it
%%% repeatedly with varying threshold values set as range for exon and
%%% intron channels.

%%% indicate the type and order of the RNA channel to be analyzed.
%%% Put '1' for nucleus, '2' for exon, '3' for intron, '4' for protein IF, '5' for DTC channel.
%%% and '0' for skipping that channel.
%%% For example, if you have three channels in the order of exon Ch1, intron,
%%% and exon Ch2 then use: ChOrder = [ 2 1 2 ].

inputPath = 'E:\Mutations\N2 Control'; % <-- change folder as desired.
masterOutputPath = 'E:\Mutations\N2 Control\04012025 Results'; %%% Specify the master output folder location.
ExonThresRange = 2:0.5:5 ;            %%% Set range for exon threshold.
IntronThresRange = 1.5:0.01:1.5 ; %%% Set range for intron threshold.
ChOrder = [  3 2 1  ]; %%% Change according to images used

cnE = length(ExonThresRange);
cnI = length(IntronThresRange);

for i = 1:cnE
    thresForExon = ExonThresRange(i);
    for j = 1:cnI
        thresForIntron = IntronThresRange(j);
        subPath = strcat('Ex_',num2str(thresForExon),'_Int_',num2str(thresForIntron));
        Batch_smFISH_v2_10(inputPath, masterOutputPath, subPath, thresForExon, thresForIntron, ChOrder);
    end
end
