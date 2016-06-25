% THIS CODE IS FOR MODEL 2 Genome Wide

% ************* EMPTY MATRIX FOR DATA STORAGE *************
AverageLOPT=[];
Individual_LOPT=[];

AverageMIN_8Pathway=[];
MIN_8Pathway=[];

AverageMAX_8Pathway=[];
MAX_8Pathway=[];

FluxLOPT=[];
FVA_MIN_All=[];
FVA_MAX_All=[];
SPAN_ALL=[];

span_Avg=0.05;   
array=[];
zoo=[];
too=[];
foo =[];
all=[];
Pearson=[];
Spearman=[];
doo=[];
% Experimental Data Reader for 18 Time Frames
Exp_Data = importdata('/Users/debjit/Documents/Ymet/MyModel/SubNetwork/Model2/1.Protein Data Preparation/Experimental Data/5.Experimental Data.txt');   % Use Exp_Data.data for using the data part of it, ROW1 corresponds to First TF
TF_Exp = Exp_Data.data;



%$$$$$$$$$$$$$$$$$$       LOOP    For the first 2 differenet objective functions    $$$$$$$$$$$$$$$$%
for k=1  
    
	this_Model=modelAllGW(k,1);   % Each Model is chosen in a serial order
    
	new=[];
    
    
	for m=1:6 % For 18 Time Frames

        
        my =[];
        h = waitbar(0,'Initializing ...');
        s2 = num2str(m);
        s1= 'Time     ';
        Str = strcat(s1, s2);
        
        for def=1:1   % Loop for deleting 905 + 1 genes one by one    
            if (def <=1)
                fileName = ['Finale_ExpressionGW_' num2str(def)];
                data.(fileName)=load([fileName '.txt']);
                TF_Protein=data.(fileName); 
                fileName
            else
                Protein_Data = importdata('/Users/debjit/Documents/Ymet/MyModel/GenomeWide/1.Protein Data Preparation/9.Expression_GW_Finale.txt');   % Use Protein_Data.data for using the data part of it, COLUMN1 corresponds to First TF
                TF_Protein = Protein_Data.data;
            end;
            finalModel=changeRxnBounds(this_Model,this_Model.rxns,TF_Protein(:,m),'u');          % Upper bounds in the model are fixed to the expression values for that TF
            finalModel=changeRxnBounds(finalModel,'SUCFUMtm',0,'u');  %% PUTTING b6_b7 as zero
        
            for g=1:length(finalModel.lb)  % Setting the lower bounds to the gene expression data
                if((finalModel.lb(g,1) <0))
                    finalModel.lb(g,1) = -(TF_Protein(g,m));
                end;
            end;
        
            finalModel=changeRxnBounds(finalModel,'EX_glc(e)',0,'l');            %% PUTTING THE GLUCOSE UPTAKE TO 0
            finalModel=changeRxnBounds(finalModel,'EX_nh4(e)',0,'l');           %% PUTTING THE nh4 UPTAKE TO 0
            finalModel=changeRxnBounds(finalModel,'EX_ac(e)',-1000,'l');     %% PUTTING THE ACETATE UPTAKE TO -1000 
            finalModel=changeRxnBounds(finalModel,'EX_o2(e)',-1000,'l');     %% PUTTING THE O2 UPTAKE TO -1000 

            waitbar(def/906,h,Str);
  %------------------------------------------------------------------------------------------------------------------------------------------------------%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             MAXIMIZATION        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %------------------------------------------------------------------------------------------------------------------------------------------------------%
  
            % OPTIMIZING THE MODEL'
            FBAmax=optimizeCbModel(finalModel,'max');
            didi=FBAmax.stat
            
            my=[my,didi];
            
            Z=FBAmax.f;
            zValues(m,def)=Z;
            LOPT=FBAmax.x;
                    
% Pathway Average

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gluconeogenesis = v7 --> v18
%   Glycogenolysis = v19 --> v30 +b10( v52)+ v1m(v34)
%   TCA = v35(v2m) --> v42(v9m) 
%   Glyoxylate = v4 + v5 
%   Glutamate = v32+v33+v51(b9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (numel(LOPT)==0)  % Taking care of the no solution problem when the LOPT is empty numel(LOPT) gives the size of the matrix
                array=[k,m,def];
                new=[new;array];

                LOPT=zeros(1577,1);
            end;
    
            SGluconeogenesisLOPT= ([LOPT(1307),-GtZero(LOPT(463)),ltZero(LOPT(1257)),ltZero(LOPT(1255)),-GtZero(LOPT(800)),-GtZero(LOPT(1503)),-GtZero(LOPT(734)),LOPT(737),-GtZero(LOPT(1254)),-GtZero(LOPT(1258)),LOPT(1509),LOPT(827)]);

            SGlyocgenolysisLOPT= ([LOPT(826),ltZero(LOPT(1258)),ltZero(LOPT(1254)),LOPT(1248),ltZero(LOPT(734)),ltZero(LOPT(1503)),ltZero(LOPT(800)),-GtZero(LOPT(1255)),-GtZero(LOPT(1257)),ltZero(LOPT(463)),LOPT(1373),LOPT(1230),ltZero(LOPT(1380)),LOPT(1238)]);

            STCALOPT= ([LOPT(354),ltZero(LOPT(110)),LOPT(961),LOPT(1237),-GtZero(LOPT(1459)),ltZero(LOPT(1455)),ltZero(LOPT(771)),ltZero(LOPT(1060))]);

            SGlyoxylateLOPT=([LOPT(965),LOPT(1044)]);

            SGlutamateLOPT=([LOPT(962),LOPT(841),ltZero(LOPT(1217))]);


            TLOPT=[];
            TLOPT(1)=mean(SGluconeogenesisLOPT);

            TLOPT(2)=mean(SGlyocgenolysisLOPT);

            TLOPT(3)=mean(STCALOPT);

            TLOPT(4)=mean(SGlyoxylateLOPT);

            TLOPT(5)=mean(SGlutamateLOPT);
            
         
% Product Average
 
            SNucleotideLOPT = [LOPT(788)];
    
            SAminoLOPT =[LOPT(1253),-GtZero(LOPT(244)),LOPT(97),-GtZero(LOPT(245)),LOPT(1331), LOPT(841)];

            SLipidLOPT =[ltZero(LOPT(85)),ltZero(LOPT(86)),LOPT(778),LOPT(1088)];
        

            TLOPT(6)= mean(SNucleotideLOPT);
    
            TLOPT(7) = mean(SAminoLOPT);

            TLOPT(8)= mean(SLipidLOPT);


            AverageLOPT=[AverageLOPT;TLOPT];  

            Individual_LOPT=[Individual_LOPT;SGluconeogenesisLOPT, SGlyocgenolysisLOPT, STCALOPT, SGlyoxylateLOPT, SGlutamateLOPT, SNucleotideLOPT, SAminoLOPT, SLipidLOPT]; 





% All Optimal solutions
            FluxLOPT=[FluxLOPT;LOPT'];
            
            
  
 % COMPARISON STUDIES
 
             Rexp1=TF_Exp(m,:);
            

% ******************  PEARSON AND SPEARMAN CORRELATION CALCULATIONS  ******************

               P_Time1=0;
                    P_Time1=corr(Rexp1',(TLOPT/1000)');
                    Pearson(m,def)=P_Time1;
        
               S_Time1=0;
                    S_Time1=corr(Rexp1',(TLOPT/1000)','type','Spearman');
                    Spearman(m,def)=S_Time1;
            end;
            close(h);
            foo=[foo;zValues];
            doo=[doo;my];
    end;

    zoo=[zoo;Pearson];
    too=[too;Spearman];
    all=[all;new];
end;


dlmwrite('Results_M2/1A.Average_Pathway_Optimal.csv',AverageLOPT);

dlmwrite('Results_M2/2A.Individual_8Pathway_Optimal.csv',Individual_LOPT);

dlmwrite('Results_M2/3A.Individual_All_Optimal.csv',FluxLOPT);

dlmwrite('Results_M2/Pearson.csv',zoo);

dlmwrite('Results_M2/Spearman.csv',too);

dlmwrite('Results_M2/zValues.csv',foo);

dlmwrite('Results_M2/ERROR.csv',all);

dlmwrite('Results_M2/Status.csv',doo);

%h = waitbar(0,'Initializing ...');
%s2 = num2str(k);
%s1= 'Objective     ';
%Str = strcat(s1, s2);
%waitbar(k/10,h,Str);








