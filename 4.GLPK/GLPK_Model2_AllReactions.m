% THIS CODE IS FOR MODEL 2

keep3 modelAll;   %% Clears all variables in the workspace except modelAll

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

xoo=[];
zoo=[];
too=[];
foo =[];

% Experimental Data Reader for 18 Time Frames
Exp_Data = importdata('/Users/debjit/Desktop/MyModel/SubNetwork/Model2/1.Protein Data Preparation/Experimental Data/5.Experimental Data.txt');   % Use Exp_Data.data for using the data part of it, ROW1 corresponds to First TF
TF_Exp = Exp_Data.data;


% Expression Data Reader for 18 Time Frames
Protein_Data = importdata('/Users/debjit/Desktop/MyModel/SubNetwork/Model2/1.Protein Data Preparation/8.62Finale_Expression.txt');   % Use Protein_Data.data for using the data part of it, COLUMN1 corresponds to First TF
TF_Protein = Protein_Data.data;



%$$$$$$$$$$$$$$$$$$       LOOP    For the 10 differenet objective functions    $$$$$$$$$$$$$$$$%
for k=1:10
    
    this_Model=modelAll(k,1);   % Each Model is chosen in a serial order
    
    for m=1:18
        
        finalModel=changeRxnBounds(this_Model,this_Model.rxns,TF_Protein(:,m),'u');          % Upper bounds in the model are fixed to the expression values for that TF
        finalModel=changeRxnBounds(finalModel,'b6_b7',0,'u');   
        temp = finalModel;
        Deficient=this_Model.rxns;
        
        for def=1:(length(Deficient)+1)
            
            if (def <=length(Deficient))
                finalModel=changeRxnBounds(temp,Deficient(def),0,'u'); 
            else
                finalModel=temp;
            end;
                
  %------------------------------------------------------------------------------------------------------------------------------------------------------%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             MAXIMIZATION        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %------------------------------------------------------------------------------------------------------------------------------------------------------%
  
            % OPTIMIZING THE MODEL'
             FBAmax=optimizeCbModel(finalModel,'max');
             Z=FBAmax.f
             zValues(m,def)=Z;
             
             LOPT=FBAmax.x;
             
             solutionmodelAll(k) =FBAmax;
             
             [MIN_FVA,MAX_FVA]=fluxVariability(finalModel,100,'max');
             
% Pathway Average

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gluconeogenesis = v7 --> v18
%   Glycogenolysis = v19 --> v30 +b10( v52)+ v1m(v34)
%   TCA = v35(v2m) --> v42(v9m) 
%   Glyoxylate = v4 + v5 
%   Glutamate = v32+v33+v51(b9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
SGluconeogenesisLOPT= [LOPT(7),LOPT(8),LOPT(9),LOPT(10),LOPT(11),LOPT(12),LOPT(13),LOPT(14),LOPT(15),LOPT(16),LOPT(17),LOPT(18)];

SGlyocgenolysisLOPT= [LOPT(19),LOPT(20),LOPT(21),LOPT(22),LOPT(23),LOPT(24),LOPT(25),LOPT(26),LOPT(27),LOPT(28),LOPT(29),LOPT(30),LOPT(52),LOPT(34)];

STCALOPT= [LOPT(35),LOPT(36),LOPT(37),LOPT(38),LOPT(39),LOPT(40),LOPT(41),LOPT(42)];

SGlyoxylateLOPT=[LOPT(4),LOPT(5)];

SGlutamateLOPT=[LOPT(32),LOPT(33),LOPT(51)];


            TLOPT=[];
            TLOPT(1)=mean(SGluconeogenesisLOPT);

            TLOPT(2)=mean(SGlyocgenolysisLOPT);

            TLOPT(3)=mean(STCALOPT);

            TLOPT(4)=mean(SGlyoxylateLOPT);

            TLOPT(5)=mean(SGlutamateLOPT);
            
         
% Product Average
 
SNucleotideLOPT = [LOPT(53)];
    
SAminoLOPT =[LOPT(54),LOPT(55),LOPT(56),LOPT(57),LOPT(58), LOPT(33)];

SLipidLOPT =[LOPT(59),LOPT(60),LOPT(61),LOPT(62)];
        

            TLOPT(6)= mean(SNucleotideLOPT);
    
            TLOPT(7) = mean(SAminoLOPT);

            TLOPT(8)= mean(SLipidLOPT);


AverageLOPT=[AverageLOPT;TLOPT];  

Individual_LOPT=[Individual_LOPT;SGluconeogenesisLOPT, SGlyocgenolysisLOPT, STCALOPT, SGlyoxylateLOPT, SGlutamateLOPT, SNucleotideLOPT, SAminoLOPT, SLipidLOPT]; 



%************************************* FLUX VARIABILITY ALONG THE 8 DIFFERENT PATHWAYS (Individual and average)***************************************

% ---------------------------- MIN -----------------------------------------%

% FLUX VARIABILITY (MIN) along each of the 5 pathways%
 
SGluconeogenesisMIN_FVA= [MIN_FVA(7),MIN_FVA(8),MIN_FVA(9),MIN_FVA(10),MIN_FVA(11),MIN_FVA(12),MIN_FVA(13),MIN_FVA(14),MIN_FVA(15),MIN_FVA(16),MIN_FVA(17),MIN_FVA(18)];

SGlyocgenolysisMIN_FVA= [MIN_FVA(19),MIN_FVA(20),MIN_FVA(21),MIN_FVA(22),MIN_FVA(23),MIN_FVA(24),MIN_FVA(25),MIN_FVA(26),MIN_FVA(27),MIN_FVA(28),MIN_FVA(29),MIN_FVA(30),MIN_FVA(52),MIN_FVA(34)];

STCAMIN_FVA= [MIN_FVA(35),MIN_FVA(36),MIN_FVA(37),MIN_FVA(38),MIN_FVA(39),MIN_FVA(40),MIN_FVA(41),MIN_FVA(42)];

SGlyoxylateMIN_FVA=[MIN_FVA(4),MIN_FVA(5)];

SGlutamateMIN_FVA=[MIN_FVA(32),MIN_FVA(33),MIN_FVA(51)];

            TMIN_FVA=[];
            
            TMIN_FVA(1)=mean(SGluconeogenesisMIN_FVA);

            TMIN_FVA(2)=mean(SGlyocgenolysisMIN_FVA);

            TMIN_FVA(3)=mean(STCAMIN_FVA);

            TMIN_FVA(4)=mean(SGlyoxylateMIN_FVA);

            TMIN_FVA(5)=mean(SGlutamateMIN_FVA);
            

% FLUX VARIABILITY  (MIN) for the 3 products
 
SNucleotideMIN_FVA = [MIN_FVA(53)];
    
SAminoMIN_FVA  =[MIN_FVA(54),MIN_FVA(55),MIN_FVA(56),MIN_FVA(57),MIN_FVA(58), MIN_FVA(33)];

SLipidMIN_FVA  =[MIN_FVA(59),MIN_FVA(60),MIN_FVA(61),MIN_FVA(62)];
   
        TMIN_FVA(6) = mean(SNucleotideMIN_FVA );
    
        TMIN_FVA(7)  = mean(SAminoMIN_FVA );

        TMIN_FVA(8) = mean(SLipidMIN_FVA );


 AverageMIN_8Pathway=[AverageMIN_8Pathway;TMIN_FVA]; 

MIN_8Pathway=[MIN_8Pathway;SGluconeogenesisMIN_FVA, SGlyocgenolysisMIN_FVA, STCAMIN_FVA, SGlyoxylateMIN_FVA, SGlutamateMIN_FVA, SNucleotideMIN_FVA, SAminoMIN_FVA, SLipidMIN_FVA];



% ---------------------------- MAX -----------------------------------------%

% FLUX VARIABILITY (MAX) along each of the 5 pathways%
 
SGluconeogenesisMAX_FVA= [MAX_FVA(7),MAX_FVA(8),MAX_FVA(9),MAX_FVA(10),MAX_FVA(11),MAX_FVA(12),MAX_FVA(13),MAX_FVA(14),MAX_FVA(15),MAX_FVA(16),MAX_FVA(17),MAX_FVA(18)];

SGlyocgenolysisMAX_FVA= [MAX_FVA(19),MAX_FVA(20),MAX_FVA(21),MAX_FVA(22),MAX_FVA(23),MAX_FVA(24),MAX_FVA(25),MAX_FVA(26),MAX_FVA(27),MAX_FVA(28),MAX_FVA(29),MAX_FVA(30),MAX_FVA(52),MAX_FVA(34)];

STCAMAX_FVA= [MAX_FVA(35),MAX_FVA(36),MAX_FVA(37),MAX_FVA(38),MAX_FVA(39),MAX_FVA(40),MAX_FVA(41),MAX_FVA(42)];

SGlyoxylateMAX_FVA=[MAX_FVA(4),MAX_FVA(5)];

SGlutamateMAX_FVA=[MAX_FVA(32),MAX_FVA(33),MAX_FVA(51)];

            TMAX_FVA=[];
            
            TMAX_FVA(1)=mean(SGluconeogenesisMAX_FVA);

            TMAX_FVA(2)=mean(SGlyocgenolysisMAX_FVA);

            TMAX_FVA(3)=mean(STCAMAX_FVA);

            TMAX_FVA(4)=mean(SGlyoxylateMAX_FVA);

            TMAX_FVA(5)=mean(SGlutamateMAX_FVA);
            

% FLUX VARIABILITY  (MAX) for the 3 products
 
SNucleotideMAX_FVA = [MAX_FVA(53)];
    
SAminoMAX_FVA  =[MAX_FVA(54),MAX_FVA(55),MAX_FVA(56),MAX_FVA(57),MAX_FVA(58), MAX_FVA(33)];

SLipidMAX_FVA  =[MAX_FVA(59),MAX_FVA(60),MAX_FVA(61),MAX_FVA(62)];
   
        TMAX_FVA(6) = mean(SNucleotideMAX_FVA );
    
        TMAX_FVA(7)  = mean(SAminoMAX_FVA );

        TMAX_FVA(8) = mean(SLipidMAX_FVA );


 AverageMAX_8Pathway=[AverageMAX_8Pathway;TMAX_FVA]; 

MAX_8Pathway=[MAX_8Pathway;SGluconeogenesisMAX_FVA, SGlyocgenolysisMAX_FVA, STCAMAX_FVA, SGlyoxylateMAX_FVA, SGlutamateMAX_FVA, SNucleotideMAX_FVA, SAminoMAX_FVA, SLipidMAX_FVA];


% All Optimal solutions
            FluxLOPT=[FluxLOPT;LOPT'];
            
 % All Flux Variability
            FVA_MIN_All=[FVA_MIN_All;MIN_FVA'];
            FVA_MAX_All=[FVA_MAX_All;MAX_FVA'];
            SPAN=MAX_FVA-MIN_FVA;
            SPAN_ALL=[SPAN_ALL;SPAN'];
            
  
 % COMPARISON STUDIES
 
 Rexp1=TF_Exp(m,:);
            
% ******************  EUCLIDEAN DISTANCE CALCULATIONS  ******************

% FOR Time Frame 1 
        Exp_Pathway1=[];
        Comp_Pathway1=[];
        Matrix1=[];
        
        for d=1:8
                TSPAN_FVA(d)=(TMAX_FVA(d)-TMIN_FVA(d))/1000;
                if (TSPAN_FVA(d)<=span_Avg) TSPAN_FVA(d)=span_Avg; end;
                Exp_Pathway1(d)=Rexp1(d)*TSPAN_FVA(d);
                Comp_Pathway1(d)=(TLOPT(d)/1000)*TSPAN_FVA(d);
        end;
        Matrix1(1,:)= Exp_Pathway1;
        Matrix1(2,:)= Comp_Pathway1;
        
        Distance_Time1= pdist(Matrix1,'euclidean');
        Euclidean(m,def)= Distance_Time1;
        
                        
% ******************  PEARSON AND SPEARMAN CORRELATION CALCULATIONS  ******************

   P_Time1=0;
        P_Time1=corr(Rexp1',(TLOPT/1000)');
        Pearson(m,def)=P_Time1;
        
   S_Time1=0;
        S_Time1=corr(Rexp1',(TLOPT/1000)','type','Spearman');
        Spearman(m,def)=S_Time1;
        
        end;
        
    end;
    xoo=[xoo;Euclidean];
    zoo=[zoo;Pearson];
    too=[too;Spearman];
    foo=[foo;zValues];
end;


dlmwrite('Results_M2/1A.Average_Pathway_Optimal.csv',AverageLOPT);
dlmwrite('Results_M2/1B.Average_Pathway_MIN.csv', AverageMIN_8Pathway);
dlmwrite('Results_M2/1C.Average_Pathway_MAX.csv', AverageMAX_8Pathway);

dlmwrite('Results_M2/2A.Individual_8Pathway_Optimal.csv',Individual_LOPT);
dlmwrite('Results_M2/2B.Individual_8Pathway_MIN.csv',MIN_8Pathway);
dlmwrite('Results_M2/2C.Individual_8Pathway_MAX.csv',MAX_8Pathway);

dlmwrite('Results_M2/3A.Individual_All_Optimal.csv',FluxLOPT);
dlmwrite('Results_M2/3B.Individual_All_MIN.csv',FVA_MIN_All);
dlmwrite('Results_M2/3C.Individual_All_MAX.csv',FVA_MAX_All);

dlmwrite('Results_M2/Euclidean.csv',xoo);

dlmwrite('Results_M2/Pearson.csv',zoo);

dlmwrite('Results_M2/Spearman.csv',too);

dlmwrite('Results_M2/zValues.csv',foo);



%h = waitbar(0,'Initializing ...');
%s2 = num2str(k);
%s1= 'Objective     ';
%Str = strcat(s1, s2);
%waitbar(k/10,h,Str);








