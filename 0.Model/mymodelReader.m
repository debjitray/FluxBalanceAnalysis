% THIS CODE FIRST READS THE MODEL (.XML) FILE
% THEN IT CREATES A 10 NEW MODELS FOR THE 10 DIFFERENT OBJECTIVE FUNCTIONS
% THESE MODELS ARE STORED SERIALLY IN THE 'modelAll'

function modelAll=mymodelReader(filename)

model=readCbModel(filename);
clc;

disp '<><><><><><><><><><><><><><><><><><><><><><><><><><>' 
disp '-----------------------------------------------------------------' 
disp 'net ATP model initialized and saved as --> 'modelAll(1,1)'' 
rxnListnetATP= {'V6' 'V25' 'V32' 'V4m' 'V5m' 'V9m' 'Nucleotide1' 'Amino1' 'V7m' 'V26' 'V29' 'V6m' 'V1' 'V7' 'V10' 'V22' 'V30' 'Lipid1' 'Lipid2' 'V11' 'V33' 'Lipid3'};
coEffnetATP=   [3	3	3	3	3	3	3	3	2	1	1	1	-1	-1	-1	-1	-1	-1	-1	-3	-3	-3];
modelnetATP=changeObjective(model,rxnListnetATP,coEffnetATP);
disp '-----------------------------------------------------------------' 


disp 'ATP production model initialized and saved as --> 'modelAll(2,1)''
rxnListATPp= {'V6' 'V25' 'V32' 'V4m' 'V5m' 'V9m' 'Nucleotide1' 'Amino1' 'V7m' 'V26' 'V29' 'V6m'};
coEffATPp=   [3		3	3	3	3	3	3	3	2	1	1	1];
modelATPp=changeObjective(model,rxnListATPp,coEffATPp);
disp '-----------------------------------------------------------------' 


disp 'ATP consumption model initialized and saved as --> 'modelAll(3,1)''
rxnListATPc= {'V11' 'V33' 'Lipid3' 'V1' 'V7' 'V10' 'V22' 'V30' 'Lipid1' 'Lipid2'};
coEffATPc=   [3		3	3	1	1	1	1	1	1	1];
modelATPc=changeObjective(model,rxnListATPc,coEffATPc);
disp '-----------------------------------------------------------------' 


disp 'Nucleotide model initialized and saved as --> 'modelAll(4,1)''
rxnListNucleotide= {'Nucleotide1'};
coEffNucleotide=   [1];
modelNucleotide=changeObjective(model,rxnListNucleotide,coEffNucleotide);
disp '-----------------------------------------------------------------' 


disp 'Amino model initialized and saved as --> 'modelAll(5,1)''
rxnListAmino= {'Amino1' 'Amino2' 'Amino3' 'Amino4' 'Amino5'};
coEffAmino=   [1	1	1	1	1];
modelAmino=changeObjective(model,rxnListAmino,coEffAmino);
disp '-----------------------------------------------------------------' 


disp 'Lipid model initialized and saved as --> 'modelAll(6,1)''
rxnListLipid= {'Lipid1' 'Lipid2' 'Lipid3' 'Lipid4'};
coEffLipid=   [1	1	1	1];
modelLipid=changeObjective(model,rxnListLipid,coEffLipid);
disp '-----------------------------------------------------------------' 


disp 'Carbohydrate Production model initialized and saved as --> 'modelAll(7,1)''
rxnListCarboProduction= {'V17' 'V18'};
coEffCarboProduction=   [1	1];
modelCarboProduction=changeObjective(model,rxnListCarboProduction,coEffCarboProduction);
disp '-----------------------------------------------------------------' 



disp 'Carbohydrate Uptake model initialized and saved as --> 'modelAll(8,1)''
rxnListCarboUptake= {'V19'};
coEffCarboUptake=   [1];
modelCarboUptake=changeObjective(model,rxnListCarboUptake,coEffCarboUptake);
disp '-----------------------------------------------------------------' 


disp 'Acetate Uptake model initialized and saved as --> 'modelAll(9,1)''
rxnListAcetate= {'bext'};
coEffAcetate=   [1];
modelAcetate=changeObjective(model,rxnListAcetate,coEffAcetate);
disp '-----------------------------------------------------------------' 


disp 'Glutamate Uptake model initialized and saved as --> 'modelAll(10,1)''
rxnListGlutamate= {'V33'};
coEffGlutamate=   [1];
modelGlutamate=changeObjective(model,rxnListGlutamate,coEffGlutamate);
disp '-----------------------------------------------------------------' 




modelAll(1,1)=modelnetATP;
modelAll(2,1)=modelATPp;
modelAll(3,1)=modelATPc;
modelAll(4,1)=modelNucleotide;
modelAll(5,1)=modelAmino;
modelAll(6,1)=modelLipid;
modelAll(7,1)=modelCarboProduction;
modelAll(8,1)=modelCarboUptake;
modelAll(9,1)=modelAcetate;
modelAll(10,1)=modelGlutamate;

disp '*********************************************************************************';
disp 'All 10 models saved in modelAll in their order';
disp '*********************************************************************************';

end











