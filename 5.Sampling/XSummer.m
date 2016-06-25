% SUBNETWORK MODEL

% ************* EMPTY MATRIX FOR DATA STORAGE *************


array_Gluconeogenesis=[];
array_Glyocgenolysis=[];
array_TCA=[];
array_Glyoxylate=[];
array_Glutamate=[];

array_Nucleotide=[];
array_Amino=[];
array_Lipid=[];

for i=1:18;

    TIME = ['Time' num2str(i)];

    a=load([TIME '/DR_Yeast_Flux4_1.mat']);
    b=load([TIME '/DR_Yeast_Flux4_2.mat']);
    c=load([TIME '/DR_Yeast_Flux4_3.mat']);
    d=load([TIME '/DR_Yeast_Flux4_4.mat']);
    e=load([TIME '/DR_Yeast_Flux4_5.mat']);
    
    
    GOOL=[a.points,b.points,c.points,d.points,e.points];

    
    if (i~=14)
    

        SGluconeogenesis= (GOOL(7,:)+ GOOL(8,:)+GOOL(9,:)+GOOL(10,:)+GOOL(11,:)+GOOL(12,:)+GOOL(13,:)+GOOL(14,:)+GOOL(15,:)+GOOL(16,:)+GOOL(17,:)+GOOL(18,:))/12;
    
        SGlyocgenolysis= (GOOL(19,:)+GOOL(20,:)+GOOL(21,:)+GOOL(22,:)+GOOL(23,:)+GOOL(24,:)+GOOL(25,:)+GOOL(26,:)+GOOL(27,:)+GOOL(28,:)+GOOL(29,:)+GOOL(30,:)+GOOL(50,:)+GOOL(33,:))/14;
    
        STCA= (GOOL(34,:)+GOOL(35,:)+GOOL(36,:)+GOOL(37,:)+GOOL(38,:)+GOOL(39,:)+GOOL(40,:)+GOOL(41,:))/8;

        SGlyoxylate=(GOOL(4,:)+GOOL(5,:))/2;

        SGlutamate=(GOOL(31,:)+GOOL(32,:)+GOOL(49,:))/3;


 
        SNucleotide = (GOOL(51,:))/1;
    
        SAmino =(GOOL(52,:)+GOOL(53,:)+GOOL(54,:)+GOOL(55,:)+GOOL(56,:)+ GOOL(32,:))/6;

        SLipid =(GOOL(57,:)+GOOL(58,:)+GOOL(59,:)+GOOL(60,:))/4;
    
    
    
    else
        
        SGluconeogenesis= (GOOL(7,:)+ GOOL(8,:)+GOOL(9,:)+GOOL(10,:)+GOOL(11,:)+GOOL(12,:)+GOOL(13,:)+GOOL(14,:)+GOOL(15,:)+GOOL(16,:)+GOOL(17,:)+GOOL(18,:))/12;
    
        SGlyocgenolysis= (GOOL(19,:)+GOOL(20,:)+GOOL(21,:)+GOOL(22,:)+GOOL(23,:)+GOOL(24,:)+GOOL(25,:)+GOOL(26,:)+GOOL(27,:)+GOOL(28,:)+GOOL(29,:)+GOOL(30,:)+GOOL(50,:)+GOOL(33,:))/14;
    
        STCA= (GOOL(34,:)+GOOL(35,:)+GOOL(36,:)+GOOL(37,:)+GOOL(38,:)+GOOL(39,:)+GOOL(40,:)+GOOL(41,:))/8;

        SGlyoxylate=(GOOL(4,:)+GOOL(5,:))/2;

        SGlutamate=(GOOL(31,:)+GOOL(32,:)+GOOL(49,:))/3;


 
        SNucleotide = (GOOL(51,:))/1;
    
        SAmino =(GOOL(52,:)+GOOL(53,:)+GOOL(54,:)+GOOL(55,:)+GOOL(56,:)+ GOOL(32,:))/6;

        SLipid =(GOOL(57,:)+GOOL(58,:)+GOOL(59,:))/3;
    
    end;
    
    
    
    array_Gluconeogenesis=[array_Gluconeogenesis;SGluconeogenesis];
    array_Glyocgenolysis=[array_Glyocgenolysis;SGlyocgenolysis];
    array_TCA=[array_TCA;STCA];
    array_Glyoxylate=[array_Glyoxylate;SGlyoxylate];
    array_Glutamate=[array_Glutamate;SGlutamate];
    
    array_Nucleotide=[array_Nucleotide;SNucleotide];
    array_Amino=[array_Amino;SAmino];
    array_Lipid=[array_Lipid;SLipid];
    
end;

    
    
dlmwrite('/Users/debjit/Desktop/sampleCbModel/Test1/5000/A.Gluconeogenesis.csv',array_Gluconeogenesis');

dlmwrite('/Users/debjit/Desktop/sampleCbModel/Test1/5000/B.Glyocgenolysis.csv',array_Glyocgenolysis');

dlmwrite('/Users/debjit/Desktop/sampleCbModel/Test1/5000/C.TCA.csv',array_TCA');

dlmwrite('/Users/debjit/Desktop/sampleCbModel/Test1/5000/D.Glyoxylate.csv',array_Glyoxylate');

dlmwrite('/Users/debjit/Desktop/sampleCbModel/Test1/5000/E.Glutamate.csv',array_Glutamate');



dlmwrite('/Users/debjit/Desktop/sampleCbModel/Test1/5000/F.Nucleotide.csv',array_Nucleotide');

dlmwrite('/Users/debjit/Desktop/sampleCbModel/Test1/5000/G.Amino.csv',array_Amino');

dlmwrite('/Users/debjit/Desktop/sampleCbModel/Test1/5000/H.Lipid.csv',array_Lipid');



