
Exp_Data = importdata('/Users/debjit/Documents/Ymet/MyModel/SubNetwork/Model2/1.Protein Data Preparation/Experimental Data/5.Experimental Data.txt');   % Use Protein_Data.data for using the data part of it, COLUMN1 corresponds to First TF
Exp = Exp_Data.data;


Gluco = importdata('A.Gluconeogenesis.csv');
Glyco = importdata('B.Glyocgenolysis.csv');
TCA = importdata('C.TCA.csv');
Gly = importdata('D.Glyoxylate.csv');
Gluta = importdata('E.Glutamate.csv');
N = importdata('F.Nucleotide.csv');
A = importdata('G.Amino.csv');
L = importdata('H.Lipid.csv');

i=1;
array=[];
test=[];


for m=1:18
    
    for i=1:5000

        LOPT = [Gluco(i,m),Glyco(i,m),TCA(i,m),Gly(i,m),Gluta(i,m),N(i,m),A(i,m),L(i,m)];
        
        test=[test;LOPT];

        array(m,i)=corr(Exp(m,:)',(LOPT/1000)');

    end;
    
end;

final=array';

dlmwrite('Pearson_5000.csv',final);

Mean_Values=[];
Max_Values=[];
Min_Values=[];
Std_Values=[];





for i=1:18
    Mean_Values=[Mean_Values;mean(final(:,i))];
    Max_Values=[Max_Values;max(final(:,i))];
    Min_Values=[Min_Values;min(final(:,i))];
    Std_Values=[Std_Values;std(final(:,i))];
         
 end;

All=[Mean_Values';Max_Values';Min_Values';Std_Values'];

dlmwrite('Analysis_5000.csv',All);




