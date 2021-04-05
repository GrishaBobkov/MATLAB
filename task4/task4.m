xyz=[0 0 1 1; 0 1 0 1; 0 0 0 0];
m=[1 1 1 1];
q=[1 1 1 1];
k=[1 2 3 4; 2 4 1 3; 194.9068812 194.9068812 194.9068812 194.9068812; 0.9 0.9 0.9 0.9];
[Fr, Dr]=VibraStates(xyz, q, m, k)
C=SpecificHeatVib(Fr, [10000; 20000; 30000; 50000; 100000; 200000])
I=IR_Spectra(Fr, Dr, q, m, [10000; 20000; 30000; 50000; 100000; 200000])

xyz=[0 1; 0 0; 0 0];
m=[1 1];
q=[1 -1];
k=[1 ; 2 ; 143.9964485; 1.1];
[Fr, Dr]=VibraStates(xyz, q, m, k)
C=SpecificHeatVib(Fr, [10000; 20000; 30000; 50000; 100000; 200000])
I=IR_Spectra(Fr, Dr, q, m, [10000; 20000; 30000; 50000; 100000; 200000])