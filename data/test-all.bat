REM test-all.bat

time /T
   ..\bin\MakeArray ATest1
   ..\bin\MakeArray ATest2
time /T

time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 -a 59001 -b 62000 -f -g U89959 > out1.gbS
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest2 -a 59001 -b 62000 -f -g U89959 > out2.gbS
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 ATest2 -a 59001 -b 62000 -f -g U89959  > out12.gbS
time /T
time /T
   ..\bin\GeneSeqer -h -s Arabidopsis -d ATest1 ATest2 -a 59001 -b 62000 -f -g U89959  > out12.gbS.html
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 ATest2 -g U89959 > out.gbA
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 ATest2 -o outfile.gbA -g U89959 > stdout.gbA
time /T
time /T
   ..\bin\GeneSeqer -h -s Arabidopsis -d ATest1 ATest2 -o outfile.gbA.html -g U89959 > stdout.gbA.html
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 ATest2 -w 0.8 -g U89959 > out.gbA8w
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 ATest2 -M 60000 -g U89959 > outM60.gbA
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 ATest2 -R -l L89959 > out.faA
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 ATest2 -x 30 -y 45 -z 60 -l L89959 > out.faAX
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -d ATest1 ATest2 -x 30 -y 45 -z 60 -p prmfileHQ -l L89959 > out.faAXHQ
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -E est3 -g U89959 > out.est3
time /T
time /T
   ..\bin\GeneSeqer -s rice -e rice_ests -g U89959 > out.rice
time /T
time /T
   ..\bin\GeneSeqer -s rice -e rice_ests -p prmfileU -g U89959 > out.riceU
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -q QP-S23728 -a 72001 -b 76000 -r -g AC002396 > out.qp1
time /T
time /T
   ..\bin\GeneSeqer -h -s Arabidopsis -q QP-S23728 -a 72001 -b 76000 -r -g AC002396 > out.qp1.html
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -Q QP-U89959 -a 101 -b 13000 -g U89959 > out.qp2
time /T
time /T
   ..\bin\GeneSeqer -h -s Arabidopsis -d ATest1 ATest2 -Q QP-U89959 -a 101 -b 13000 -g U89959 > out.estqp2.html
time /T
time /T
   ..\bin\GeneSeqer -s Arabidopsis -Q QP-U89959 -I ..\include\PPAM250 -a 101 -b 13000 -g U89959 > out.qp3
time /T
time /T
   ..\bin\SplicePredictor -s Arabidopsis -p 1 -g U89959 > sp.out
time /T
time /T
   ..\bin\SplicePredictorLL -s Arabidopsis -p 1 -g U89959 > spLL.out
time /T
