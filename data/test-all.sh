#test-all.sh

date
   ../bin/MakeArray ATest1; ../bin/MakeArray ATest2
date

date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 -a 59001 -b 62000 -f -g U89959 \
	> out1.gbS
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest2 -a 59001 -b 62000 -f -g U89959 \
	> out2.gbS
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 ATest2 -a 59001 -b 62000 -f \
	-g U89959  > out12.gbS
date
date
   ../bin/GeneSeqer -h -s Arabidopsis -d ATest1 ATest2 -a 59001 -b 62000 -f \
	-g U89959  > out12.gbS.html
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 ATest2 -g U89959 > out.gbA
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 ATest2 -o outfile.gbA \
	-g U89959 > stdout.gbA
date
date
   ../bin/GeneSeqer -h -s Arabidopsis -d ATest1 ATest2 -o outfile.gbA.html \
	-g U89959 > stdout.gbA.html
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 ATest2 -w 0.8 -g U89959 > out.gbAw8
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 ATest2 -M 60000 -g U89959 > outM60.gbA
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 ATest2 -R -l L89959 > out.faA
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 ATest2 -x 30 -y 45 -z 60 \
	-l L89959 > out.faAX
date
date
   ../bin/GeneSeqer -s Arabidopsis -d ATest1 ATest2 -x 30 -y 45 -z 60 \
	-p prmfileHQ -l L89959 > out.faAXHQ
date
date
   ../bin/GeneSeqer -s Arabidopsis -E est3 -g U89959 > out.est3
date
date
   ../bin/GeneSeqer -s rice -e rice_ests -g U89959 > out.rice
date
date
   ../bin/GeneSeqer -s rice -e rice_ests -p prmfileU -g U89959 > out.riceU
date
date
   ../bin/GeneSeqer -s Arabidopsis -q QP-S23728 -a 72001 -b 76000 -r -g \
	AC002396 > out.qp1
date
date
   ../bin/GeneSeqer -h -s Arabidopsis -q QP-S23728 -a 72001 -b 76000 -r -g \
	AC002396 > out.qp1.html
date
date
   ../bin/GeneSeqer -s Arabidopsis -Q QP-U89959 -a 101 -b 13000 \
	-g U89959 > out.qp2
date
date
   ../bin/GeneSeqer -h -s Arabidopsis -d ATest1 ATest2 -Q QP-U89959 -a 101 -b 13000 \
	-g U89959 > out.estqp2.html
date
date
   ../bin/GeneSeqer -s Arabidopsis -Q QP-U89959 -I ../include/PPAM250 \
	-a 101 -b 13000 -g U89959 > out.qp3
date
date
   ../bin/SplicePredictor -s Arabidopsis -p 1 -g U89959 > sp.out
date
date
   ../bin/SplicePredictorLL -s Arabidopsis -p 1 -g U89959 > spLL.out
date
