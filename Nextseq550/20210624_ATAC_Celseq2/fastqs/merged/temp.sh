#!bin/bash
f=F63M_PBS.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F63M_IFN.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F61K_PBS.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F63K_PBS.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F22F_PBS.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F62K_PBS.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F61K_IFN.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F63K_IFN.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F62M_PBS.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F62M_IFN.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F22F_IFN.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
f=F62K_IFN.p1.fq.gz; zcat ${f} | wc -l | awk -v n=${f} '{print n,$1/4}';
