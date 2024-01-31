# pig-nosed-turtle
A chromosome-level genome assembly of the pig-nosed turtle (Carettochelys insculpta)
## soapec.sh
``` bash
/public/home/zhuchenglong/software/soapec/SOAPec_bin_v2.03/bin/KmerFreq_HA -k 21 -t 50 -p PigNoseTurtle -l reads.list -L 150 >soapec.kmerfreq.log 2>soapec.kmerfreq.err
/public/home/zhuchenglong/software/soapec/SOAPec_bin_v2.03/bin/KmerFreq_HA -k 23 -t 50 -p PigNoseTurtle -l reads.list -L 150 >soapec.kmerfreq.log 2>soapec.kmerfreq.err
/public/home/zhuchenglong/software/soapec/SOAPec_bin_v2.03/bin/KmerFreq_HA -k 25 -t 50 -p PigNoseTurtle -l reads.list -L 150 >soapec.kmerfreq.log 2>soapec.kmerfreq.err
/public/home/zhuchenglong/software/soapec/SOAPec_bin_v2.03/bin/KmerFreq_HA -k 27 -t 50 -p PigNoseTurtle -l reads.list -L 150 >soapec.kmerfreq.log 2>soapec.kmerfreq.err
```
