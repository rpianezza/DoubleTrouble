03-redo-gypsy-7-sim1
================
roko
8/29/2023

# Intro

Previously we realized that gypsy-7-sim1 is good representative sequence
of gypsy-7. So we need to redo the repeatmasking with gypsy-29 and
gypsy-7-sim1

# Prepare reference

``` bash
cat ../rawseqs-v2/gypsy-ext.fa|fasta-reader.py |fasta-subsequences.py --fasta-ids gypsy-29,gypsy-7-sim1 |fasta-writter.py > gypsy-final.fa
# I truly start to love my bioawk library ;)
```

# RepeatMask

``` bash
for i in *.fa; do RepeatMasker -pa 20 -no_is -s -nolow -dir out-gypsy-v3 -lib replib-gypsy-v3/gypsy-final.fa $i;done > /dev/null 

# selfmask
RepeatMasker -pa 20 -no_is -s -nolow -dir selfmask -lib gypsy-final.fa gypsy-final.fa   
```

# Prepare

``` bash
awk '{print $0,"self"}' selfmask/gypsy-final.fa.ori.out |perl -pe 's/\.fa\.ori\.out//' >tmp-merged-self.sum
```
