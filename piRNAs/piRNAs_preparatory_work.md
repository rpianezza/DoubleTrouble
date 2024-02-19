piRNAs preparatory work
================
Almorò Scarpa

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.4.0     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.5.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

Remove the adapter sequence and put then new files in a folder

``` bash
for file in *.fastq.gz; do filename=$(basename "$file" .fastq); cutadapt -j 10 -a TGGAATTCTCGGGTGCCAAGG -o "./noadapt/$filename" "$file"; done
```

In parallel remove the extra column and name the new files .clean:

``` bash
ls *.fastq.gz | parallel 'output_file=$(basename {} .fastq.gz); gzip -cd {} | awk "{print \$1}" | gzip -c > "$output_file.clean.fastq.gz"'
```

Size filter, this will keep only small RNAs:

``` bash
ls *.clean.fastq.gz | parallel 'output_file=$(basename {} .clean.fastq.gz); gzip -cd {} | paste - - - - | awk "length(\$2) > 17 && length(\$2) < 36" | tr "\t" "\n" | gzip -c > "$output_file.trim.fastq.gz"'
```

Novoalign to index the fast file with all the RNAs:

``` bash
novoindex newfilename.nvi fastafile
```

The fasta files with all the transcript, non-coding RNAs and miRNAs can
be dowloaded from Flybase But they need to be modified changing the
header, keeping the ID and adding “\_type”, type can be either mRNA,
miRNA or ncRNA. Then cat with a te library, TEs header must end with
“\_te”

Map wit Novoalign:

``` bash
for file in *.trim.fastq.gz; do filename=$(basename "$file" .trim.fastq.gz); gzip -cd "$file" | /Applications/novocraft/novoalign -d /path/to/folder/TEs_and_RNAs.nvi -f - -F STDFQ -o SAM -o FullNW -r RANDOM > "${filename}.sam"; done
```

Make directory bam:

``` bash
mkdir bam
```

Turn .sam into bam:

``` bash
for i in *.sam ; do n=${i%.sam} ; samtools view -Sb $i > bam/$n.bam ; done
```

Sort the .bam file:

``` bash
for i in *.bam ; do n=${i%.bam} ; samtools sort $i -o $n.sort ; done
```

piRNAs on TEs:

``` bash
for i in *.sort; do n=${i%.sort}; samtools view $i |python /path/to/folder/piRNA-distribution-onTE.py --sam - --sample-id $n > $n-graph-on-TE.forR ; done
```

Ping-pong:

``` bash
for i in *.sort; do n=${i%.sort}; samtools view $i | python /path/to/folder/ping-pong-signature.py --sam - --max-mm 2 --sample-id $n > $n.pps ; done
```
