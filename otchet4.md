## Otchet 4
### sed
```bash
cp /tmp/reads.fasta reads.fasta
sed -r 's/T/U/g' < reads.fasta >  rna_sequence.txt
sed -ri 's/T/U/3g' reads.fasta
sed -ri 's/read/seq_read/g' reads.fasta
```

### awk
```bash
$ awk -F ',' < score.csv '{ if ($3 == 5) {print} }'
Alice,22,5
Charlie,23,5
Eve,20,5
Grace,18,5
$ awk -F ',' < score.csv '{ if ($2 > 20 ) {print} }' | wc  -l
5
$ sed < score.csv -r "s/,/$(printf '\t')/g" | tee score.tsv | awk -F "$(printf '\t')" '{ if ($2 > 20 && $3 == 5) {print} }'
Alice	22	5
Charlie	23	5
$ cat score.tsv
Name	Age	Score
Alice	22	5
Bob	19	4
Charlie	23	5
David	21	3
Eve	20	5
Frank	25	4
Grace	18	5
```