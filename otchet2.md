## Otchet 2
```bash
$ file -i otchet2.md
otchet.md: regular file
```
### 1. Разметка генов lncРНК человека в формате GTF

```bash
$ pwd
/home/students/ai24/gind.alex/Bioinf_term1/practice2
```

1. Скачайте файл
```bash
wget 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.long_noncoding_RNAs.gtf.gz' -O 'human_lncRNA_genes.gtf.gz'
```


2. Распакуйте скачанный файл
```bash
$ gunzip human_lncRNA_genes.gtf.gz && ls -lah 
total 85M
drwxr-xr-x. 2 gind.alex ai-24  36 Oct  9 20:28 .
drwxr-xr-x. 4 gind.alex ai-24  40 Oct  9 20:26 ..
-rw-r--r--. 1 gind.alex ai-24 85M Aug  4  2020 human_lncRNA_genes.gtf
```


3. Удалите первые 5 строчек с заголовком
```bash
$ head -n 5 human_lncRNA_genes.gtf
##description: evidence-based annotation of the human genome (GRCh38), version 35 (Ensembl 101) - long non-coding RNAs
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2020-06-03
$ tail -n +6 human_lncRNA_genes.gtf > human_lncRNA_genes_no_header.gtf
$ head -n 1 human_lncRNA_genes_no_header.gtf
chr1	HAVANA	gene	29554	31109	.	+	.	gene_id "ENSG00000243485.5"; gene_type "lncRNA"; gene_name "MIR1302-2HG"; level 2; hgnc_id "HGNC:52482"; tag "ncRNA_host"; havana_gene "OTTHUMG00000000959.2";
```


4. Найдите все уникальные типы  описаний (например, “exon”) в третьем столбце файла human_lncRNA_genes_no_header.gtf
так как не сказано, что вывести надо в каком-то строгом формате, то проще посчитать через `-c` у `uniq`
```bash
$ cut -f 3 human_lncRNA_genes_no_header.gtf | sort | uniq -c
 178049 exon
  17957 gene
  48684 transcript
```


5. Сохраните только 1, 3, 4, 5 и 7 столбцы в файл lncRNA_columns.tsv. Какая информация о генах осталась в этом наборе колонок?
`cut -f 1,3,4,5,7 human_lncRNA_genes_no_header.gtf > lncRNA_columns.tsv`


таблица у нас будет `seqname,feature,start,end,strand`


6. Найдите все строки, содержащие "transcript" в lncRNA_columns.tsv.
```bash
$ grep 'transcript' lncRNA_columns.tsv | wc -l
48684
```

### 2 FASTA файл с аминокислотными последовательностями

1. Извлеките все строки, начинающиеся с символа >

```bash
$ grep -e '^>' /tmp/PIWI* > seq_info.txt && head seq_info.txt
>PIWL2_DANRE Piwi-like protein 2 OS=Danio rerio
>PIWL2_XENTR Piwi-like protein 2 OS=Xenopus tropicalis
>PIWL4_HUMAN Piwi-like protein 4 OS=Homo sapiens
>PIWL2_MOUSE Piwi-like protein 2 OS=Mus musculus
>PIWL4_MOUSE Piwi-like protein 4 OS=Mus musculus
>PIWL2_HUMAN Piwi-like protein 2 OS=Homo sapiens
>PIWL1_DANRE Piwi-like protein 1 OS=Danio rerio
>PIWL1_HUMAN Piwi-like protein 1 OS=Homo sapiens
>PIWL1_MOUSE Piwi-like protein 1 OS=Mus musculus
>PIWL1_CHICK Piwi-like protein 1 OS=Gallus gallus
```

2. Создайте конвейер, который извлекает информацию о последовательностях из /tmp/PIWI_proteins.fasta, преобразует ее в нижний регистр и сохраняет в файл lowercase_seq_info.txt.
```bash
$ grep -v -e '^>' /tmp/PIWI*  | tr '[:upper:]' '[:lower:]'  > lowercase_seq_info.txt && head lowercase_seq_info.txt
mdpkrptfpsppgvirapwqqstedqsqlldqpslgrarglimpideplpgrgrafsvpg
mdptrppfrgspfhtplgvrppvletkeegphgravllprgrallgasapssdttqrdps
msgrarvkargiarspsatevgriqasplprsvdlsnneasssngflgtsristndkygi
mdpvrplfrgptpvhpsqcvrmpgcwpqaprplepawgragpagrglvfrkpedsspplq
msgrarvrargittghsarevgrssrdlmvtsaspgdseagggtsvisqpyelgvssgdg
mdpfrpsfrgqspihpsqcqavrmpgcwpqaskpldpalgrgapagrghvfgkpeepstq
mtgrararsrgrgrgqepaapgaqppvsqeaakpvvstpsegqlvgrgrqkpapgamsee
mtgrararargrargqetaqlvgstasqqpgyiqprpqpppaegelfgrgrqrgtaggta
mtgrararargrargqetvqhvgaaasqqpgyipprpqqsptegdlvgrgrqrgmvvgat
mtgrararargrppgqeaaippvgaasaqktlpshpseqrqslqpchppplteepggrgr
```

3. Последовательно используйте команды tr и cut, чтобы извлечь только имена белков из файла seq_info.txt
`tr < seq_info.txt -d '[=>=]' | cut -d ' ' -f 1`

4. Добавьте к полученному в 2.3. конвейеру команды sort и затем uniq, чтобы подсчитать количество уникальных белков

```bash
$ tr < seq_info.txt -d '[=>=]' | cut -d ' ' -f 1 | sort | uniq -c > unique_proteins.txt
$ cat unique_proteins.txt | sort
      1 PIWL1_CHICK
      1 PIWL2_DANRE
      1 PIWL2_MOUSE
      1 PIWL2_XENTR
      1 PIWL3_HUMAN
      1 PIWL4_RAT
      2 PIWL1_DANRE
      2 PIWL1_MOUSE
      2 PIWL2_HUMAN
      3 PIWL1_HUMAN
      3 PIWL4_HUMAN
      3 PIWL4_MOUSE
```

5. Используйте команду cut, чтобы извлечь названия организмов из файла seq_info.txt. 
```bash
$ head seq_info.txt  | cut -d '=' -f 2
Danio rerio
Xenopus tropicalis
Homo sapiens
Mus musculus
Mus musculus
Homo sapiens
Danio rerio
Homo sapiens
Mus musculus
Gallus gallus
```

6. К команде cut из 2.5. добавьте команду sort и затем uniq -c, чтобы подсчитать количество последовательностей для каждого организма. 
```bash
$ cut < seq_info.txt -d '=' -f 2 | sort | uniq -c  | sort | tee organism_counts.txt
      1 Gallus gallus
      1 Rattus norvegicus
      1 Xenopus tropicalis
      3 Danio rerio
      6 Mus musculus
      9 Homo sapiens
```