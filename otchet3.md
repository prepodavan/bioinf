## Otchet 3
### sum_arguments.sh
```bash
$ touch sum_arguments.sh && chmod u+x sum_arguments.sh && cat > sum_arguments.sh <<EOL
#!/bin/bash
echo \$((\$1 + \$2))
EOL
$ ./sum_arguments.sh 3 4 
7
```

### count_lines.sh
```bash
$ touch count_lines.sh  && chmod u+x count_lines.sh && cat > count_lines.sh <<EOL
#!/bin/bash


SRC="/tmp/ai24_term1_pr3_files"
FILES="\$(ls \$SRC/*.txt)"
for f in \$FILES
do
        echo "Имя файла: \$f"
        echo "Количество строк: \$(wc -l \$f | cut -d ' ' -f 1)"
        echo ""
done
EOL
$ ./count_lines.sh
Имя файла: /tmp/ai24_term1_pr3_files/1.txt
Количество строк: 9

Имя файла: /tmp/ai24_term1_pr3_files/cat.txt
Количество строк: 2

Имя файла: /tmp/ai24_term1_pr3_files/example2.txt
Количество строк: 9

Имя файла: /tmp/ai24_term1_pr3_files/example.txt
Количество строк: 9

Имя файла: /tmp/ai24_term1_pr3_files/fruits.txt
Количество строк: 10

Имя файла: /tmp/ai24_term1_pr3_files/organisms.txt
Количество строк: 23

Имя файла: /tmp/ai24_term1_pr3_files/test.txt
Количество строк: 2
```

### make_dirs.sh
```bash
$ touch make_dirs.sh  && chmod u+x make_dirs.sh && cat > make_dirs.sh <<EOL
#!/bin/bash


mkdir ./dirs
for letter in {A..H}
do
	mkdir ./dirs/old\$letter
done
EOL
$ ./make_dirs.sh
$ tree
.
├── count_lines.sh
├── dirs
│   ├── oldA
│   ├── oldB
│   ├── oldC
│   ├── oldD
│   ├── oldE
│   ├── oldF
│   ├── oldG
│   └── oldH
├── make_dirs.sh
└── sum_arguments.sh

10 directories, 3 files
```

```bash
$ touch rename_dirs.sh  && chmod u+x rename_dirs.sh && cat > rename_dirs.sh <<EOL
#!/bin/bash


for letter in {A..H}
do
	old="./dirs/old\$letter"
	new="./dirs/new\$letter"
	if [ -d "\$old" ]; then
		mv "\$old" "\$new"
	fi
done
EOL
$ tree
.
├── count_lines.sh
├── dirs
│   ├── oldA
│   ├── oldB
│   ├── oldC
│   ├── oldD
│   ├── oldE
│   ├── oldF
│   ├── oldG
│   └── oldH
├── make_dirs.sh
├── rename_dirs.sh
└── sum_arguments.sh

10 directories, 4 files
$ ./rename_dirs.sh
$ tree
.
├── count_lines.sh
├── dirs
│   ├── newA
│   ├── newB
│   ├── newC
│   ├── newD
│   ├── newE
│   ├── newF
│   ├── newG
│   └── newH
├── make_dirs.sh
├── rename_dirs.sh
└── sum_arguments.sh

10 directories, 4 files
```

### range_check.sh
```bash
$ touch range_check.sh  && chmod u+x range_check.sh && cat > range_check.sh <<EOL
#!/bin/bash


read -p "Enter number: " number
if ((number >= 1 && number <= 100)); then
	echo "valid"
else
	echo "invalid"
fi
EOL
$ /bin/bash range_check.sh
Enter number: 33
valid
$ /bin/bash range_check.sh
Enter number: -10
invalid
```

### compile.sh
```bash
$ touch compile.sh  && chmod u+x compile.sh && cat > compile.sh <<EOL
#!/bin/bash


for f in \$(ls ~/Bioinf_term1/practice3/*.sh)
do
	echo "= > \$f < ="
	cat \$f
	echo ""
done
EOL
$ ./compile.sh > otchet3.txt
```