## Otchet
```bash
$ file -i otchet.md
otchet.md: regular file
```
### 1. Создание и работа с файлами и директориями


1. Создадим сразу обе директории для домашки через `-p`
```bash
mkdir -p Bioinf_term1/practice1
```


2. Перейдите в practice1 и проверьте путь до рабочей директории.
Для этого можно использоваться пвд
```bash
cd Bioinf_term1/practice1 && pwd
```
outputs: `/home/students/ai24/gind.alex/Bioinf_term1/practice1`


3. Создать три файла и проверить. Для проверки можно использовать cat чтобы воспользоваться глобом.
```bash
$ echo "Это первый файл" > fle1.txt
$ echo "Это второй файл" > fle2.txt
$ echo "Это третий файл" > fle3.txt
$ cat fl*
Это первый файл
Это второй файл
Это третий файл
```


4. Создайте вложенные директории data/raw
Воспользуемся -p чтобы два раза не вставать
```bash
$ mkdir -p data/raw
$ tree
.
├── data
│   └── raw
├── fle1.txt
├── fle2.txt
├── fle3.txt
├── otchet.md
└── otchet.txt -> otchet.md

3 directories, 5 files
```


5. Переместите все файлы созданные в задании 1.3 в директорию data/raw одной командой.
Воспользуемся глобом и для перемещения и для проверки катом
```bash
gind.alex@kodomo:~/Bioinf_term1/practice1$ mv fl* ./data/raw/
gind.alex@kodomo:~/Bioinf_term1/practice1$ tree
.
├── data
│   └── raw
│       ├── fle1.txt
│       ├── fle2.txt
│       └── fle3.txt
├── otchet.md
└── otchet.txt -> otchet.md

3 directories, 5 files
gind.alex@kodomo:~/Bioinf_term1/practice1$ cat data/raw/fl*
Это первый файл
Это второй файл
Это третий файл
```


6. Используя относительный путь, перейдите в директорию raw . Предложите вариант перехода, используя абсолютный путь.
```bash
$ cd ./data/raw
$ pwd
/home/students/ai24/gind.alex/Bioinf_term1/practice1/data/raw
```

через абсолютный тоже можно воспользоваться пвд, так как в задании не упоминалось, что я должен ввести абсолютный путь вручную
`cd $(pwd)/data/raw`


7. Соберите строки всех файлов в один файл
`cat *.txt > all.txt`


8. Просмотрите первые две строки файла all.txt
```bash
$ head -n 2 all.txt
Это первый файл
Это второй файл
```


9. Подсчитайте количество строк, слов и символов в файле all.txt
```bash
$ wc -m -l -w all.txt
3  9 48 all.txt
```


10. Вернитесь в домашнюю директорию
`cd ~`


### 2. Права доступа файлов и директорий


1. Просмотрите права доступа к директории practice1 и ее содержимому.
```bash
$ ls -lah Bioinf_term1
total 4.0K
drwxr-xr-x.  3 gind.alex ai-24   23 Oct  9 16:35 .
drwxr-x---+ 14 gind.alex ai-24 4.0K Oct  9 17:11 ..
drwxr-xr-x.  3 gind.alex ai-24   75 Oct  9 17:00 practice1
```


2. Измените права доступа к файлу all.txt, чтобы только владелец мог его читать и писать, а остальные пользователи не имели доступа.
6 - чтение и запись
0 - запрет на все операции
сначала владелец, потом группа, потом остальные
```bash
$ chmod 0600 ./Bioinf_term1/practice1/data/raw/all.txt
```


3. Проверьте, что права доступа изменились.
проверим, что права изменились у нужного файла лсом по файлу. лсом по папке практике проверим, что права не сделались рекурсивно у всех файлов и папок
```bash
$ ls -lah Bioinf_term1/practice1/data/raw/all.txt
-rw-------. 1 gind.alex ai-24 87 Oct  9 17:05 Bioinf_term1/practice1/data/raw/all.txt
$ ls -lah Bioinf_term1/practice1
total 20K
drwxr-xr-x. 3 gind.alex ai-24  75 Oct  9 17:00 .
drwxr-xr-x. 3 gind.alex ai-24  23 Oct  9 16:35 ..
drwxr-xr-x. 3 gind.alex ai-24  17 Oct  9 16:37 data
-rw-r--r--. 1 gind.alex ai-24   4 Oct  9 16:41 otchet.md
-rw-r--r--. 1 gind.alex ai-24 16K Oct  9 17:18 .otchet.md.swp
lrwxrwxrwx. 1 gind.alex ai-24   9 Oct  9 16:40 otchet.txt -> otchet.md
```


### 3. Архивация и сжатие
1. Создайте архив practice1.tar из директории practice1
`tar -c -f practice1.tar Bioinf_term1/practice1`


2. Сожмите архив с помощью gzip.
`gzip -k practice1.tar`
-k нужен, чтобы не удалять текущий файл


3. Создайте директорию archive в домашней директории и скопируйте practice1.tar.gz в нее.
`mkdir archive && cp practice1.tar.gz archive/practice1.tar.gz`
проверим, что файлы одинаковые
```bash
$ md5sum practice1.tar.gz
c4ac4a18850f5ab82528ce671916b536  practice1.tar.gz
$ md5sum archive/practice1.tar.gz
c4ac4a18850f5ab82528ce671916b536  archive/practice1.tar.gz
```


4. Перейдите в директорию ~/archive и распакуйте архив
```bash
$ cd archive/ && gunzip practice1.tar.gz && tar -x -f practice1.tar && tree
.
├── Bioinf_term1
│   └── practice1
│       ├── data
│       │   └── raw
│       │       ├── all.txt
│       │       ├── fle1.txt
│       │       ├── fle2.txt
│       │       └── fle3.txt
│       ├── otchet.md
│       └── otchet.txt -> otchet.md
└── practice1.tar

5 directories, 7 files
```


### 4. Скачивание файлов

1. `wget 'https://today.uconn.edu/wp-content/uploads/2017/07/GettyImages-146798910-CubanTreeFrog-1024x683.jpg' -O download.jpg`
2. `scp gind.alex@kodomo:/home/students/ai24/gind.alex/Bioinf_term1/practice1/download.jpg download.jpg`
3. 
```bash
$ ls download*
download.jpg
```
