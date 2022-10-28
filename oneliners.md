### convert NCBI genome print out to fasta
 1. cut and paste into text file
 ```
        1 gggggaggag ccaagatggc cgaataggaa cagctccggt ctacagctcc cagcgtgagc
       61 gacgcagaag acggtgattt ctgcatttcc atctgaggta ccgggttcat ctcactaggg
      121 agtgccagac agtgggcgca ggccagtgtg tgtgcgcacc gtgcgcgagc cgaagcaggg
      181 cgaggcattg cctcacctgg gaagcgcaag gggtcaggga gttccctttc tgagtcaaag
 ```
 2. and add header line
 ```
 >L19088 Human LINE1 (L1.3) repetitive element DNA sequence. https://www.ncbi.nlm.nih.gov/nuccore/L19088.1
 ```
 3 run oneliner
```
awk '{ if (NR > 1) { astr = $2 $3 $4 $5 $6 $7; printf("%s",toupper(astr)); } else { print $0 } }' L1.3 > L1.3.fa
```

### fastq to fasta
 awk '{ if ( NR % 4 == 1 ) print ">" $0; if ( NR % 4 == 2 ) print $0; }' file.fastq > out.fasta

### Get the first region from a fasta file.

awk '{ if (($1 ~ /^>/) && ( NR > 1 )) { exit } print $0 }' hap1.fa | less -S

Get the first region from a fasta file.
  ($1 ~ /^>/) If the first field begins with '>'
  NR > 1      if the current line read from the file is greater than 1
  exit        This will exit at the start of the second region 
  print $0    Print the first line, and subsequent lines if they do not begin with '>'
  hap1.fa     The file we are looking at.
  | less -S   Pipe results to less, -S option prevents word wrap and allows horizontal scrolling with arrow keys.
              Alternately redirect output(the first region) to file with "> hap1_region1.fa"
 
#### All region and count of nucleotide: (first line blank)
gunzip -c hap1.fa.gz | awk '{ if ( $1 ~ /^>/ ) { printf("%s\n%s, ",nCount,$0); nCount = 0 } else { nCount = nCount + length($0) } } END { print nCount }' | sort -k5,5 -n | less -S

#### List nucleotide count for each region:
* gunzip -c hap1.fa.gz | awk '{ if ( $1 ~ /^>/ ) { if ( NR > 1 ) { printf("%s\n",nCount); nCount = 0 } } else { nCount = nCount + length($0) } } END { print nCount }' | sort -k1,1 -n | less -S

  * + *region and count output:* awk '{ if ( $1 ~ /^>/ ) { if ( NR > 1 ) { printf("%s\n",nCount); nCount = 0; print $1 } } else { nCount = nCount + length($0) } } END { print nCount }' heart.all_size.5merge.collapsed.longest_rep.fa |  less

  * + *Each count and name on lone line and sorted by size:* awk '{ if ( $1 ~ /^>/ ) { if ( NR > 1 ) { printf("%s : %s\n",nCount,lastR); nCount = 0;} lastR = $1 }  else { nCount = nCount + length($0) } } END { printf("%s : %s\n",nCount,lastR) }' heart.all_size.5merge.collapsed.longest_rep.fa |  sort -n | less

  * + *'size tab name' => "23    ALU" on each line:*  awk '{ if ( $1 ~ /^>/ ) { if ( NR > 1 ) { printf("%s\t%s\n",nCount,lastR); nCount = 0;} lastR = substr($1,2,length($1)-1) } else { nCount = nCount + length($0) } } END { printf("%s\t%s\n",nCount,lastR) }'  ../biolib/humsubrep.ref | less

### Get specified region from a fasta file:
By name (or part name) [will exit after first region match]
```
awk -v reg="1260"  '{ if (($1 ~ /^>/) && ($1 ~ reg)) { fnd = 1; aRow = NR } if ( fnd == 1 ) { if (($1 ~ /^>/) && (NR > aRow)) { exit } print $0 }}' Odioica_reference_v3.0.fa | less -S
```

#### Get region number from a fasta file:
```
awk -v reg="5"  '{ if ($1 ~ /^>/) { regCnt += 1} if (regCnt == reg) { fnd = 1; aRow = NR } if ( fnd == 1 ) { if (($1 ~ /^>/) && (NR > aRow)) { exit } print $0 }}' Odioica_reference_v3.0.fa | less -S
```

### Scaffold/Contig/Region sizes in a fasta assembly file.

awk '{ if ( $0 ~ /^>/ ) { if (a > 0) { printf("%s\n%s : ",a,$0); a = 0 } else { printf("%s : ",$0) } }  else  { a = a + length($0) } } END { print a } ' Odioica_reference_v3.0.fa | less

Ouput Sample:
```
>scaffold_1 : 3167015
>scaffold_2 : 2850711
>scaffold_3 : 2515153
>scaffold_4 : 2289853
>scaffold_5 : 1792778
>scaffold_6 : 1516148
>scaffold_7 : 1487711
.....
```

### Split fasta:

output each with 1.fa, 2.fa 3.fa etc... </br>
awk '/^>/{s=++d".fa"} {print > s}' multi.fa  </br>
( re: https://awesomeopensource.com/project/crazyhottommy/bioinformatics-one-liners )

Make it more useful: </br>
```
awk '/^>/{s=++d".fa"; fname = substr($1,2,length($1)); n = index(fname,"GL");   if (n == 0) {print fname "-" $2 ".fa";} } ' hs37d5.fa
awk '/^>/{s=++d".fa"; fname = substr($1,2,length($1)); n = index(fname,"GL");   if ($1 !~ /GL|NC/) {print fname "-" $2 ".fa";} } ' hs37d5.fa
```

FINAL:
Split a genome into a file for each of its chromosomes, or filter for chromosomes required.

First test that the generated filenames will be suitable and unique, and filtered are as you require: (adjust $1 $2 $3 etc and specify filter for 
 chromosomes not required eg: 'GL' 'NC' [add another with '|MT' ] ) 
``` 
awk '/^>/{fname = substr($1,2,length($1)) "-" $2 ".fa"; if ($1 !~ /GL|NC/) {faOut = 1} else {faOut = 0} if (faOut == 1) print fname }' hs37d5.fa
```

Now update this command to output chromosomes requried to the filenames, based on your filename requirements, and this will create the output file 
for each chromosome:
```
awk '/^>/{fname = substr($1,2,length($1)) "-" $2 ".fa"; if ($1 !~ /GL|NC/) {faOut = 1} else {faOut = 0} } {if (faOut == 1) {print > fname} }' hs37d5.fa
```
hs37d5.fa : the  genome name.

$1 $2 $3 ... are the space delimited column data from the current line.  AWK process line by line and extract fields to $1, $2 etc. Special case $0 is 
  the whole line (all fields)
  
This part : {fname = substr($1,2,length($1)) "-" $2 ".fa"; if ($1 !~ /GL|NC/) {faOut = 1} else {faOut = 0} }  will only be processed if the  
 current line for the genome begins with ">"
 
fname => create a filename from the chromosome line data in the fasta file, (the substr functions removes the ">" from the start of the first field) 
  build you filename based on the data in your chromosome /region  identifier, ie $3 could be a better filename, make sure it is unique, otherwise
   you will write over previous output.
 
the first 'if' will not print any chromosome where the first field contians  a GL or a NC  (substitute for your chromosome/region not required)
 
The final if will print to fname if faOut is set to 1. ( here the '>' is a redirect to file not a greater than)

Maybe alter above based on this?  
```
awk ' { if ( $0 ~ /^>/ ) { if (i == 1) { exit }  i = i + 1; } print $0 }' Odioica_reference_v3.0.fa > 0.fa
```

### fastq.gz to fasta:
This is the fastest bash one-liner to convert fastq to fasta, assuming 4 lines per FASTQ record, more on https://github.com/stephenturner/oneliners

sed -n '1~4s/^@/>/p;2~4p' test.fastq > test.fasta


### count N's fasta in 80 cols or ??  under construction...
awk -FN ' BEGIN { i = 0; j = 0; bg = 0; ed = 0; } {  if ( $1 ~ /^>/  ) {if ( i > 1 )  { ed = j - 80; print bg,ed, (ed - bg), "line ", NR } j = j - 80; i = 0; print $0 } j = j + 80; if ( ( NF-1 ) == 80 ) { i = i + 80 } else { if ( i > 1 )  { ed = j - 80; i = 0; print bg,ed, (ed - bg), " line ", NR } else { i = 0 } } if ( i == 80 ) { bg = j - 80; print NR } } END { }' GCA_000001405.28_GRCh38.p13_genomic.fna > testN.txt

Better : https://superuser.com/questions/485800/whats-the-quickest-way-to-count-the-number-of-each-character-in-a-file/485811  
echo 'int cache[256],x,y;char buf[4096],letters[]="tagcnTAGCN-"; int main(){while((x=read(0,buf,sizeof buf))>0)for(y=0;y<x;y++)cache[(unsigned char)buf[y]]++;for(x=0;x<sizeof letters-1;x++)printf("%c: %d\n",letters[x],cache[letters[x]]);}' | gcc -w -xc -; ./a.out < chr6GRCh38p13.fna; rm a.out;

#### Split a Fasta (ref/contig/scaffold) file into reads
echo 'int x,y,d;char buf[4096];int main() {while ((x=read(0,buf,sizeof buf))>0) { if (d > 0) {printf(">R%d\n",d); for(y=0;y<x;y++) printf("%c",buf[y]); printf("\n");} d++;} }' | gcc -w -xc -; ./a.out < hap1-9.fa | less -S

read size = 4096
first read will read line 1 of file = scaffold name we skip this. 
Reads area labelled R1 ... Rn     n will be how ever many reads are split up.

### View a fasta file which has fixed line width, to one line of sequence per read/id
awk '{ if ( substr($0,1,1) == ">" ) printf("\n%s\n",$0); else printf("%s",$0) }' /home/humrep/repbase_hum/humsub-labelDup.ref | less -S

### work with samtools and gzipped file
index fastq:
 gzip -d -c your.fastq.gz | samtools faidx - -o your.index    The '-' will make samtools read from stdin.


### subsample fastq.gz:

Get every fourth read: ( IE divide by 4)

gunzip -c /reads/illumina.fastq.gz | awk '{ if ((NR % 16) == 1) {ot = 1}; if (ot ==1) { print $0; a = a + 1 }; if( a == 4) { a = 0; ot = 0}; }' | gzip > IlluminaDiv4.fastq.gz

IF PAIRED reads and pairs are consecutive:   (some tools need /1 /2 or .1 .2 removed! )

gunzip -c /reads/illumina.fastq.gz | awk '{ if ((NR % 32) == 1) {ot = 1}; if (ot ==1) { print $0; a = a + 1 }; if( a == 8) { a = 0; ot = 0}; }' | gzip > IlluminaDiv4.fastq.gz

### SUB SAMPLE From BAM
IF YOU NEED TO KEEP READ PAIRS use SAMTOOLS "-s":

samtools view -bs 42.0005 $illumina
```ruby
samtools has a subsampling option: 

-s FLOAT: 
 Integer part is used to seed the random number generator [0]. 
 Part after the decimal point sets the fraction of templates/pairs to subsample [no subsampling] 

samtools view -bs 42.1 in.bam > subsampled.bam 
 
will subsample 10 percent mapped reads with 42 as the seed for the random number generator. 
```

### get large inserts from VCF file:
less -S HG002.m84005_220919_232112_s2.GRCh38.pbsv.vcf.gz | awk 'BEGIN {rx="SVLEN.[0-9][0-9][0-9][0-9]"} {if ($0 ~ rx) print }' | less -S
 note the regex rx does not need to be inside // in the comparison $0 ~ rx  [not $1 ~ /rx/ ] 


### Find unique ? under development ( ie from two files output from above "inserts from VCF file" )
awk '{ ind = $1 ":" $2; if (FNR == NR) { k[ind]=1; } else { if (ind in k) {k[ind]=k[ind]+1;} else {k[ind]=1;}  }'
if k[ind] = 1 then it is unique.  But need a window typqe function too?

ie two vcf's find  non duplicates (or duplicates k[i] ==2)
```ruby
awk '{ ind = $1 ":" $2; if (ind in k) {k[ind]=k[ind]+1;} else {k[ind]=1;} } END { for ( i in k) {if (k[i] == 1) print i}  }'  HG02_S2.vcf HG02_S1.vcf | less -S
```
=> neeed FNR == NR to identify file containing unique.
Add file id: 
```ruby
 awk '{ ind = $1 ":" $2; if (FNR == NR) {file = "A"} else {file = "B"} kf[ind]=file; if (ind in k) {k[ind]=k[ind]+1;} else {k[ind]=1;} } END { for ( i in k) {if (k[i] == 1) print i, kf[i]}  }'  HG02_S2.vcf HG02_S1.vcf |  sort -V | less -S
```

```ruby 
chr1:191341 A
chr1:977117 A
chr1:977247 B
chr1:1289035 A
chr1:1929384 A
chr1:2663027 A
chr1:2675143 B
```
