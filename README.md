#abismal
./abismalidx chr6.fa index38_6

parallel-fastq-dump --sra-id SRR1042908 --threads 46 --outdir . --split-files

/labs/mpsnyder/moqri/soft/abismal/bin/abismal -i /labs/mpsnyder/moqri/soft/abismal/bin/index38_6 -o SRR1042914_6.sam SRR1042914_1.fastq SRR1042914_2.fastq -t 46

#methypipe
../configure CPPFLAGS='-I/labs/mpsnyder/moqri/soft/htslib/install/include -I/scg/apps/software/gsl/2.6/include' LDFLAGS='-L/labs/mpsnyder/moqri/soft/htslib/install/lib -L/scg/apps/software/gsl/2.6/lib'
#openssl
./config --prefix=/labs/mpsnyder/moqri/soft/openssl --openssldir=/labs/mpsnyder/moqri/soft/openssl no-ssl2


export PATH=/labs/mpsnyder/moqri/soft/openssl/bin:$PATH
export LD_LIBRARY_PATH=/labs/mpsnyder/moqri/soft/openssl/lib
export LC_ALL="en_US.UTF-8"
export LDFLAGS="-L /labs/mpsnyder/moqri/soft/openssl/lib -Wl,-rpath,/labs/mpsnyder/moqri/soft/openssl/lib"

module load samtools
samtools sort -O sam -o SRR1042908_6s.sam SRR1042908_6.sam -@ 16 -m 4000000

/labs/mpsnyder/moqri/soft/methpipe-5.0.0/inst/bin/duplicate-remover -S stats SRR1042908_6s.sam SRR1042908_6sd.samm

/labs/mpsnyder/moqri/soft/methpipe-5.0.0/inst/bin/methcounts -c /labs/mpsnyder/moqri/soft/abismal/bin/chr6.fa -o SRR1042908.meth /labs/mpsnyder/moqri/data/meth/wg/skin/SRR1042908_6sd.sam -v -n
