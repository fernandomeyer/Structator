TESTINDEX:=./testdata/index/test
TEMPINDEX:=./testdata/index/temp
VERSION:=Structator1.1

all:
	mkdir ./bin/
	cd ./construct/ && make && mv afconstruct ../bin/
	cd ./search/ && make && mv afsearch ../bin/

tar:./bin/afconstruct ./bin/afsearch
	mkdir ./$(VERSION)/
	mkdir ./$(VERSION)/bin/
	mkdir ./$(VERSION)/testdata/
	mkdir ./$(VERSION)/testdata/index/
	cp ./bin/afconstruct ./$(VERSION)/bin/
	cp ./bin/afsearch ./$(VERSION)/bin/
	cp $(TESTINDEX).aflk ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).aflkr ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).alph ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).base ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).des ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).lcp ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).lcpr ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).seq ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).skp ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).skpr ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).suf ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).sufr ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).tseq ./$(VERSION)/testdata/index/
	cp ./testdata/dna.alphab ./$(VERSION)/testdata/
	cp ./testdata/dna_rna.comp ./$(VERSION)/testdata/
	cp ./testdata/indexed.out ./$(VERSION)/testdata/
	cp ./testdata/online.out ./$(VERSION)/testdata/
	cp ./testdata/patterns.pat ./$(VERSION)/testdata/
	cp ./testdata/RF00193.pat ./$(VERSION)/testdata/
	cp ./testdata/RF00635.pat ./$(VERSION)/testdata/
	cp ./testdata/rna.alphab ./$(VERSION)/testdata/
	cp ./testdata/test.fa ./$(VERSION)/testdata/
	cp ./testdata/test ./$(VERSION)/testdata/
	cp -R ./doc/ ./$(VERSION)/
	tar zcvf $(VERSION).tar.gz ./$(VERSION)/
	rm -R ./$(VERSION)/
	
test:./bin/afconstruct ./bin/afsearch ./testdata/test
	cd ./testdata/ && chmod 755 ./test && ./test

clean:
	cd ./construct/ && make clean
	cd ./search/ && make clean
	rm -R -f ./bin/ ./$(VERSION)/
	rm -f $(TEMPINDEX).* ./testdata/temp_*

