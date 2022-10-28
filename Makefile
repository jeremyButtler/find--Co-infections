all: fastqGrep cluster scoreReads

fastqGrep:
	$(MAKE) -C source/fastqGrepSrc;
	mv source/fastqGrepSrc/fastqGrep coInfectScripts/fastqGrep;
	chmod a+x coInfectScripts/fastqGrep;

cluster:
	$(MAKE) -C source/clusterSrc;
	mv source/clusterSrc/cluster coInfectScripts/cluster;
	chmod a+x coInfectScripts/cluster;

scoreReads:
	$(MAKE) source/scoreReadsSrc;
	mv source/scoreReadsSrc/scoreReads coInfectScripts/scoreReads;
	chmod a+x coInfectScripts/scoreReads;

egcc: fastqGrepEgcc clusterEgcc scoreReadsEgcc

fastqGrepEgcc:
	$(MAKE) egcc -C source/fastqGrepSrc;
	mv source/fastqGrepSrc/fastqGrep coInfectScripts/fastqGrep;
	chmod a+x coInfectScripts/fastqGrep;

clusterEgcc:
	$(MAKE) egcc -C source/clusterSrc;
	mv source/clusterSrc/cluster coInfectScripts/cluster;
	chmod a+x coInfectScripts/cluster;

scoreReadsEgcc:
	$(MAKE) egcc -C source/scoreReadsSrc;
	mv source/scoreReadsSrc/scoreReads coInfectScripts/scoreReads;
	chmod a+x coInfectScripts/scoreReads;
