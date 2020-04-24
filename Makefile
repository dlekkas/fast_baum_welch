OUTDIR = bin
OBJDIR = obj
SRCDIR = src

INCLUDEDIR = include
CC = /usr/bin/g++

CFLAGS = -c -std=c++17 -Wall -O3

DBGFLAGS = -g -DDEBUG


all: $(OUTDIR)/baum_welch

$(OUTDIR)/baum_welch: $(OBJDIR)/bw.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o
	$(CC) $(OBJDIR)/bw.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o -o $@

$(OBJDIR)/test.o: $(SRCDIR)/test.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/test.cc

$(OBJDIR)/bw.o: $(SRCDIR)/bw.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/bw.cc

$(OBJDIR)/hmm.o: $(SRCDIR)/hmm.cc $(INCLUDEDIR)/hmm.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/hmm.cc

$(OBJDIR)/benchmark.o: $(SRCDIR)/benchmark.cc $(INCLUDEDIR)/benchmark.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/benchmark.cc


debug: $(OUTDIR)/baum_welch_dbg

$(OUTDIR)/baum_welch_dbg: $(OBJDIR)/bw_dbg.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o
	$(CC) $(OBJDIR)/bw_dbg.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o -o $@

$(OBJDIR)/bw_dbg.o: $(SRCDIR)/bw.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(DBGFLAGS) $(CFLAGS) $(SRCDIR)/bw.cc

clean:
	rm -f $(OBJDIR)/*.o $(OUTDIR)/*
