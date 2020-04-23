OUTDIR = bin
OBJDIR = obj
SRCDIR = src

INCLUDEDIR = include
CC = /usr/bin/g++
CFLAGS = -c -std=c++17 -Wall -O3 -g

all: $(OUTDIR)/baum_welch

$(OUTDIR)/baum_welch: $(OBJDIR)/bw.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o
	$(CC) $(OBJDIR)/bw.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o -o $@

$(OBJDIR)/test.o: $(SRCDIR)/test.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/test.cc

$(OBJDIR)/bw.o: $(SRCDIR)/bw.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/bw.cc

$(OBJDIR)/hmm.o: $(SRCDIR)/hmm.cc $(INCLUDEDIR)/hmm.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/hmm.cc

clean:
	rm -f $(OBJDIR)/*.o $(OUTDIR)/*
