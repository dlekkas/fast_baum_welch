OUTDIR = bin
OBJDIR = obj
SRCDIR = src

INCLUDEDIR = include
CC = /usr/bin/g++

CFLAGS = -c -std=c++17 -Wall -O3 -DDEBUG

DBGFLAGS = -g -DDEBUG


all: $(OUTDIR)/baum_welch $(OUTDIR)/baum_welch_basic_opts $(OUTDIR)/baum_welch_baseline
#all: $(OUTDIR)/baum_welch_baseline

$(OUTDIR)/baum_welch: $(OBJDIR)/baum_baseline_impl.o $(OBJDIR)/bw.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o $(OBJDIR)/generator.o
	$(CC) $(OBJDIR)/baum_baseline_impl.o $(OBJDIR)/bw.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o $(OBJDIR)/generator.o -o $@

$(OUTDIR)/baum_welch_basic_opts: $(OBJDIR)/baum_baseline_impl.o $(OBJDIR)/bw_basic_opts.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o $(OBJDIR)/generator.o
	$(CC) $(OBJDIR)/bw_basic_opts.o $(OBJDIR)/baum_baseline_impl.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o $(OBJDIR)/generator.o -o $@


$(OUTDIR)/baum_welch_baseline: $(OBJDIR)/baum_baseline_impl.o $(OBJDIR)/bw.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o $(OBJDIR)/generator.o
	$(CC) $(OBJDIR)/baum_baseline_impl.o $(OBJDIR)/bw.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o $(OBJDIR)/generator.o -o $@

$(OBJDIR)/test.o: $(SRCDIR)/test.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/test.cc

$(OBJDIR)/bw.o: $(SRCDIR)/bw.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/bw.cc

$(OBJDIR)/bw_basic_opts.o: $(SRCDIR)/bw_basic_opts.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/bw_basic_opts.cc

$(OBJDIR)/hmm.o: $(SRCDIR)/hmm.cc $(INCLUDEDIR)/hmm.h $(OBJDIR)/generator.o
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/hmm.cc

$(OBJDIR)/benchmark.o: $(SRCDIR)/benchmark.cc $(INCLUDEDIR)/benchmark.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/benchmark.cc

$(OBJDIR)/generator.o: $(SRCDIR)/generator.cc $(INCLUDEDIR)/generator.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/generator.cc

$(OBJDIR)/baum_baseline_impl.o: $(SRCDIR)/baum_baseline_impl.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/baum_baseline_impl.cc


debug: $(OUTDIR)/baum_welch_dbg $(OUTDIR)/baum_welch_basic_opts_dbg

$(OUTDIR)/baum_welch_dbg: $(OBJDIR)/bw_dbg.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o
	$(CC) $(OBJDIR)/bw_dbg.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o $(OBJDIR)/generator.o -o $@

$(OBJDIR)/bw_dbg.o: $(SRCDIR)/bw.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(DBGFLAGS) $(CFLAGS) $(SRCDIR)/bw.cc

$(OUTDIR)/baum_welch_basic_opts_dbg: $(OBJDIR)/bw_dbg_basic_opts.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o
	$(CC) $(OBJDIR)/bw_dbg_basic_opts.o $(OBJDIR)/hmm.o $(OBJDIR)/test.o $(OBJDIR)/benchmark.o $(OBJDIR)/generator.o -o $@

$(OBJDIR)/bw_dbg_basic_opts.o: $(SRCDIR)/bw_basic_opts.cc $(INCLUDEDIR)/bw.h
	$(CC) -o $@ $(DBGFLAGS) $(CFLAGS) $(SRCDIR)/bw_basic_opts.cc

clean:
	rm -f $(OBJDIR)/*.o $(OUTDIR)/*
