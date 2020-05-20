#include <cstdint>

#define VOLATILE __volatile__
#define ASM __asm__


/* RDTSCP is an instruction on several Intel and compatible CPUs that Reads the
 * Time Stamp Counter. The Intel manuals contain more information. */
#define COUNTER_LO(a) ((a).int32.lo)
#define COUNTER_HI(a) ((a).int32.hi)
#define COUNTER_VAL(a) ((a).int64)


#define COUNTER(a) \
	((uint64_t) COUNTER_VAL(a))

#define COUNTER_DIFF(a,b) \
	(COUNTER(a)-COUNTER(b))


typedef union
{       uint64_t int64;
		struct {uint32_t lo, hi;} int32;
} tsc_counter;


#define RDTSC2(cpu_c) \
	  ASM VOLATILE ("rdtscp" : "=a" ((cpu_c).int32.lo), "=d"((cpu_c).int32.hi)) \

#define CPUID() \
		ASM VOLATILE ("cpuid" : : "a" (0) : "bx", "cx", "dx" )

#define LFENCE() \
		asm __volatile__ ("lfence \n": )


static void init_tsc() {
	; // no need to initialize anything for x86
}

static uint64_t start_tsc(void) {
    tsc_counter start;
	CPUID();
    RDTSC2(start);
    return COUNTER_VAL(start);
}

static uint64_t stop_tsc() {
	tsc_counter end;
	RDTSC2(end);
	CPUID();
	return COUNTER_VAL(end);
}
