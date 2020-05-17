### bw

**forward_backward:** `3(T-1)M(M+1) + 2M` muls, `T` divs, `(T-1)M(2*M + 1) + M` adds

**update and check:** `2M((T-1)*3*M + T)` muls, `M*(T*(M+1) + N)` divs, `M*((T-1)*(2M+1) + T(N+2))` adds


### bw basic opts

**forward_backward:** `3(T-1)M(M+1) + 2M` muls, `T` divs, `(T-1)M(2*M + 1) + M` adds

**update and check:** `M((T-1)*3*M + T)` muls, `M*(T + M + N)` divs, `M*((T-1)*(M+1) + T(N+1))` adds


### bw opts v2

**forward_backward:** `(T-1)M(3M+2) + 2M` muls, `(T-1)(M+1) + 1` divs, `(T-1)M(2*M + 1) + M` adds  (SAME COUNT AS BASIC_OPTS)

**update and check:** `M*((N+2)*T + M*(2T+1)` muls, `M*(N*(T+1) + M + 2T)` divs, `M*(T*(N+1) + (T-1)*(M+1))` adds)

#### Analysis of counts can be found within each corresponding source file
