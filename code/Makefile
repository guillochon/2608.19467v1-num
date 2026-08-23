CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -Wpedantic
LDLIBS ?= -lgmp

all: likelihood

likelihood: likelihood.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

check: likelihood
	@set -eu; \
	got="$$(./likelihood c5blowup 3 2>/dev/null)"; \
	expect='n=15 likelihood=63977511069907/43503039787261205233506826321920000000000'; \
	test "$$got" = "$$expect"; \
	got="$$(./likelihood bipartite 7 8 2>/dev/null)"; \
	expect='n=15 likelihood=7628328998218493/1044072954894268925604163831726080000000000'; \
	test "$$got" = "$$expect"; \
	echo 'published n=15 exact checks passed'

clean:
	rm -f likelihood

.PHONY: all check clean
