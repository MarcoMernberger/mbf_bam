#!/bin/bash
#a) you need to run anysnake shell --with-perf
#b) your cargo.toml must contain
#[profile.release]
#debug = true
#lto = false
#c) build your stuff
#cargo build --release
perf record --call-graph dwarf,16384 -e cpu-clock -F 997 python benchmark_read_counting.py
perf script | ./FlameGraph/stackcollapse-perf.pl | ./FlameGraph/stackcollapse-recursive.pl | c++filt | ./rust-unmangle/rust-unmangle | ./FlameGraph/flamegraph.pl > flame.svg
echo "result in flame.svg"


