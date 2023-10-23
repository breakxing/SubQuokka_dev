#!/bin/sh
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit/ori/25/cpu.ini -c ../circuit/ori/25/qaoa25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit/sub/25/cpu.ini -c ../circuit/sub/25/qaoa25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit/ideal/25/cpu.ini -c ../circuit/ideal/25/qaoa25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit_io/ori/25/cpu.ini -c ../circuit_io/ori/25/qaoa25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit_io/sub/25/cpu.ini -c ../circuit_io/sub/25/qaoa25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit_io/ideal/25/cpu.ini -c ../circuit_io/ideal/25/qaoa25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit/ori/25/cpu.ini -c ../circuit/ori/25/qft25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit/sub/25/cpu.ini -c ../circuit/sub/25/qft25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit/ideal/25/cpu.ini -c ../circuit/ideal/25/qft25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit_io/ori/25/cpu.ini -c ../circuit_io/ori/25/qft25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit_io/sub/25/cpu.ini -c ../circuit_io/sub/25/qft25.txt
perf stat -a -e cache-misses,cache-references ../../Quokka -i ../circuit_io/ideal/25/cpu.ini -c ../circuit_io/ideal/25/qft25.txt