#!/bin/bash
#
exe=../maureparticle.x
testList=( Euler-Ts30-Steady Euler-Ts30-DYN Euler-Ts1-Steady RK2-Ts30-Steady RK2-Ts30-DYN Euler-Ts1-DYN )
testList5=( Euler-Ts30-Steady-5 Euler-Ts30-DYN-5 Euler-Ts1-Steady-5 RK2-Ts30-Steady-5 RK2-Ts30-DYN Euler-Ts1-DYN-5 )
logFile="runTests.log"
#
# run test inputs and produce test outputs - all file names
# must be specified on the command line because the Maureparticle
# code will look for file name in all caps by default (which
# works correctly on case insensitive filesystems, i.e., Windows)
# then compare test outputs with expected outputs
passed=0
failed=0
echo "[$(date +'%Y-%h-%d-T%H:%M:%S%z')] start maureparticle $logFile" > $logFile
for t in ${testList[@]}; do
    $exe --meshfile fort.14                 \
         --maureparticleoutputfile $t.xyz   \
         --maureparameterinputfile $t.inp   \
         --elementlookuptablefile el2el.tbl \
         --nodelookuptablefile node2el.tbl  \
         --velocityfile fort.64             \
         > $t.log 2>&1
    echo "[$(date +'%Y-%h-%d-T%H:%M:%S%z')] $t" >> $logFile
    cat $t.log >> $logFile
    if [[ $(diff --report-identical-files $t.xyz __$t.xyz) == *"identical" ]]; then
        echo $t passed | tee --append $t.log
        passed=$(( $passed + 1 ))
    else
        echo $t failed | tee --append $t.log
        failed=$(( $failed + 1 ))
    fi
done
# test the use of a separate file with initial particle locations
# that automatically uses all particles found in the file without
# the number of particles being explicitly provided
for t in ${testList5[@]}; do
    $exe --meshfile fort.14                 \
         --maureparticleoutputfile $t.xyz   \
         --maureparameterinputfile "none"   \
         --maureparticleinputfile  initial-particle-locations.inp \
         --elementlookuptablefile el2el.tbl \
         --nodelookuptablefile node2el.tbl  \
         --velocityfile fort.64             \
         --slam0 -71.0                      \
         --sfea0 40.0                       \
         --initial-search-only              \
         --initiallocationoutputfile $t.lct \
         > $t.log 2>&1
    echo "[$(date +'%Y-%h-%d-T%H:%M:%S%z')] $t" >> $logFile
    cat $t.log >> $logFile
    if [[ $(diff --report-identical-files $t.xyz __$t.xyz) == *"identical" ]]; then
        echo $t passed | tee --append $t.log
        passed=$(( $passed + 1 ))
    else
        echo $t failed | tee --append $t.log
        failed=$(( $failed + 1 ))
    fi
done
echo $passed passed | tee --append $logFile
echo $failed failed | tee --append $logFile