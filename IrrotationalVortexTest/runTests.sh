#!/bin/bash
#
exe=../maureparticle.x
testList=(  Euler-Ts30-Steady   Euler-Ts30-DYN   Euler-Ts1-Steady   RK2-Ts30-Steady   RK2-Ts30-DYN   Euler-Ts1-DYN )
declare -a testList5
for t in ${testList[@]}; do
    testList5+=( ${t}-5 )
done
logFile="runTests.log"
#
# run test inputs and produce test outputs - all file names
# must be specified on the command line because the Maureparticle
# code will look for file name in all caps by default (which
# works correctly on case insensitive filesystems, i.e., Windows)
# then compare test outputs with expected outputs
tests=0
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
    tests=$(( $tests + 1 ))
done
# test the use of a separate file with initial particle locations
# that automatically uses all particles found in the file without
# the number of particles being explicitly provided
for t in ${testList5[@]}; do
    $exe --meshfile fort.14                 \
         --maureparticleoutputfile $t.xyz   \
         --maureparameterinputfile $t.inp   \
         --maureparticleinputfile  initial-particle-locations.inp \
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
    tests=$(( $tests + 1 ))
done
# test the use of a separate file with initial particle locations
# that automatically uses all particles found in the file without
# the number of particles being explicitly provided and the parameter
# file not being provided
t=${testList5[0]}
$exe --meshfile fort.14                 \
     --maureparameterinputfile "none"   \
     --maureparticleinputfile  initial-particle-locations.inp \
     --elementlookuptablefile el2el.tbl \
     --nodelookuptablefile node2el.tbl  \
     --slam0 -71.0                      \
     --sfea0 40.0                       \
     --initial-search-only              \
     --initiallocationoutputfile $t.lct \
     > ${t}-lct.log 2>&1
echo "[$(date +'%Y-%h-%d-T%H:%M:%S%z')] $t" >> $logFile
cat ${t}-lct.log >> $logFile
if [[ $(diff --report-identical-files $t.lct __$t.lct) == *"identical" ]]; then
    echo ${t}-lct passed | tee --append ${t}-lct.log
    passed=$(( $passed + 1 ))
else
    echo ${t}-lct failed | tee --append ${t}-lct.log
    failed=$(( $failed + 1 ))
fi
tests=$(( $tests + 1 ))
# report number of tests, how many passed and failed
echo $tests  tests  | tee --append $logFile
echo $passed passed | tee --append $logFile
echo $failed failed | tee --append $logFile