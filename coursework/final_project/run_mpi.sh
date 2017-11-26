#!/bin/bash

for WORKFLOW in 13 99
do
    if [ "${WORKFLOW}" -eq 13 ]; then
        WORKFLOW_NAME="CyberShake"
    else
        WORKFLOW_NAME="Montage"
    fi
    for WORKFLOW_SIZE in 50 100 200
    do
        FILENAME="${WORKFLOW_NAME}_${WORKFLOW_SIZE}.xml.dag"
        if [ "${WORKFLOW}" -eq "99" ]; then
            if [ "${WORKFLOW_SIZE}" -eq 50 ]; then
                BEST=171
            elif [ "${WORKFLOW_SIZE}" -eq 100 ]; then
                BEST=362
            else
                BEST=781
                FILENAME="MONTAGE.n.200.0.dax.dag"
            fi
        else
            if [ "${WORKFLOW_SIZE}" -eq 50 ]; then
                BEST=799
            elif [ "${WORKFLOW_SIZE}" -eq 100 ]; then
                BEST=1773
            else
                BEST=2326
                FILENAME="CYBERSHAKE.n.200.19.dax.dag"
            fi
        fi
        for POPULATION in 25 50 100
        do
            PRINT_TEXT="${WORKFLOW_NAME}_${WORKFLOW_SIZE} | ${POPULATION}"
            echo "${PRINT_TEXT} BEGIN"
            for TEST_NUMBER in 1 2 3 4 5
            do
                echo -n "${PRINT_TEXT} | $TEST_NUMBER/5: "
                /opt/openmpi-2.1.1/bin/mpirun -n 4 --machinefile machines bin/HEA_mpi.o -c cluster/cluster.vcl -w "${FILENAME}" -p "${POPULATION}" -g 100 -b "${BEST}"
            done
            echo "${PRINT_TEXT} END"
            date
            echo ""
        done
    done
done
