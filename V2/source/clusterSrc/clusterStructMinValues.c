/*##############################################################################
# Name: clustGraphStructMinValues.h
# Use: Has intalization function for minValues stucter (stores user parameters)
##############################################################################*/

#include "clusterStructMinValues.h"

/*##############################################################################
# Output: Modifies: initStruct to have variables set to default values
##############################################################################*/
void initMinValuesStruct(
    struct minValues *initStruct /*minValues struct to initalize*/
) /*sets all variables in minValue struct to default values*/
{ /*initMinValuesStruct*/
    initStruct->minSharedEdgesUInt = 50; /*node shares 3 edges with cluster*/
    initStruct->minReadsPerClustULng = 100; /*min # reads to keep a cluster*/
    initStruct->minReadsULng = 100;        /*min @ of reads to a single read*/
    initStruct->minMapqDbl = 13;          /*min mapq to keep an alignemnt*/

    return;
} /*initMinValuesStruct*/
