/*##############################################################################
# Name: clustGraphStructMinValues.h
# Use: Has stucter to store user input parameters
##############################################################################*/

#ifndef CLUSTERSTRUCTMINVALUES_H
#define CLUSTERSTRUCTMINVALUES_H

typedef struct minValues
{
    unsigned long minSharedEdgesUInt;
    unsigned long minReadsPerClustULng,
                  minReadsULng;
    double minMapqDbl;
}minValues;

void initMinValuesStruct(
    struct minValues *initStruct /*minValues struct to initalize*/
); /*sets all variables in minValue struct to default values*/

#endif
