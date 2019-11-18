// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2019  MDOODZ Developper team
//
// This file is part of MDOODZ.
//
// MDOODZ is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MDOODZ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MDOODZ.  If not, see <http://www.gnu.org/licenses/>.
// =========================================================================

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "ctype.h"
#include "header_MDOODZ.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MinMaxArrayI( int* array, int scale, int size, char* text ) {
    // Search min/max values of an array and prints it to standard out
    int min=array[0], max=array[0];
    int k;
    for (k=1; k<size; k++) {
        if (array[k]>max) max = array[k];
        if (array[k]<min) min = array[k];
    }
    printf( "min(%s) = %d max(%s) = %d\n", text, min*scale, text, max*scale);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  ArrayPlusArray( DoodzFP* arr1, DoodzFP* arr2, int size ) {
    int k;
#pragma omp parallel for shared( arr1, arr2) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] += arr2[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  Initialise1DArrayChar( char* arr, int n, char val ) {
    int k;
#pragma omp parallel for shared( arr ) private( k ) firstprivate( val, n ) schedule( static )
    for(k=0;k<n;k++) {
        arr[k] = val;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double SumArray( double* array, double scale, int size, char* text ) {
    // Search min/max values of an array and prints it to standard out
    double sum=0.0;
    int k;
    for (k=0; k<size; k++) {
        sum += (array[k]);
    }
    printf( "sum(%s) = %2.12e\n", text, sum*scale );
    return sum;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateInputFile( char fin_name[], int NewNumber ) {
    
    FILE* fin;
    char  FieldName[] = "istep";
    int h, find=0, bufmax=50, length;
    char line[bufmax];
    
    // Get the length of the demanded string
    length = strlen( FieldName );
    
    // Buffer array to contain the string to compare with
    char *param1, *param2;
    param2 = malloc( (length+1)*sizeof(char));
    asprintf(&param1, "%s", FieldName);
    
    // Open File
    if (fopen(fin_name,"rwt") == NULL) {
        printf("Setup file '%s' does not exist\nExiting...\n", fin_name);
        exit(1);
    }
    else {
        fin = fopen(fin_name,"r+");
    }
    // Seach for sign 'equal' in the lines
    while ( find == 0 && feof(fin) !=1 ) {
        
        // Read new line
        fgets ( line, sizeof(line), fin );
        if (feof(fin)) {
            printf("Warning : Parameter '%s' not found in the setup file\n", FieldName);
        }
        
        // Get the first 'length' characters of the line
        for (h=0;h<length;h++) {
            param2[h] = line[h];
        }
        param2[length] = '\0';
        
        // printf("p1 %s p2 %s\n", param1, param2);
        
        // Check if we found the right parameter name
        if ( strcmp(param1, param2) == 0 ) {
            find = 1;
            // Search for equal sign and replace by the current step
            for (h=0;h<bufmax;h++) {
                
                //printf(" %c ", line[h]);
                
                if(strlen(line)> 0 && line[h]=='=') {
                    
                    //printf("pos =%d SEEK_CUR=%d NewNumber=%d\n", pos, SEEK_CUR, NewNumber);
                    fseek(fin, -6, SEEK_CUR);
                    fprintf( fin, "%05d", NewNumber);
                    break;
                }
            }
        }
    }
    
    // Close file
    fclose(fin);
    free(param1);
    free(param2);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  ArrayMinusArray( DoodzFP* arr1, DoodzFP* arr2, int size ) {
    int k;
#pragma omp parallel for shared( arr1, arr2) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] -= arr2[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  Initialise1DArrayInt( int* arr, int n, int val ) {
    int k;
#pragma omp parallel for shared( arr ) private( k ) firstprivate( val, n ) schedule( static )
    for(k=0;k<n;k++) {
        arr[k] = val;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  Initialise1DArrayDouble( double* arr, int n, double val ) {
    int k;
#pragma omp parallel for shared( arr ) private( k ) firstprivate( val, n ) schedule( static )
    for(k=0;k<n;k++) {
        arr[k] = val;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void ArrayTimesScalar( double* arr1, double scalar, int size ) {
    int k;
#pragma omp parallel for shared( arr1, scalar) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] *= scalar;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MinMaxArray( DoodzFP* array, double scale, int size, char* text ) {
    
    // Search min/max values of an array and prints it to standard out
    double min=array[0], max=array[0];
    int k;
    
#pragma omp parallel
    {
        double pmin=array[0], pmax=array[0];
#pragma omp for
        for (k=0; k<size; k++) {
            if (array[k]>pmax) pmax = array[k];
            if (array[k]<pmin) pmin = array[k];
        }
#pragma omp flush (max)
        if (pmax>max) {
#pragma omp critical
            {
                if (pmax>max) max = pmax;
            }
        }
        
#pragma omp flush (min)
        if (pmin<min) {
#pragma omp critical
            {
                if (pmin<min) min = pmin;
            }
        }
    }
    printf( "min(%s) = %2.6e max(%s) = %2.6e\n", text, min*scale, text, max*scale);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MinMaxArrayTag( DoodzFP* array, double scale, int size, char* text, char* tag ) {
    
    // Search min/max values of an array and prints it to standard out
    double min=array[0], max=array[0];
    int k;
    
#pragma omp parallel
    {
        double pmin=array[0], pmax=array[0];
#pragma omp for
        for (k=0; k<size; k++) {
            if (tag[k] < 30) {
                if (array[k]>pmax) pmax = array[k];
                if (array[k]<pmin) pmin = array[k];
            }
        }
#pragma omp flush (max)
        if (pmax>max) {
#pragma omp critical
            {
                if (pmax>max) max = pmax;
            }
        }
        
#pragma omp flush (min)
        if (pmin<min) {
#pragma omp critical
            {
                if (pmin<min) min = pmin;
            }
        }
    }
    printf( "min(%s) = %2.6e max(%s) = %2.6e\n", text, min*scale, text, max*scale);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MinMaxArrayVal( DoodzFP* array, int size, double* min, double* max ) {
    
    // Search min/max values of an array and prints it to standard out
    *min=array[0];
    *max=array[0];
    int k;
    
#pragma omp parallel
    {
        double pmin=array[0], pmax=array[0];
#pragma omp for
        for (k=0; k<size; k++) {
            if (array[k]>pmax) pmax = array[k];
            if (array[k]<pmin) pmin = array[k];
        }
#pragma omp flush (max)
        if (pmax>*max) {
#pragma omp critical
            {
                if (pmax>*max) *max = pmax;
            }
        }
        
#pragma omp flush (min)
        if (pmin<*min) {
#pragma omp critical
            {
                if (pmin<*min) *min = pmin;
            }
        }
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  ArrayEqualArray( DoodzFP* arr1, DoodzFP* arr2, int size ) {
    int k;
#pragma omp parallel for shared( arr1, arr2) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] = arr2[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  ArrayEqualArrayI( int* arr1, int* arr2, int size ) {
    int k;
#pragma omp parallel for shared( arr1, arr2) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] = arr2[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ArrayTimesScalarArray( double* arr1, double scalar, double* arr2, int size ) {
    int k;
#pragma omp parallel for shared( arr1, arr2, scalar) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] *= scalar*arr2[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ArrayDividedScalarArray( double* arr1, double scalar, double* arr2, int size ) {
    int k;
#pragma omp parallel for shared( arr1, arr2, scalar) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] /= scalar*arr2[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ArrayPlusScalarArray( double* arr1, double scalar, double* arr2, int size ) {
    int k;
#pragma omp parallel for shared( arr1, arr2, scalar) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] += scalar*arr2[k];
        //        printf("%2.2e %2.2e %2.2e\n", scalar*arr2[k], scalar, arr2[k]);
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  IsNanArray2DFP( DoodzFP* arr1, int size1 ) {
    int k;
    //#pragma omp parallel for shared( arr1, mesh ) private(k) schedule( static )
    for(k=0;k<size1;k++) {
        
        if ( isnan(arr1[k])!=0 ) {
            printf("Nan Scheisse!\n");
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  IsInfArray2DFP( DoodzFP* arr1, int size1 ) {
    int k;
    //#pragma omp parallel for shared( arr1, mesh ) private(k) schedule( static )
    for(k=0;k<size1;k++) {
        
        if ( isinf(arr1[k])!=0 ) {
            printf("Inf Scheisse!\n");
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  Initialise2DArray( DoodzFP* arr, int nx, int nz, double val ) {
    int c;
#pragma omp parallel for shared( arr ) private( c ) firstprivate( val, nx, nz ) schedule( static )
    for(c=0;c<nx*nz;c++) {
        arr[c] = val;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void  Initialise2DArrayInt( int* arr, int nx, int nz, int val ) {
    int c;
#pragma omp parallel for shared( arr ) private( c ) firstprivate( val, nx, nz ) schedule( static )
    for(c=0;c<nx*nz;c++) {
        arr[c] = val;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MinMaxArrayPart( DoodzFP* array, double scale, int size, char* text, int* phase ) {
    
    // Search min/max values of an array and prints it to standard out
    double min=array[0], max=array[0];
    int k;
    
#pragma omp parallel
    {
        double pmin=array[0], pmax=array[0];
#pragma omp for
        for (k=0; k<size; k++) {
            if (array[k]>pmax && phase[k]!=-1 ) pmax = array[k];
            if (array[k]<pmin && phase[k]!=-1 ) pmin = array[k];
        }
#pragma omp flush (max)
        if (pmax>max) {
#pragma omp critical
            {
                if (pmax>max) max = pmax;
            }
        }
        
#pragma omp flush (min)
        if (pmin<min) {
#pragma omp critical
            {
                if (pmin<min) min = pmin;
            }
        }
    }
    printf( "min(%s) = %2.12e max(%s) = %2.12e\n", text, min*scale, text, max*scale);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

