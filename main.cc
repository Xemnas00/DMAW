/**
    MAW: Minimal Absent Words
    Copyright (C) 2014 Alice Heliou and Solon P. Pissis.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "mawdefs.h"

void save_on_file_simple_matrix_double(char * filename, double ** matrix, INT num_seqs, unsigned char ** seqs_id);
void save_on_file_INT(char * filename, INT ** matrix, INT num_seqs, unsigned char ** seqs_id);
void save_on_file_double(char * filename, double ** matrix, INT num_seqs, unsigned char ** seqs_id);
void save_on_file_INT_histograms(char * filename, INT ** hist_matrix, INT num_seqs, INT k, INT K, unsigned char ** seqs_id);
void save_on_file_INT_histograms_complex(char * filename, INT ** matrix, INT num_seqs, INT k, INT K, unsigned char ** seqs_id);


int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
	FILE *          out_fd;                 // the output file descriptor
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name
        unsigned char * seq    = NULL;         // the sequence in memory
        unsigned char * seq_id = NULL;         // the sequence id in memory
        unsigned char ** seqs    = NULL;         // the sequences in memory
        unsigned char ** seqs_id = NULL;         // the sequences id in memory
        unsigned int    num_seqs = 0;           // the total number of sequences considered
	char *          alphabet;               // the alphabet
	unsigned int    i, j;

	double ** D, ** SC;
	double * buf, * bufSC;

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );


	/* Check the arguments */
        if ( i < 9 )
        {
                usage ();
                return ( 1 );
        }
        else
        {
                if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet = ( char * ) DNA;
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet = ( char * ) PROT;
                else if ( ! strcmp ( "SIXTEEN", sw . alphabet ) ) alphabet = ( char * ) SIXTEEN;
                else
                {
                        for(i = 0; i < num_seqs; i++) fprintf( out_fd, "\t%s", seqs_id[i] );

    fprintf( out_fd, "\n");

	for ( i = 0; i < num_seqs; i++ )
	{
		fprintf( out_fd, "%s", seqs_id[i] );
		for ( j = 0; j < num_seqs; j++ )
		{
			fprintf( out_fd, "\t%.2lf", D[i][j]);
		}
		fprintf( out_fd, "\n");
	}

	fclose ( out_fd );    fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
                        return ( 1 );
                }

                input_filename          = sw . input_filename;
                output_filename         = sw . output_filename;
        }


	double start = gettime();

	/* Read the (Multi)FASTA file in memory */
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

	char c;
	c = fgetc( in_fd );
	do
	{
		unsigned int seq_id_len = 0;
		unsigned int max_alloc_seq_id = 0;
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
			return ( 1 );
		}
		else
		{
			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id )
				{
					seq_id = ( unsigned char * ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id += ALLOC_SIZE;
				}
				if( c == ' ' ) continue;
				else	seq_id[ seq_id_len++ ] = c;
			}
			seq_id[ seq_id_len ] = '\0';

		}

		INT max_alloc_seq = 0;
		INT seq_len = 0;

		while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
		{
			if( seq_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
				c = fgetc( in_fd );
				break;
			}
			if( c == '\n' ) continue;

			c = toupper( c );

			if ( seq_len >= max_alloc_seq )
			{
				seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
			{
				seq[ seq_len++ ] = c;
			}
			else
			{
				if ( strchr ( IUPAC, c ) )
				{
					seq[ seq_len++ ] = 'N';
				}
				else if ( c == ' ' )
				{
				}
				else
				{
					fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
					return ( 1 );
				}
			}

		}
		if( seq_len != 0 )
		{
			num_seqs++;
			if ( seq_len >= max_alloc_seq )
			{
				seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq += ALLOC_SIZE;
			}
			seq[ seq_len ] = '\0';

			if ( sw . c )
			{
				if ( sw . K > seq_len ) sw . K = seq_len;

				seq = ( unsigned char * ) realloc ( seq,   ( seq_len + seq_len + 1 ) * sizeof ( unsigned char ) );
				unsigned char * seq2 = ( unsigned char * ) strdup ( ( const char * ) seq );
				strcat ( ( char * ) seq, ( char * ) seq2 );
				seq_len += seq_len;
				free ( seq2 );

				#ifdef _USE_32
				fprintf( stderr, "Processing circular sequence %s of length %d.\n", ( char * ) seq_id, seq_len );
				#endif

				#ifdef _USE_64
				fprintf( stderr, "Processing circular sequence %s of length %ld.\n", ( char * ) seq_id, seq_len );
				#endif
			}
			else
			{
				#ifdef _USE_32
				fprintf( stderr, "Processing linear sequence %s of length %d.\n", ( char * ) seq_id, seq_len );
				#endif

				#ifdef _USE_64
				fprintf( stderr, "Processing linear sequence %s of length %ld.\n", ( char * ) seq_id, seq_len );
				#endif
			}
			seq[ seq_len ] = '\0';

			seqs = ( unsigned char ** ) realloc ( seqs,   ( num_seqs ) * sizeof ( unsigned char * ) );
			seqs[num_seqs - 1] = ( unsigned char * ) calloc ( seq_len + 1,   sizeof ( unsigned char ) );
			strcat ( ( char * ) seqs[ num_seqs - 1 ], ( char * ) seq );
			seqs[num_seqs - 1][ seq_len ] = '\0';
			seqs_id = ( unsigned char ** ) realloc ( seqs_id,   ( num_seqs ) * sizeof ( unsigned char * ) );
			seqs_id[num_seqs - 1] = ( unsigned char * ) calloc ( seq_id_len + 1,   sizeof ( unsigned char ) );
			strcat ( ( char * ) seqs_id[ num_seqs - 1 ], ( char * ) seq_id );
			seqs_id[num_seqs - 1][ seq_id_len ] = '\0';

		}
		free ( seq );
		seq = NULL;
		free ( seq_id );
		seq_id = NULL;

	} while( c != EOF );

#if 0
	for ( i = 0; i < NmawX; i ++ )
	  fprintf( stderr, "<%c, %d, %d>\n", ( char ) mawX[i] . letter, mawX[i] . pos, mawX[i] . size );

	fprintf( stderr, "-----------\n" );

	for ( i = 0; i < NmawY; i ++ )
	  fprintf( stderr, "<%c, %d, %d>\n", ( char ) mawY[i] . letter, mawY[i] . pos, mawY[i] . size );
#endif
	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	if ( num_seqs == 1 )
	{
		fprintf ( stderr, " Error: the number of input sequences must be greater than 1!\n" );
		return ( 1 );
	}


	/* 2d dynamic memory allocation of the SCMAW Distance Matrix */
	SC = ( double ** ) malloc ( ( num_seqs ) * sizeof ( double * ) );
	if ( SC == NULL )
	{
		fprintf ( stderr, " ERROR: SC distance matrix could not be allocated!!!\n" );
		return ( 1 );
	}
	bufSC =  ( double * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( double ) );
	if ( bufSC == NULL )
	{
		fprintf ( stderr, " ERROR: SC distance could not be allocated!!!\n" );
		return ( 1 );
	}
	for ( i = 0; i < num_seqs; ++ i ) 	SC[i] = &bufSC[( size_t ) i * ( size_t ) ( num_seqs ) ];

	/* 2d dynamic memory allocation of the DMAW Distance Matrix */
	D = ( double ** ) malloc ( ( num_seqs ) * sizeof ( double * ) );
   	if ( D == NULL )
    	{
      		fprintf ( stderr, " ERROR: D distance matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
	buf =  ( double * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( double ) );
   	if ( buf == NULL )
    	{
      		fprintf ( stderr, " ERROR: D distance matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
        for ( i = 0; i < num_seqs; ++ i ) 	D[i] = &buf[( size_t ) i * ( size_t ) ( num_seqs ) ];

	fprintf( stderr, "Computing minimal absent words and making the comparison.\n" );


	/*Matrices*/
	//SCMAW set highest occurences matrix
	INT ** scMawHighestOccMat;
	INT * bufScMawHighestOccMat;
	scMawHighestOccMat = (INT **) malloc ( (num_seqs) * sizeof(INT *) );

	bufScMawHighestOccMat =  ( INT * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( INT ) );
   	if ( bufScMawHighestOccMat == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs; ++ i ) 	scMawHighestOccMat[i] = &bufScMawHighestOccMat[( size_t ) i * ( size_t ) ( num_seqs ) ];

    //DMAW set highest occurences Matrix
    INT ** dMawHighestOccMat;
	INT * bufDMawHighestOccMat;
	dMawHighestOccMat = (INT **) malloc ( (num_seqs) * sizeof(INT *) );

	bufDMawHighestOccMat =  ( INT * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( INT ) );
   	if ( bufDMawHighestOccMat == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs; ++ i ) 	dMawHighestOccMat[i] = &bufDMawHighestOccMat[( size_t ) i * ( size_t ) ( num_seqs ) ];

    //UNION set highest occurences matrix
    INT ** UnionHighestOccMat;
	INT * UnionbufHighestOccMat;
	UnionHighestOccMat = (INT **) malloc ( (num_seqs) * sizeof(INT *) );

	UnionbufHighestOccMat =  ( INT * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( INT ) );
   	if ( UnionbufHighestOccMat == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs; ++ i ) 	UnionHighestOccMat[i] = &UnionbufHighestOccMat[( size_t ) i * ( size_t ) ( num_seqs ) ];

    //DMAW/SCMAW cardinality ratio matrix
    double ** cardRatioMatrix;
	double * bufCardRatioMatrix;
	cardRatioMatrix = (double **) malloc ( (num_seqs) * sizeof(double *) );

	bufCardRatioMatrix =  ( double * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( double ) );
   	if ( bufCardRatioMatrix == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs; ++ i ) 	cardRatioMatrix[i] = &bufCardRatioMatrix[( size_t ) i * ( size_t ) ( num_seqs ) ];

    //DMAW/Union cardinality ratio matrix
    double ** UnionCardRatioMatrix;
	double * UnionbufCardRatioMatrix;
	UnionCardRatioMatrix = (double **) malloc ( (num_seqs) * sizeof(double *) );

	UnionbufCardRatioMatrix =  ( double * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( double ) );
   	if ( UnionbufCardRatioMatrix == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs; ++ i ) 	UnionCardRatioMatrix[i] = &UnionbufCardRatioMatrix[( size_t ) i * ( size_t ) ( num_seqs ) ];


    //DMAW/SCMAW total length ratio matrix
    double ** totalLengthRatioMatrix;
	double * bufTotalLengthMatrix;
	totalLengthRatioMatrix = (double **) malloc ( (num_seqs) * sizeof(double *) );

	bufTotalLengthMatrix =  ( double * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( double ) );
   	if ( bufTotalLengthMatrix == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs; ++ i ) 	totalLengthRatioMatrix[i] = &bufTotalLengthMatrix[( size_t ) i * ( size_t ) ( num_seqs ) ];

    //DMAW/Union total length ratio matrix
    double ** UnionTotalLengthRatioMatrix;
	double * UnionBufTotalLengthMatrix;
	UnionTotalLengthRatioMatrix = (double **) malloc ( (num_seqs) * sizeof(double *) );

	UnionBufTotalLengthMatrix =  ( double * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( double ) );
   	if ( UnionBufTotalLengthMatrix == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs; ++ i ) 	UnionTotalLengthRatioMatrix[i] = &UnionBufTotalLengthMatrix[( size_t ) i * ( size_t ) ( num_seqs ) ];

    //Total MAWs histogram matrix
    INT ** histMatrix;
	INT * bufHistMatrix;
	histMatrix = (INT **) malloc ( (num_seqs) * sizeof(INT *) );

	bufHistMatrix =  ( INT * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( sw . K - sw . k + 1 ) ), sizeof ( INT ) );
   	if ( bufHistMatrix == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs; ++ i ) 	histMatrix[i] = &bufHistMatrix[( size_t ) i * ( size_t ) ( sw . K - sw . k + 1 ) ];

    //DMAW MAWs histogram matrix
    INT ** DMAWhistMatrix;
	INT * DMAWbufHistMatrix;
	DMAWhistMatrix = (INT **) malloc ( (num_seqs * num_seqs ) * sizeof(INT *) );

	DMAWbufHistMatrix =  ( INT * ) calloc ( ( ( size_t ) ( num_seqs * num_seqs ) ) * ( ( size_t ) ( sw . K - sw . k + 1 ) ), sizeof ( INT ) );
   	if ( DMAWbufHistMatrix == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs * num_seqs; ++ i ) 	DMAWhistMatrix[i] = &DMAWbufHistMatrix[( size_t ) i * ( size_t ) ( sw . K - sw . k + 1 ) ];

    //SCMAW MAWs histogram matrix
    INT ** SCMAWhistMatrix;
	INT * SCMAWbufHistMatrix;
	SCMAWhistMatrix = (INT **) malloc ( (num_seqs * num_seqs) * sizeof(INT *) );

	SCMAWbufHistMatrix =  ( INT * ) calloc ( ( ( size_t ) ( num_seqs * num_seqs ) ) * ( ( size_t ) ( sw . K - sw . k + 1 ) ), sizeof ( INT ) );
   	if ( SCMAWbufHistMatrix == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs * num_seqs; ++ i ) 	SCMAWhistMatrix[i] = &SCMAWbufHistMatrix[( size_t ) i * ( size_t ) ( sw . K - sw . k + 1 ) ];

    //UNION MAWs histogram matrix
    INT ** UnionhistMatrix;
	INT * UnionbufHistMatrix;
	UnionhistMatrix = (INT **) malloc ( (num_seqs * num_seqs) * sizeof(INT *) );

	UnionbufHistMatrix =  ( INT * ) calloc ( ( ( size_t ) ( num_seqs * num_seqs ) ) * ( ( size_t ) ( sw . K - sw . k + 1 ) ), sizeof ( INT ) );
   	if ( UnionbufHistMatrix == NULL )
    	{
      		fprintf ( stderr, " ERROR: scOcc matrix could not be allocated!!!\n" );
      		return ( 1 );
    	}
    for ( i = 0; i < num_seqs * num_seqs; ++ i ) 	UnionhistMatrix[i] = &UnionbufHistMatrix[( size_t ) i * ( size_t ) ( sw . K - sw . k + 1 ) ];

    FILE *mawFile;
    if ( ! ( mawFile = fopen ( "./RESULTS/MAWs.txt", "w") ) ) {
		fprintf ( stderr, " Error: Cannot open file %s!\n", "MAWs.txt" );
		return ( 1 );
	}

	#pragma omp parallel for
	for ( int i = 0; i < num_seqs; i++ )
	{
       
		TMaw * mawX = NULL;
		unsigned int NmawX = 0;
      
		compute_maw ( seqs[i], seqs_id[i], sw, &mawX, &NmawX );

		//printf("%ld", NmawX);
        fprintf(mawFile, "MAWs of %s=\"%s\":\n", seqs_id[i], seqs[i]);
		for(int maw_num = 0; maw_num < NmawX; maw_num++) {
            //printf("Lunghezza %ld rilevata!", mawX[maw_num].size);
            //if(mawX[maw_num].size - sw.k < 0) printf("NOPE");
            histMatrix[i][mawX[maw_num].size + 1 - sw.k]++;

            fprintf(mawFile, " <%c,%ld,%ld> = \"%c", mawX[maw_num] . letter, mawX[maw_num] . pos, mawX[maw_num] . size, mawX[maw_num].letter );
            for(INT k = 0; k < mawX[maw_num].size; k++) fprintf(mawFile, "%c", seqs[i][mawX[maw_num].pos + k]);
            fprintf(mawFile, "\"\n");
		}


		for ( int j = i; j < num_seqs; j++ )
		{
			TMaw * mawY = NULL;
			unsigned int NmawY = 0;

			compute_maw ( seqs[j], seqs_id[j], sw, &mawY, &NmawY );

			double XY_D = 0, XY_SC = 0;
			if ( i == j )	{
                D[i][j] = 0.0;
                scMawHighestOccMat[i][j] = 0;
                dMawHighestOccMat[i][j] = 0;
                cardRatioMatrix[i][j] = 0.0;
            }
			else
			{
				INT scMawHighestOcc = 0, dMawHighestOcc = 0, unionHighestOcc = 0, scNumber = 0, dNumber = 0, unionNumber = 0, scTotalLength = 0, dTotalLength = 0, unionTotalLength = 0;
 
				maw_seq_comp( seqs[i], mawX,  &NmawX, seqs[j], mawY, &NmawY, &XY_SC, sw . k, sw . K, &scMawHighestOcc, &scNumber, &scTotalLength, i, j, SCMAWhistMatrix, num_seqs );
				maw_seq_comp_factors( seqs[i], mawX, &NmawX, seqs[j], mawY, &NmawY, &XY_D, sw . k, sw . K, &dMawHighestOcc, &dNumber, &dTotalLength, i, j, DMAWhistMatrix, num_seqs, &unionHighestOcc, &unionNumber, &unionTotalLength, UnionhistMatrix );
				SC[i][j] = SC[i][j] = XY_SC;
				D[i][j] = D[j][i] = XY_D;
				scMawHighestOccMat[i][j] = scMawHighestOccMat[j][i] = scMawHighestOcc;
				dMawHighestOccMat[i][j] = dMawHighestOccMat[j][i] = dMawHighestOcc;
				UnionHighestOccMat[i][j] = UnionHighestOccMat[j][i] = unionHighestOcc;
				cardRatioMatrix[i][j] = cardRatioMatrix[j][i] = (double) dNumber / scNumber;
				UnionCardRatioMatrix[i][j] = UnionCardRatioMatrix[j][i] = (double) dNumber / unionNumber;
				totalLengthRatioMatrix[i][j] = totalLengthRatioMatrix[j][i] = (double) dTotalLength / scTotalLength;
				UnionTotalLengthRatioMatrix[i][j] = UnionTotalLengthRatioMatrix[j][i] = (double) dTotalLength / unionTotalLength;
			}
			free ( mawY );
		}
		free ( mawX );
	}

    fclose(mawFile);

	double end = gettime();

//    char * output_dist_filename = (char *) malloc(100 * sizeof(char));
//	sprintf(output_dist_filename, "./RESULTS/%s", output_filename);
//
//	if ( ! ( out_fd = fopen ( output_dist_filename, "w" ) ) )
//	{
//		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
//		return ( 1 );
//	}

	//fprintf( out_fd, "%d\n", num_seqs );
    //fprintf( out_fd, "\t");

//    for(i = 0; i < num_seqs; i++) fprintf( out_fd, "\t%s", seqs_id[i] );
//
//    fprintf( out_fd, "\n");
//
//	for ( i = 0; i < num_seqs; i++ )
//	{
//		fprintf( out_fd, "%s", seqs_id[i] );
//		for ( j = 0; j < num_seqs; j++ )
//		{
//			fprintf( out_fd, "\t%.2lf", D[i][j]);
//		}
//		fprintf( out_fd, "\n");
//	}
//
//	fclose ( out_fd );


	char * dmaw_output_dist_filename = (char *) malloc(100 * sizeof(char));
	sprintf(dmaw_output_dist_filename, "./RESULTS/DMAW_%s.csv", output_filename);
    save_on_file_simple_matrix_double(dmaw_output_dist_filename, D, num_seqs, seqs_id);

	char * scmaw_output_dist_filename = (char *) malloc(100 * sizeof(char));
	sprintf(scmaw_output_dist_filename, "./RESULTS/SCMAW_%s.csv", output_filename);
	save_on_file_simple_matrix_double(scmaw_output_dist_filename, SC, num_seqs, seqs_id);

	char * scmaw_most_common_lengths = (char *) malloc(100 * sizeof(char));
	sprintf(scmaw_most_common_lengths, "./RESULTS/SCMAW_most_common_lengths_%s.csv", output_filename);
	save_on_file_INT(scmaw_most_common_lengths, scMawHighestOccMat, num_seqs, seqs_id);

	char * dmaw_most_common_lengths = (char *) malloc(100 * sizeof(char));
	sprintf(dmaw_most_common_lengths, "./RESULTS/DMAW_most_common_lengths_%s.csv", output_filename);
	save_on_file_INT(dmaw_most_common_lengths, dMawHighestOccMat, num_seqs, seqs_id);

	char * union_most_common_lengths = (char *) malloc(100 * sizeof(char));
	sprintf(union_most_common_lengths, "./RESULTS/UNION_most_common_lengths_%s.csv", output_filename);
	save_on_file_INT(union_most_common_lengths, UnionHighestOccMat, num_seqs, seqs_id);

	char * card_ratios = (char *) malloc(100 * sizeof(char));
	sprintf(card_ratios, "./RESULTS/DMAW_symmetricDifference_cardinality_ratios_%s.csv", output_filename);
	save_on_file_double(card_ratios, cardRatioMatrix, num_seqs, seqs_id);

	char * union_card_ratios = (char *) malloc(100 * sizeof(char));
	sprintf(union_card_ratios, "./RESULTS/DMAW_UNION_cardinality_ratios_%s.csv", output_filename);
	save_on_file_double(union_card_ratios, UnionCardRatioMatrix, num_seqs, seqs_id);

	char * length_ratios = (char *) malloc(100 * sizeof(char));
	sprintf(length_ratios, "./RESULTS/DMAW_symmetricDifference_length_ratios_%s.csv", output_filename);
	save_on_file_double(length_ratios, totalLengthRatioMatrix, num_seqs, seqs_id);

	char * union_length_ratios = (char *) malloc(100 * sizeof(char));
	sprintf(union_length_ratios, "./RESULTS/DMAW_UNION_length_ratios_%s.csv", output_filename);
	save_on_file_double(union_length_ratios, UnionTotalLengthRatioMatrix, num_seqs, seqs_id);

	char * hist_seqs = (char *) malloc(100 * sizeof(char));
	sprintf(hist_seqs, "./RESULTS/Total_maws_histograms_%s.csv", output_filename);
	save_on_file_INT_histograms(hist_seqs, histMatrix, num_seqs, sw.k, sw.K, seqs_id);

	char * SCMAW_hist_seqs = (char *) malloc(100 * sizeof(char));
	sprintf(SCMAW_hist_seqs, "./RESULTS/SCMAW_maws_histograms_%s.csv", output_filename);
	save_on_file_INT_histograms_complex(SCMAW_hist_seqs, SCMAWhistMatrix, num_seqs, sw.k, sw.K, seqs_id);

	char * DMAW_hist_seqs = (char *) malloc(100 * sizeof(char));
	sprintf(DMAW_hist_seqs, "./RESULTS/DMAW_maws_histograms_%s.csv", output_filename);
	save_on_file_INT_histograms_complex(DMAW_hist_seqs, DMAWhistMatrix, num_seqs, sw.k, sw.K, seqs_id);

	char * UNION_hist_seqs = (char *) malloc(100 * sizeof(char));
	sprintf(UNION_hist_seqs, "./RESULTS/UNION_maws_histograms_%s.csv", output_filename);
	save_on_file_INT_histograms_complex(UNION_hist_seqs, UnionhistMatrix, num_seqs, sw.k, sw.K, seqs_id);


        //fprintf( stderr, "Elapsed time for processing %d sequence(s): %lf secs.\n", num_seqs, ( end - start ) );

	for ( i = 0; i < num_seqs; i ++ )
		free ( seqs[i] );
	free ( seqs);
	for ( i = 0; i < num_seqs; i ++ )
		free ( seqs_id[i] );
	free ( seqs_id );

        free ( sw . input_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );
	free ( buf );
        free ( D );
        free(bufSC);
        free(SC);
        //Buffers
        free(bufDMawHighestOccMat);
        free(bufScMawHighestOccMat);
        free(bufHistMatrix);
        free(DMAWbufHistMatrix);
        free(SCMAWbufHistMatrix);
        //Matrices
        free(dMawHighestOccMat);
        free(scMawHighestOccMat);
        free(histMatrix);
        free(DMAWhistMatrix);
        free(SCMAWhistMatrix);
        //Filenames 
        free(dmaw_output_dist_filename);
        free(scmaw_output_dist_filename);
        free(scmaw_most_common_lengths);
        free(dmaw_most_common_lengths);
        free(union_most_common_lengths);
        free(card_ratios);
        free(union_card_ratios);
        free(length_ratios);
        free(union_length_ratios);
        free(hist_seqs);
        free(SCMAW_hist_seqs);
        free(DMAW_hist_seqs);
        free(UNION_hist_seqs);

	return ( 0 );
}

void save_on_file_simple_matrix_double(char * filename, double ** matrix, INT num_seqs, unsigned char ** seqs_id) {

	FILE *file;

	if ( ! ( file = fopen ( filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", filename );
		return ;
	}

	fprintf(file, ",");
	for(int i = 0; i < num_seqs; i++) fprintf(file, "%s,", seqs_id[i]);
	for(int i = 0; i < num_seqs; i++) {
		fprintf(file, "\n");
		fprintf(file, "%s,", seqs_id[i]);
		for(int j = 0; j < num_seqs; j++) {
			fprintf(file, "%.2f,", matrix[i][j]);
		}

	}

	fclose(file);
}

void save_on_file_INT(char * filename, INT ** matrix, INT num_seqs, unsigned char ** seqs_id) {

    FILE *file;

    if ( ! ( file = fopen ( filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", filename );
		return ;
	}

	fprintf(file, ",");
	for(int i = 0; i < num_seqs; i++) fprintf(file, "%s,", seqs_id[i]);
	for(int i = 0; i < num_seqs; i++) {
        fprintf(file, "\n");
        fprintf(file, "%s,", seqs_id[i]);
        for(int j = 0; j < num_seqs; j++) {
            fprintf(file, "%ld,", matrix[i][j]);
        }

	}

	fclose(file);

}

void save_on_file_double(char * filename, double ** matrix, INT num_seqs, unsigned char ** seqs_id) {

    FILE *file;

    if ( ! ( file = fopen ( filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", filename );
		return ;
	}

	fprintf(file, ",");
	for(int i = 0; i < num_seqs; i++) fprintf(file, "%s,", seqs_id[i]);
	for(int i = 0; i < num_seqs; i++) {
        fprintf(file, "\n");
        fprintf(file, "%s,", seqs_id[i]);
        for(int j = 0; j < num_seqs; j++) {
            fprintf(file, "%.2f%%,", matrix[i][j] * 100);
        }

	}

    int cont = 0;
	double media = 0.0, deviazione_standard = 0.0;

	for(int i = 0; i < num_seqs; i++) {
        for(int j = i + 1; j < num_seqs; j++) {
            media += matrix[i][j] * 100;
            cont++;
        }
	}

	media /= cont;

	for(int i = 0; i < num_seqs; i++) {
        for(int j = i + 1; j < num_seqs; j++) {
            double numeratore = matrix[i][j] * 100 - media;
            numeratore = numeratore * numeratore;
            deviazione_standard += numeratore;
        }
	}

	deviazione_standard /= cont;

	deviazione_standard = sqrt(deviazione_standard);

	fprintf(file, "\n\n\nMean (above the main diagonal): %.2f", media);
	fprintf(file, "\nStandard deviation (above the main diagonal): %.2f", deviazione_standard);

	fclose(file);
}

void save_on_file_INT_histograms(char * filename, INT ** hist_matrix, INT num_seqs, INT k, INT K, unsigned char ** seqs_id) {

    FILE *file;

    if ( ! ( file = fopen ( filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", filename );
		return ;
	}

    fprintf(file, ",");
    for(int i = 0; i < (K - k + 1); i++) fprintf(file, "Len-%ld,", i + k);
	for(int i = 0; i < num_seqs; i++) {
        fprintf(file, "\n%s,", seqs_id[i]);
        for(int j = 0; j < (K - k + 1); j++) {
            fprintf(file, "%ld,", hist_matrix[i][j]);
        }
	}

	fclose(file);
}

void save_on_file_INT_histograms_complex(char * filename, INT ** matrix, INT num_seqs, INT k, INT K, unsigned char ** seqs_id) {

    FILE *file;

    if ( ! ( file = fopen ( filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", filename );
		return ;
	}

    fprintf(file, ",");
    for(int i = 0; i < (K - k + 1); i++) fprintf(file, "Len-%ld,", i + k);
	for(int i = 0; i < num_seqs; i++) {
        for(int X = i + 1; X < num_seqs; X++) {
            fprintf(file, "\n%s_%s,", seqs_id[i], seqs_id[X]);
            for(int j = 0; j < (K - k + 1); j++) {
                fprintf(file, "%ld,", matrix[i * num_seqs + X][j]);
            }
        }
	}

	fclose(file);
}
