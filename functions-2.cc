 /**
    MAW: Minimal Absent Words  v. 1.0
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>

#include "mawdefs.h"

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                           // include header for suffix sort
#endif

#include <sdsl/bit_vectors.hpp>					  // include header for bit vectors
#include <sdsl/rmq_support.hpp>					  //include header for range minimum queries
#include "stack.h"        					  // include header for stack structure


#include <stdbool.h>

using namespace sdsl;
using namespace std;

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};

unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP )
{
	INT i=0, j=0;

	LCP[0] = 0;
	for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
		if ( ISA[i] != 0 )  //// i non indica la posizione nell'iSA del suffisso lessicograficamente minore
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )
				j++;
			LCP[ISA[i]] = j;
		}

	return ( 1 );
}

unsigned int compute_maw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, TMaw ** Occ, unsigned int * NOcc )
{
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT n = strlen ( ( char * ) seq );
	int sigma;
    	bit_vector * Before;
    	bit_vector * Beforelcp;

	if      ( ! strcmp ( "DNA", sw . alphabet ) )   sigma = strlen ( ( char * ) DNA );
        else if ( ! strcmp ( "PROT", sw . alphabet ) )  sigma = strlen ( ( char * ) PROT );
        else if ( ! strcmp ( "SIXTEEN", sw . alphabet ) )  sigma = strlen ( ( char * ) SIXTEEN );

        /* Compute the suffix array */
        SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

	#ifdef _USE_64
        if( divsufsort64( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_32
        if( divsufsort( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

        /*Compute the inverse SA array */
        invSA = ( INT * ) calloc( n , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < n; i ++ )
        {
                invSA [SA[i]] = i;
        }

	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP array */
        if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }

	free ( invSA );

	INT v_size = 2 * n;

    	Before = new bit_vector[sigma];
        if( ( Before == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for Before.\n" );
                return ( 0 );
        }
	Beforelcp = new bit_vector[sigma];
        if( ( Beforelcp == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for BeforeLCP.\n" );
                return ( 0 );
        }
    	for ( INT i = 0; i < sigma; i++ )
	{
		Before[i] = bit_vector( v_size, 0 );
		Beforelcp[i] = bit_vector( v_size, 0 );
    	}

	GetBefore ( seq, n , sigma, SA, LCP, Before, Beforelcp );

	GetMaws( seq, seq_id, SA, n, sigma, LCP, Before, Beforelcp, sw . k, sw . K, sw . output_filename, Occ, NOcc );

	//printf("NOcc: %ld", *NOcc);
	delete [] Before;
	delete [] Beforelcp ;
	free ( SA );
	free ( LCP );

	return ( 1 );
}


unsigned char Mapping( int a )
{
	char c = DEL;
        switch ( a )
	{
            case 0:
                c = 'A';
                break;
            case 1:
                c = 'C';
                break;
            case 2:
                c = 'G';
                break;
            case 3:
                c = 'T';
                break;
            case 4:
                c = 'N';
                break;
            case 5:
                c = 'R';
                break;
            case 6:
                c = 'D';
                break;
            case 7:
                c = 'Q';
                break;
            case 8:
                c = 'E';
                break;
            case 9:
                c = 'H';
                break;
            case 10:
                c = 'I';
                break;
            case 11:
                c = 'L';
                break;
            case 12:
                c = 'K';
                break;
            case 13:
                c = 'M';
                break;
            case 14:
                c = 'F';
                break;
            case 15:
                c = 'P';
                break;
            case 16:
                c = 'S';
                break;
            case 17:
                c = 'W';
                break;
            case 18:
                c = 'Y';
                break;
            case 19:
                c = 'V';
                break;
                //Aggiunte da me
            case 20:
                c = 'B';
                break;
            case 21:
                c = 'J';
                break;
            case 22:
                c = 'O';
                break;
        }
	return ( c );
}

int RevMapping ( unsigned char b )
{
	int a = -1;
        switch ( b )
	{
            case 'A':
                a = 0;
                break;
            case 'C':
                a = 1;
                break;
            case 'G':
                a = 2;
                break;
            case 'T':
                a = 3;
                break;
            case 'N':
                a = 4;
                break;
            case 'R':
                a = 5;
                break;
            case 'D':
                a = 6;
                break;
            case 'Q':
                a = 7;
                break;
            case 'E':
                a = 8;
                break;
            case 'H':
                a = 9;
                break;
            case 'I':
                a = 10;
                break;
            case 'L':
                a = 11;
                break;
            case 'K':
                a = 12;
                break;
            case 'M':
                a = 13;
                break;
            case 'F':
                a = 14;
                break;
            case 'P':
                a = 15;
                break;
            case 'S':
                a = 16;
                break;
            case 'W':
                a = 17;
                break;
            case 'Y':
                a = 18;
                break;
            case 'V':
                a = 19;
                break;
            case 'B':
                a = 20;
                break;
            case 'J':
                a = 21;
                break;
            case 'O':
                a = 22;
                break;
        }
	return ( a );
}

/* computes the reverse complement of str */
unsigned int RevComStr ( unsigned char * str, unsigned char * str2, INT iLen )
{
   INT i = 0;
   while ( iLen -- )
    {
      switch ( str[iLen] )
       {
         case 'A':
           str2[i++] = 'T';
           break;
         case 'C':
           str2[i++] = 'G';
           break;
         case 'G':
           str2[i++] = 'C';
           break;
         case 'T':
           str2[i++] = 'A';
           break;
         case 'N':
           str2[i++] = 'N';
           break;
         default:
           return ( 0 );
       }
    }
   return ( 1 );
}


unsigned int GetBefore (
				unsigned char * seq,
				INT n,
				int sigma,
				INT * SA,
				INT * LCP,
				bit_vector * Before,
				bit_vector * Beforelcp )
{
        INT hm = 0;
        INT k = 0;
        INT lcp;
        INT mem;
        INT proxa;
        INT proxb;

        TStack lifo_lcp;
        StackNew ( &lifo_lcp, sizeof( INT ) );
        TStack lifo_mem;
        StackNew ( &lifo_mem, sizeof( INT ) );
        TStack lifo_rem;
        StackNew ( &lifo_rem, sizeof( INT ) );

        lcp = 0;
        StackPush(&lifo_lcp, &lcp);

        /* Max LCP value */
        for ( int i = 0; i < n; i++ )
                if( LCP[i] > hm )
                        hm = LCP[i];
        hm = hm + 2;

        bit_vector* interval = new bit_vector[sigma];
        for ( INT i = 0; i < sigma; i++)
                interval[i]=bit_vector(hm,0);

	interval[RevMapping(seq[n- 1])][0]=1;

        // First pass : top-down
        for ( INT i = 0; i < n; i++ )
        {
                // first we update the interval table
                // we empty the interval that corresponds to a higher lcp value
                if ( i > 0 && LCP[i] < LCP[i-1])
                {
                        StackPop(&lifo_lcp,&lcp);
                        while(!StackEmpty(&lifo_lcp)&&lcp>LCP[i])
                        {
                                StackPop(&lifo_lcp,&mem);
                                if (mem <=LCP[i])
                                {
                                        for (int j=0; j<sigma; j++)
                                        {
                                                if (mem!=LCP[i]){interval[j][LCP[i]]=interval[j][lcp];}     //initialisation of the next intervals if it hasn't been open
                                                Before[j][2*i-1]=interval[j][lcp];
                                                Beforelcp[j][2*i-1]=interval[j][lcp];
                                                if (mem==LCP[i]){Beforelcp[j][2*i-1]=interval[j][mem];}
                                        }
                                }

                                for (int j =0; j<sigma ;j++){interval[j][lcp]=0;}
                                        lcp=mem;
                        }
                        StackPush(&lifo_lcp,&lcp);
                }

        	// we update those having a lower lcp
                if ( SA[i] - 1 >= 0 )
                        k = RevMapping(seq[SA[i] - 1]);
                else
                        k = - 1;
                if ( k != -1 )
                {
                        while(!StackEmpty(&lifo_lcp))
                        {
                                StackPop(&lifo_lcp,&lcp);
                                StackPush(&lifo_mem, &lcp);
                                if (interval[k][lcp]==1){break;}
                                interval[k][lcp]=1;
                        }

                        while (!StackEmpty(&lifo_mem))
                        {
                                StackPop(&lifo_mem,&lcp);
                                StackPush(&lifo_lcp, &lcp);
                        }

                        interval[k][LCP[i]]=1;
                }

                if ( ( i - 1 ) >= 0 && SA[i - 1] - 1 >= 0 )
                        if (i>0 && LCP[i]>0 && RevMapping(seq[SA[i - 1] - 1])!=-1) // in this case we also add the letter preceding the last suffix
                        {
                                interval[RevMapping(seq[SA[i - 1] - 1])][LCP[i]]=1;
                        }

                // Before, Before_LCP
                for (int j =0; j<sigma ;j++)
                {
                        Beforelcp[j][2*i]=interval[j][LCP[i]];
                }

                if (k!=-1)
                {
                        Before[k][2*i+1]=1;
                        Before[k][2*i]=1;
                        Beforelcp[k][2*i+1]=1;
                        Beforelcp[k][2*i]=1;
                }

                StackPop(&lifo_lcp,&lcp);
                if(lcp!=LCP[i])  // no duplicates
                {
                        StackPush(&lifo_lcp,&lcp);
                        lcp=LCP[i];
                }
                StackPush(&lifo_lcp,&lcp);
        }

        //second pass : bottom-up
        //we empty the interval table

        while(!StackEmpty(&lifo_lcp))
        {
                StackPop(&lifo_lcp,&lcp);
                for (int j =0; j<sigma ;j++){interval[j][lcp]=0;}
        }
        lcp=0;
        StackPush(&lifo_lcp,&lcp);

        for (INT i=n-1; i>-1; i--)
        {
                StackPop(&lifo_lcp,&lcp);
                proxa=LCP[i]+1;   //proxb is the lcp-value that is just higher than LCP[i]
                while(!StackEmpty(&lifo_lcp) && lcp>LCP[i])
                {
                        StackPush(&lifo_rem,&lcp);
                        StackPop(&lifo_lcp,&mem);
                        if (mem <LCP[i])            //initialisation of the interval if it hasn't been open
                        {
                                for (int j=0; j<sigma; j++){interval[j][LCP[i]]=interval[j][lcp];}
                                proxa=lcp;
                        }
                        if (mem ==LCP[i]){proxa=lcp;}
                        lcp=mem;
                }
                StackPush(&lifo_lcp,&lcp);

                // we update the lower intervals
                for (int k=0; k<sigma;k++)
                {
                        if(Before[k][2*i]==1)
                        {
                                while(!StackEmpty(&lifo_lcp))
                                {
                                        StackPop(&lifo_lcp,&lcp);
                                        StackPush(&lifo_mem, &lcp);
                                        if (interval[k][lcp]==1){break;}
                                        interval[k][lcp]=1;
                                }

                                while (!StackEmpty(&lifo_mem))
                                {
                                        StackPop(&lifo_mem,&lcp);
                                        StackPush(&lifo_lcp, &lcp);
                                }
                                interval[k][LCP[i]]=1;
                        }
                }

                for (int j =0; j<sigma ;j++)
                { //update interval and Before
                        Beforelcp[j][2*i]=Beforelcp[j][2*i] || interval[j][LCP[i]];

                        if(i<n-1)
                        {
                                Before[j][2*i+1]=Before[j][2*i+1] || interval[j][proxb];  //proxb is the lcp-value that is just higher than LCP[i+1]
                                Beforelcp[j][2*i+1]=interval[j][LCP[i+1]] || Beforelcp[j][2*i+1];
                        }
                }

                proxb=proxa;

                //we suppress higher intervals
                if (i<n-1 && LCP[i+1]>LCP[i])
                {
                        StackPop(&lifo_rem,&lcp);  // this lcp is the one that is just higher than LCP[i]
                        for (int j =0; j<sigma ;j++)
                        {
                                Before[j][2*i]=Before[j][2*i] || interval[j][lcp];
                                interval[j][lcp]=0;
                        }

                        while(!StackEmpty(&lifo_rem))
                        {
                                StackPop(&lifo_rem,&lcp);
                                for (int j =0; j<sigma ;j++){interval[j][lcp]=0;}

                        }
                }

                StackPop(&lifo_lcp,&lcp);
                if(lcp!=LCP[i])  // no duplicates
                {
                        StackPush(&lifo_lcp,&lcp);
                        lcp=LCP[i];
                }
                StackPush(&lifo_lcp,&lcp);
        }

        delete[] interval;
        StackDispose(&lifo_lcp);
        StackDispose(&lifo_mem);
        StackDispose(&lifo_rem);

	return ( 1 );
}

unsigned int GetMaws( unsigned char * seq, unsigned char * seq_id, INT * SA, INT n, int sigma, INT * LCP, bit_vector* Before, bit_vector* Beforelcp, unsigned int k, unsigned int K, char * out_file, TMaw ** Occ, unsigned int * NOcc )
{
    	FILE * out_fd;
	char * maw;

	// compute a bitvector that contains a `1', if an identical row has already been seen => to avoid duplicates.
    	bit_vector Mem = bit_vector(n,0);
    	TStack lifo_lcp;
	StackNew (&lifo_lcp, sizeof( INT ) );
	INT lcp = 0;
	INT mem;
	StackPush(&lifo_lcp, &lcp);

	for ( INT i = 0; i < n; i++ )
	{
        	StackPop(&lifo_lcp,&lcp);
            	while(!StackEmpty(&lifo_lcp)&&lcp>LCP[i])
            	{
                	StackPop(&lifo_lcp,&mem);
                	if ( mem == LCP[i] )
                	{
                    		Mem[i] = 1;
                	}
                	lcp = mem;
            	}
            	StackPush(&lifo_lcp,&lcp);
            	lcp = LCP[i];
            	StackPush(&lifo_lcp,&lcp);
	}
	StackDispose(&lifo_lcp);

	#if 0
	if ( ! ( out_fd = fopen ( out_file, "a") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", out_file );
		return ( 1 );
	}
	#endif

	/* Print the header */
	//fprintf ( out_fd, ">%s\n", ( char * ) seq_id );

        maw = ( char * ) calloc( ( K + 1 ) , sizeof( char ) );
        if( ( maw == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory.\n" );
                return ( 0 );
        }

	vector<mytuple> xtuple;

        for ( INT i = 0; i < n; i++ )
    	{
        	for( int l = 0; l < sigma; l++ )
        	{
            		bool B1 = (
					Before[l][2 * i] == 0 &&
					Beforelcp[l][2 * i] == 1 &&
					SA[i] + LCP[i] < n &&
					LCP[i] + 2 <= K	&&
					LCP[i] + 2 >= k	);

            		bool B2 = (
					i < n - 1 &&
					Before[l][2 * i + 1] == 0 &&
					Beforelcp[l][2 * i + 1] == 1 &&
					SA[i] + LCP[i + 1] < n &&
					LCP[i + 1] + 2 <= K &&
					LCP[i + 1] + 2 >= k &&
					Mem[i + 1] == 0 );

			/* Here we report MAWs */
		    	if ( B1 )
		    	{
				maw[0] = Mapping( l );
				INT start = SA[i];
				INT size = SA[i]+ LCP[i] + 1 - start;
				memcpy( &maw[1], &seq[start], size );
				maw[size + 1] = '\0';
				xtuple . push_back(make_tuple(( INT ) maw[0], i, start, size ));

		    	}
		    	else if ( B2 )
		    	{
				maw[0] = Mapping( l );
				INT start = SA[i];
				INT size = SA[i] + LCP[i + 1] + 1 - start;
				memcpy( &maw[1], &seq[start], size );
				maw[size + 1] = '\0';
				xtuple . push_back(make_tuple(( INT ) maw[0], i, start, size ));
		    	}
        	}

    	}

    	std::sort( xtuple.begin(), xtuple.end() );

	for(vector<mytuple>::iterator iter = xtuple.begin(); iter != xtuple.end(); iter++)
	  ( * NOcc ) ++;


	( * Occ ) = ( TMaw * ) realloc ( ( * Occ ),   ( ( * NOcc ) ) * sizeof ( TMaw ) );

	INT j = 0;
	for(vector<mytuple>::iterator iter = xtuple.begin(); iter != xtuple.end(); iter++)
	{
		( *Occ )[j] . letter = get<0>(*iter);
		( *Occ )[j] . pos = get<2>(*iter);
		( *Occ )[j] . size = get<3>(*iter);
		j++;
  	}

	#if 0
	fprintf( out_fd, "\n" );

	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	#endif

	free ( maw );

	return ( 1 );
}

unsigned int maw_seq_comp ( unsigned char * X, TMaw * mawX, unsigned int * NmawX, unsigned char * Y, TMaw * mawY, unsigned int * NmawY, double * XY_distance, INT k, INT K, INT * scMawHighestOcc, INT * scNumber, INT * scTotalLength, INT pos1, INT pos2, INT ** SCMAWHist, INT num_seqs )
{
	unsigned char * XY;
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT m = strlen ( ( char * ) X );
	INT n = strlen ( ( char * ) Y );
	INT N = m + n;

	XY = ( unsigned char * ) malloc( ( N + 1 ) * sizeof( unsigned char ) );

	strncpy ( ( char * ) XY, ( char * ) X, m );
	strncpy ( ( char * ) &XY[m], ( char * ) Y, n );
	XY[ m + n ] = '\0';

	/* Compute the suffix array */
        SA = ( INT * ) malloc( ( N ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

	#ifdef _USE_64
        if( divsufsort64( XY, SA,  N ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_32
        if( divsufsort( XY, SA,  N ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

        /*Compute the inverse SA array */
        invSA = ( INT * ) calloc( N , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < N; i ++ )
        {
                invSA [SA[i]] = i;
        }

	LCP = ( INT * ) calloc  ( N, sizeof( INT ) );
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP array */
        if( LCParray( XY, N, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }

        int_vector<> v( N , 0 ); // create a vector of length n and initialize it with 0s

	for ( INT i = 0; i < N; i ++ )
	{
		v[i] = LCP[i];
	}

	rmq_succinct_sct<> rmq(&v);

	util::clear(v);

#if 0
	for ( i = 0; i < (* NmawX); i++ )
	{
	  fprintf(stderr, " <%c,%d,%d>\n", mawX[i] . letter, mawX[i] . pos, mawX[i] . size );
	}
	fprintf(stderr, "\n" );
	for ( j = 0; j < (* NmawY); j++ )
	{
	  fprintf(stderr, " <%c,%d,%d>\n", mawY[j] . letter, mawY[j] . pos, mawY[j] . size );
	}
#endif
	INT i = 0;
	INT j = 0;

	INT * scMawHighestOccArr;

    scMawHighestOccArr = (INT *) calloc((K - k + 1), sizeof(INT));

	INT number = 0;

	INT totalLength = 0;

	while ( i < (* NmawX) && j < (* NmawY) )
	{
	  if ( mawX[i] . letter < mawY[j] . letter )
	  {
	      double u = mawX[i] . size + 1;
	      double one = 1;
	      ( * XY_distance ) += one / ( u * u );
	      scMawHighestOccArr[(int) u - k]++;
	      number++;
	      totalLength += (int) u;
	      SCMAWHist[pos1 * num_seqs + pos2][(int) u - k]++;
	      i++;
	  }
	  else if ( mawX[i] . letter > mawY[j] . letter )
	  {
	      double u = mawY[j] . size + 1;
	      double one = 1;
	      ( * XY_distance ) += one / ( u * u );
	      scMawHighestOccArr[(int) u - k]++;
	      number++;
	      totalLength += (int) u;
	      SCMAWHist[pos1 * num_seqs + pos2][(int) u - k]++;
	      j++;
	  }
	  else if ( mawX[i] . letter == mawY[j] . letter )
	  {
	    INT ix = mawX[i] . pos;
	    INT iy = m + mawY[j] . pos;
	    INT i_rank = invSA[ix];
	    INT j_rank = invSA[iy];
	    INT l = SOLONmin ( invSA[ ix ], invSA[ iy ] );
            INT r = SOLONmax ( invSA[ ix ], invSA[ iy ] );
	    INT LCE = LCP[rmq ( l + 1, r ) ];

	    if ( LCE >= mawX[i] . size && mawX[i] . size == mawY[j] . size )
	    {
	      i++; j++;
	    }
	    else
	    {
	      if ( XY[ix + LCE] < XY[iy + LCE] )
	      {
		double u = mawX[i] . size + 1;
		double one = 1;
	        ( * XY_distance ) += one / ( u * u );
	        scMawHighestOccArr[(int) u - k]++;
	        number++;
	        totalLength += (int) u;
	        SCMAWHist[pos1 * num_seqs + pos2][(int) u - k]++;
		i++;
	      }
	      else
	      {
	        double u = mawY[j] . size + 1;
		double one = 1;
	        ( * XY_distance ) += one / ( u * u );
	        scMawHighestOccArr[(int) u - k]++;
	        number++;
	        totalLength += (int) u;
	        SCMAWHist[pos1 * num_seqs + pos2][(int) u - k]++;
		j++;
	      }
	    }
	  }
	}

	for ( ; i < (* NmawX); i++ )
	{
	  double u = mawX[i] . size + 1;
	  double one = 1;
	  ( * XY_distance ) += one / ( u * u );
	  scMawHighestOccArr[(int) u - k]++;
	  number++;
	  totalLength += (int) u;
	  SCMAWHist[pos1 * num_seqs + pos2][(int) u - k]++;
	}

	for ( ; j < (* NmawY); j++ )
	{
	   double u = mawY[j] . size + 1;
	   double one = 1;
	   ( * XY_distance ) += one / ( u * u );
	   scMawHighestOccArr[(int) u - k]++;
	   number++;
	   totalLength += (int) u;
	   SCMAWHist[pos1 * num_seqs + pos2][(int) u - k]++;
	}

	INT maximum = 0;
    *scMawHighestOcc = 0;
    //printf("\n");
    for(INT w = 0; w < (K - k + 1); w++) {
        //printf("%ld\t", dMawHighestOccArr[w]);
        if(scMawHighestOccArr[w] > maximum) {
            maximum = scMawHighestOccArr[w];
            *scMawHighestOcc = w + k;
        }
    }

	//printf("\nSOLON scNumber: %ld\n", number);

	*scNumber = number;
	*scTotalLength = totalLength;

    free ( SA );
	free ( LCP );
	free ( invSA );
	free ( XY );
  	return ( 1 );
}

/* Edit distance */
unsigned int edit_distance ( unsigned char * x, unsigned char * y, double * XY_distance )
{
	INT i, j;
	INT sub = 1;
	INT gap = 1;

	INT m = strlen ( ( char * ) x );
	INT n = strlen ( ( char * ) y );

	INT ** D;
	D = ( INT ** ) calloc ( m + 1, sizeof ( INT * ) );
	for ( i = 0; i < m + 1; i++ )
    	{
      		D[i] = ( INT * ) calloc ( n + 1, sizeof ( INT ) );
    	}

	for ( i = 1; i < m + 1; i++ )
    	{
      		D[i][0] = i * gap;
    	}
	for ( j = 1; j < n + 1; j++ )
    	{
      		D[0][j] = j * gap;
    	}

	for ( j = 1; j < n + 1; j++ )
    	{
		for ( i = 1; i < m + 1; i++ )
        	{
          		if ( x[i - 1] == y[j - 1] )
			{
            			D[i][j] = D[i - 1][j - 1];      //a match
			}
          		else
			{
				if ( D[i-1][j-1] + sub <= D[i][j-1] + gap && D[i-1][j-1] + sub <= D[i-1][j] + gap )
				{
					D[i][j] = D[i-1][j-1] + sub;
				}
				else if ( D[i-1][j] + gap <= D[i][j-1] + gap && D[i-1][j] + gap <= D[i-1][j-1] + sub )
				{
					D[i][j] = D[i - 1][j] + gap;
				}
				else if ( D[i][j-1] + gap <= D[i - 1][j] + gap && D[i][j-1] + gap <= D[i - 1][j-1] + sub )
				{
					D[i][j] = D[i][j-1] + gap;
				}
			}
        	}
    	}

	( * XY_distance ) = ( double ) D[m][n];

	for ( i = 0; i < m + 1; i++ )
    	{
      		free ( D[i] );
    	}
	free ( D );

	return ( 1 );
}

unsigned int maw_seq_comp_factors ( unsigned char * X, TMaw * mawX, unsigned int * NmawX, unsigned char * Y, TMaw * mawY, unsigned int * NmawY, double * XY_distance, INT k, INT K, INT * dMawHighestOcc, INT * dNumber, INT * dTotalLength, INT pos1, INT pos2, INT ** DMAWHist, INT num_seqs, INT * unionHighestOcc, INT * unionNumber, INT * unionTotalLength, INT ** UNIONHist)
{
	unsigned char * XY;
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT m = strlen ( ( char * ) X );
	INT n = strlen ( ( char * ) Y );
	INT N = m + n + 1;

	XY = ( unsigned char * ) malloc( ( N + 1 ) * sizeof( unsigned char ) );

	strncpy ( ( char * ) XY, ( char * ) X, m );
    XY[m] = '#';
	strncpy ( ( char * ) &XY[m + 1], ( char * ) Y, n );
	XY[N] = '\0';

	/* Compute the suffix array */
        SA = ( INT * ) malloc( ( N ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

	#ifdef _USE_64
        if( divsufsort64( XY, SA,  N ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_32
        if( divsufsort( XY, SA,  N ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

        /*Compute the inverse SA array */
        invSA = ( INT * ) calloc( N , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < N; i ++ )
        {
                invSA [SA[i]] = i;
        }

	LCP = ( INT * ) calloc  ( N, sizeof( INT ) );
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP array */
        if( LCParray( XY, N, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }

        int_vector<> v( N , 0 ); // create a vector of length n and initialize it with 0s

	for ( INT i = 0; i < N; i ++ )
	{
		v[i] = LCP[i];
	}

	rmq_succinct_sct<> rmq(&v);

	util::clear(v);

#if 0
	for ( i = 0; i < (* NmawX); i++ )
	{
	  fprintf(stderr, " <%c,%d,%d>\n", mawX[i] . letter, mawX[i] . pos, mawX[i] . size );
	}
	fprintf(stderr, "\n" );
	for ( j = 0; j < (* NmawY); j++ )
	{
	  fprintf(stderr, " <%c,%d,%d>\n", mawY[j] . letter, mawY[j] . pos, mawY[j] . size );
	}
#endif

#if 0
    printf("MAWs of X=\"%s\":\n", X);
	for ( INT i = 0; i < (* NmawX); i++ )
	{
	  fprintf(stderr, " <%c,%ld,%ld> = \"%c", mawX[i] . letter, mawX[i] . pos, mawX[i] . size, mawX[i].letter );
      for(INT k = 0; k < mawX[i].size; k++) printf("%c", X[mawX[i].pos + k]);
      printf("\"\n");
	}
    printf("MAWs of Y=\"%s\":\n", Y);
	for ( INT j = 0; j < (* NmawY); j++ )
	{
	  fprintf(stderr, " <%c,%ld,%ld> = \"%c", mawY[j] . letter, mawY[j] . pos, mawY[j] . size, mawY[j].letter );
      for(INT k = 0; k < mawY[j].size; k++) printf("%c", Y[mawY[j].pos + k]);
      printf("\"\n");
	}

    printf("\nX=\"%s\", Y=\"%s\" \nX#Y=\"%s\"", X, Y, XY);
    printf("\nSuffix Array of X#Y: ");
    for(int i = 0; i < N; i++) printf("%ld ", SA[i]);
    printf("\nInverse Suffix Array of X#Y: ");
    for(int i = 0; i < N; i++) printf("%ld ", invSA[i]);
    printf("\nLCP Array of X#Y: ");
    for(int i = 0; i < N; i++) printf("%ld ", LCP[i]);
    printf("\n\n\n");
#endif

    /*INT * scMawHighestOccArr;

    scMawHighestOccArr = (INT *) calloc((K - k + 1), sizeof(INT));


    //printf("MAWs of X: %d\nMAWs of Y: %d\n", *NmawX, *NmawY);
	while ( i < (* NmawX) && j < (* NmawY) )
	{

	  if ( mawX[i] . letter < mawY[j] . letter )
	  {
          if(!factor_of(mawX[i], true, m, XY, N, SA, invSA, LCP)) {i++; continue;}
	      double u = mawX[i] . size + 1;
	      double one = 1;
	      ( * XY_distance ) += one / ( u * u );
	      i++;
	  }
	  else if ( mawX[i] . letter > mawY[j] . letter )
	  {
          if(!factor_of(mawY[j], false, m, XY, N, SA, invSA, LCP)) {j++; continue;}
	      double u = mawY[j] . size + 1;
	      double one = 1;
	      ( * XY_distance ) += one / ( u * u );
	      j++;
	  }
	  else if ( mawX[i] . letter == mawY[j] . letter )
	  {
	    INT ix = mawX[i] . pos;
	    INT iy = m + mawY[j] . pos;
	    INT i_rank = invSA[ix];
	    INT j_rank = invSA[iy];
	    INT l = SOLONmin ( invSA[ ix ], invSA[ iy ] );
            INT r = SOLONmax ( invSA[ ix ], invSA[ iy ] );
	    INT LCE = LCP[rmq ( l + 1, r ) ];

	    if ( LCE >= mawX[i] . size && mawX[i] . size == mawY[j] . size )
	    {
	      i++; j++;
	    }
	    else
	    {
	      if ( XY[ix + LCE] < XY[iy + LCE] )
	      {
            if(!factor_of(mawX[i], true, m, XY, N, SA, invSA, LCP)) {i++; continue;}
		double u = mawX[i] . size + 1;
		double one = 1;
	        ( * XY_distance ) += one / ( u * u );
		i++;
	      }
	      else
	      {
            if(!factor_of(mawY[j], false, m, XY, N, SA, invSA, LCP)) {j++; continue;}
	        double u = mawY[j] . size + 1;
		double one = 1;
	        ( * XY_distance ) += one / ( u * u );
		j++;
	      }
	    }
	  }
	}

	for ( ; i < (* NmawX); i++ )
	{
      if(!factor_of(mawX[i], true, m, XY, N, SA, invSA, LCP)) {i++; continue;}
	  double u = mawX[i] . size + 1;
	  double one = 1;
	  ( * XY_distance ) += one / ( u * u );
	}

	for ( ; j < (* NmawY); j++ )
	{
       if(!factor_of(mawY[j], false, m, XY, N, SA, invSA, LCP)) {j++; continue;}
	   double u = mawY[j] . size + 1;
	   double one = 1;
	   ( * XY_distance ) += one / ( u * u );
	}*/

	/////MIOOOOOOOO
	/*while ( i < (* NmawX) && j < (* NmawY) )
	{
	  if ( mawX[i] . letter < mawY[j] . letter )
	  {
	      double u = mawX[i] . size + 1;
	      double one = 1;
	      ( * XY_distance ) += one / ( u * u );
	      scMawHighestOccArr[(int) u - k]++;
	      scNumber++;
	      i++;
	  }
	  else if ( mawX[i] . letter > mawY[j] . letter )
	  {
	      double u = mawY[j] . size + 1;
	      double one = 1;
	      ( * XY_distance ) += one / ( u * u );
	      scMawHighestOccArr[(int) u - k]++;
	      scNumber++;
	      j++;
	  }
	  else if ( mawX[i] . letter == mawY[j] . letter )
	  {
	    INT ix = mawX[i] . pos;
	    INT iy = m + mawY[j] . pos;
	    INT i_rank = invSA[ix];
	    INT j_rank = invSA[iy];
	    INT l = SOLONmin ( invSA[ ix ], invSA[ iy ] );
            INT r = SOLONmax ( invSA[ ix ], invSA[ iy ] );
	    INT LCE = LCP[rmq ( l + 1, r ) ];

	    if ( LCE >= mawX[i] . size && mawX[i] . size == mawY[j] . size )
	    {
	      i++; j++;
	    }
	    else
	    {
	      if ( XY[ix + LCE] < XY[iy + LCE] )
	      {
		double u = mawX[i] . size + 1;
		double one = 1;
	        ( * XY_distance ) += one / ( u * u );
	        scMawHighestOccArr[(int) u - k]++;
	        scNumber++;
		i++;
	      }
	      else
	      {
	        double u = mawY[j] . size + 1;
		double one = 1;
	        ( * XY_distance ) += one / ( u * u );
	        scMawHighestOccArr[(int) u - k]++;
	        scNumber++;
		j++;
	      }
	    }
	  }
	}

	for ( ; i < (* NmawX); i++ )
	{
	  double u = mawX[i] . size + 1;
	  double one = 1;
	  ( * XY_distance ) += one / ( u * u );
	  scMawHighestOccArr[(int) u - k]++;
	  scNumber++;
	}

	for ( ; j < (* NmawY); j++ )
	{
	   double u = mawY[j] . size + 1;
	   double one = 1;
	   ( * XY_distance ) += one / ( u * u );
	   scMawHighestOccArr[(int) u - k]++;
	   scNumber++;
	}*/

/* Union stats */

    INT ui = 0, uj = 0;

    INT * unionHighestOccArr = (INT *) calloc((K - k + 1), sizeof(INT));

    INT uNumber = 0, uTotalLength = 0;

    //printf("\nUnion X:\n");

    while(ui < (*NmawX)) {
        /*fprintf(stderr, " <%c,%ld,%ld> = \"%c", mawX[ui] . letter, mawX[ui] . pos, mawX[ui] . size, mawX[ui].letter );
        for(INT k = 0; k < mawX[ui].size; k++) printf("%c", X[mawX[ui].pos + k]);
        printf("\"\n");*/
        uNumber++;
        double u = mawX[ui] . size + 1;
        unionHighestOccArr[(int) u - k]++;
        uTotalLength += (int) u;
        UNIONHist[pos1 * num_seqs + pos2][(int) u - k]++;
        ui++;
    }
    //printf("\nUnion Y:\n");
    while(uj < (*NmawY)) {
        if(notIn(X, Y, mawY[uj], mawX, NmawX)) {
            /*fprintf(stderr, " <%c,%ld,%ld> = \"%c", mawY[uj] . letter, mawY[uj] . pos, mawY[uj] . size, mawY[uj].letter );
            for(INT k = 0; k < mawY[uj].size; k++) printf("%c", Y[mawY[uj].pos + k]);
            printf("\"\n");*/
            double u = mawY[uj] . size + 1;
            unionHighestOccArr[(int) u - k]++;
            uNumber++;
            uTotalLength += (int) u;
            UNIONHist[pos1 * num_seqs + pos2][(int) u - k]++;
        }
        uj++;
    }

     INT uMaximum = 0;
    *unionHighestOcc = 0;
    //printf("\n");
    for(INT w = 0; w < (K - k + 1); w++) {
        //printf("%ld\t", dMawHighestOccArr[w]);
        if(unionHighestOccArr[w] > uMaximum) {
            uMaximum = unionHighestOccArr[w];
            *unionHighestOcc = w + k;
        }
    }

    *unionNumber = uNumber;
    *unionTotalLength = uTotalLength;

/* DMAW Evaluation */

    INT i = 0;
    INT j = 0;

    INT indx, shift;

    INT * dMawHighestOccArr = (INT *) calloc((K - k + 1), sizeof(INT));

    INT number = 0;

    INT totalLength = 0;

    while(i < (* NmawX)) {
        indx = invSA[mawX[i] . pos]; //shift = 0;
        if(factor_of(mawX[i], indx, XY, N, SA, invSA, LCP)) {
            double u = mawX[i] . size + 1;
	        double one = 1;
	        ( * XY_distance ) += one / ( u * u );
	        dMawHighestOccArr[(int) u - k]++;
	        number++;
	        totalLength += (int) u;
	        DMAWHist[pos1 * num_seqs + pos2][(int) u - k]++;
	        //printf("\n(%c,%ld,%ld)", mawX[i].letter, mawX[i].pos, mawX[i].size);
        }
        i++;
    }

    while(j < (*NmawY)) {
        indx = invSA[mawY[j] . pos + m + 1]; //shift = m + 1;
        if(factor_of(mawY[j], indx, XY, N, SA, invSA, LCP)) {
            double u = mawY[j] . size + 1;
	        double one = 1;
	        ( * XY_distance ) += one / ( u * u );
	        dMawHighestOccArr[(int) u - k]++;
	        number++;
	        totalLength += (int) u;
	        DMAWHist[pos1 * num_seqs + pos2][(int) u - k]++;
	        //printf("\n(%c,%ld,%ld)", mawY[j].letter, mawY[j].pos, mawY[j].size);
        }
        j++;
    }

    //printf("}\n");
    INT maximum = 0;
    *dMawHighestOcc = 0;
    //printf("\n");
    for(INT w = 0; w < (K - k + 1); w++) {
        //printf("%ld\t", dMawHighestOccArr[w]);
        if(dMawHighestOccArr[w] > maximum) {
            maximum = dMawHighestOccArr[w];
            *dMawHighestOcc = w + k;
        }
    }

    *dNumber = number;
    *dTotalLength = totalLength;
    //printf("\n");

     free ( SA );
	free ( LCP );
	free ( invSA );
	free ( XY );
	//free(scMawHighestOccArr);
	free(dMawHighestOccArr);
  	return ( 1 );
}

bool factor_of(TMaw maw, INT index, unsigned char * XY, INT N, INT * SA, INT * iSA, INT * LCP) {
    /*INT index, shift;
    if(MAW_X == true) { index = iSA[maw.pos]; shift = 0; }
    else { index = iSA[maw.pos + X_length + 1]; shift = X_length + 1; }*/
    INT up = index, down = index + 1;

    while(up >= 0 && LCP[up] >= maw.size) {
        up--;
        if(SA[up] > 0 && XY[SA[up] - 1] == maw.letter) {
            /*printf("%c", maw.letter );
            for(INT k = 0; k < maw.size; k++) printf("%c", XY[maw.pos + k + shift]);
            printf(", ");*/
            return true;
        }
    }

    while(down < N && LCP[down] >= maw.size) {
        if(SA[down] > 0 && XY[SA[down] - 1] == maw.letter) {
            /*printf("%c", maw.letter );
            for(INT k = 0; k < maw.size; k++) printf("%c", XY[maw.pos + k + shift]);
            printf(", ");*/
            return true;
        }
        down++;
    }

    return false;
}

bool notIn(unsigned char * X, unsigned char * Y, TMaw oneMawY, TMaw * mawX, unsigned int * NmawX) {
    //return true;
    INT i = 0, j, k;

    bool phase = false;

    while(i < (*NmawX)) {
        phase = false;
        if(oneMawY.letter == mawX[i].letter && oneMawY.size == mawX[i].size) {
            j = oneMawY.pos; k = mawX[i].pos;
            for(;j < oneMawY.pos + oneMawY.size && k < mawX[i].pos + mawX[i].size; j++, k++) {
                if(Y[j] != X[k]) phase = true;
            }
            if(phase == false) {
                return phase;
            }
        }

        i++;
    }

    return true;
}


