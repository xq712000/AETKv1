/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    rinv

Description
    It is used to get inverse matrix of the 5*5 matrix.


    References:
	Xu S.L.. "Commonly used algorithm set(C++ language to described)" Fourth ed. 
	TSINGHUA UNIVERSITY PRESS, 2009    

Author
    Shen C.  All rights reserved.

SourceFiles
    1RINV.H
    1RINV.C
\*---------------------------------------------------------------------------*/
#ifndef RINV_H
#define RINV_H

//#include  <iostream>
#include "IOstreams.H"
//#include  <cmath>
//#include  <fstream>
//#include <cstdlib>
//using namespace std;

	#include "SquareMatrix.H"
	#include "vector.H"

using namespace Foam;

class  rinv
{
	private: 
    		int n;
//		double  **a;
			double a[5][5];
	public:
	    rinv (int nn)
		{
//			int i;
			int i,j;
			n = nn;
//	        a = new double*[n];   //
//	        for (i=0; i<n; i++) a[i] = new double[n];
	        for (i=0; i<n; i++)
				{
					for(j=0; j<n; j++)
					{
						a[i][j] = 1;
					}
				}				
				 
		}

//	    	void input ();    //
			void input_matrix(const SquareMatrix<scalar>&);
		    void inv ();      //
  //          void output ();   //
           	void output_matrix ( SquareMatrix<scalar>&); 
		    ~rinv ()          //
		    {
//				int i;
//			    for (i=0; i<n; i++) { delete [] a[i]; }
//			    delete [] a;
		    }
  };


//#ifdef NoRepository
//#   include "1RINV.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
