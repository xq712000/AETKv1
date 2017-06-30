//#include  <iostream>
#include "IOstreams.H"
//#include  <cmath>
//#include  <fstream>
//#include <cstdlib>
//using namespace std;

	#include "1RINV.H"	
	
	#include "SquareMatrix.H"
	#include "vector.H"

using namespace Foam;
/*
  void rinv::input ()      //
  {
	  int  i, j;
	  char str1[20];
	  cout <<" \n input file name:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n file can't open " <<str1 <<endl; exit(1); }
	  for (i=0; i<n; i++)                       //
		  for (j=0; j<n; j++)  fin >>a[i][j];
	  fin.close ();
  }
*/

  void rinv::input_matrix (const SquareMatrix<scalar>& matrix )      //
  {
	  int  i, j;
//	  char str1[20];
//	  cout <<"\input file name:  ";
//	  cin >>str1;
//	  ifstream  fin (str1);
//	  if (!fin)
//	  { cout <<"\n file can't open " <<str1 <<endl; exit(1); }
		if(matrix.n() != n)
		{
			cout << "\n n!=matrix.n";
		}
	  for (i=0; i<n; i++)                       //
		  for (j=0; j<n; j++)  a[i][j]=matrix[i][j];
//	  fin.close ();
  }

  void rinv::inv ()          //
  { 
//	  int *is,*js,i,j,k;
	  int is[5],js[5],i,j,k;
      double d,p;
//      is = new int[n];
//      js = new int[n];
      for (k=0; k<=n-1; k++)
      { 
		  d=0.0;
          for (i=k; i<=n-1; i++)
          for (j=k; j<=n-1; j++)
          { 
			  p=fabs(a[i][j]);
              if (p>d) { d=p; is[k]=i; js[k]=j;}
          }
          if (d+1.0==1.0)
          { 
//			  delete [] is,js;
//              cout <<"\nsingularity matrix " <<endl;
              Info <<"\nsingularity matrix " <<endl;
			  exit(1);
          }
          if (is[k]!=k)
            for (j=0; j<=n-1; j++)
            { 
                p=a[k][j]; a[k][j]=a[is[k]][j]; a[is[k]][j]=p;
            }
          if (js[k]!=k)
            for (i=0; i<=n-1; i++)
            { 
                p=a[i][k]; a[i][k]=a[i][js[k]]; a[i][js[k]]=p;
            }
          a[k][k]=1.0/a[k][k];
          for (j=0; j<=n-1; j++)
            if (j!=k)  a[k][j]=a[k][j]*a[k][k];
          for (i=0; i<=n-1; i++)
            if (i!=k)
              for (j=0; j<=n-1; j++)
                if (j!=k) a[i][j]=a[i][j]-a[i][k]*a[k][j];
          for (i=0; i<=n-1; i++)
            if (i!=k)  a[i][k]=-a[i][k]*a[k][k];
      }
      for (k=n-1; k>=0; k--)
      { 
		  if (js[k]!=k)
            for (j=0; j<=n-1; j++)
            { 
                p=a[k][j]; a[k][j]=a[js[k]][j]; a[js[k]][j]=p;
            }
          if (is[k]!=k)
            for (i=0; i<=n-1; i++)
            { 
                p=a[i][k]; a[i][k]=a[i][is[k]]; a[i][is[k]]=p;
            }
      }
//      delete [] is, js;
  }

/*
  void rinv::output ()       //
  {
	  int  i, j;
	  char str2[20];
	  cout <<"\ninput output file name:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n file can't open " <<str2 <<endl; exit(1); }
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<"    " <<a[i][j];
			  cout <<"    " <<a[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }
*/
  void rinv::output_matrix ( SquareMatrix<scalar>& matrix)
  {
	  int  i, j;
//	  char str2[20];
//	  cout <<"\ninput output file name:  ";
//	  cin >>str2;
//	  ofstream fout (str2);
//	  if (!fout)
//	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<n; j++)
		  {
//			  fout <<"    " <<a[i][j];
//			  cout <<"    " <<a[i][j];
				matrix[i][j]=a[i][j];
		  }
//		  fout <<endl;  cout <<endl;
	  }
//	  fout.close ();
  }

//	int main ()      //
//  {

//    SquareMatrix<scalar> hmm(3);

//    hmm[0][0] = -3.0;
//    hmm[0][1] = 10.0;
//    hmm[0][2] = -4.0;
//    hmm[1][0] = 2.0;
//    hmm[1][1] = 3.0;
//    hmm[1][2] = 10.0;
//    hmm[2][0] = 2.0;
//    hmm[2][1] = 6.0;
//    hmm[2][2] = 1.0;

//    SquareMatrix<scalar> hmm2(3);

//    hmm2[0][0] = 0;
//    hmm2[0][1] = 0;
//    hmm2[0][2] = 0;
//    hmm2[1][0] = 0;
//    hmm2[1][1] = 0;
//    hmm2[1][2] = 0;
//    hmm2[2][0] = 0;
//    hmm2[2][1] = 0;
//    hmm2[2][2] = 0;

//	  rinv  c(3);
//	  c.input_matrix( hmm);         //
//		c.inv ();
//		c.output_matrix (hmm2);
//		Info<< hmm2 << endl;

//  }
