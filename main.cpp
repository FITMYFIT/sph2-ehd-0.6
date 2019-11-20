#include "SPHSolver.h"

int main()
{

  cout<<endl;
  cout<<"*****************************************"<<endl;
  cout<<"*                                       *"<<endl;
  cout<<"*          SPH_2D  C++ Edition          *"<<endl;
  cout<<"*             by Richard Liu            *"<<endl;
  cout<<"*                                       *"<<endl;
  cout<<"*****************************************"<<endl;
  cout<<endl;


  CSPHSolver sphsolver;
  sphsolver.Input();
  sphsolver.Run();
  sphsolver.Done();
}
