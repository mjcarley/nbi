#include <stdio.h>

#include <petscversion.h>

int main(int argc, char **argv)

{
  fprintf(stdout, "%d.%d\n", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR) ;
  
  return 0 ;
}
