/*

 DESCRIPTION:
   Sample program that illustrates how to use a GA to find the maximum value
of a continuous function in two variables.  This program uses a binary-to-
decimal genome.
---------------------------------------------------------------------------- */
#include <stdio.h>
#include <ga/ga.h>
#include <math.h>
#include <ga/std_stream.h>

#define cout STD_COUT

float objective(GAGenome &);

int
main(int argc, char **argv)
{
  cout << "Example 9\n\n";
  cout << "This program finds the minimum value in the function\n";
  cout << "  f(x1,x2) \n";
  cout << "with the constraints\n";
  cout << "     -a <= x1 <= b\n";
  cout << "     -c <= x2 <= d\n";
  cout << "\n\n"; cout.flush();

// See if we've been given a seed to use (for testing purposes).  When you
// specify a random seed, the evolution will be exactly the same each time
// you use that seed number.

  unsigned int seed = 0;
  for(int i=1; i<argc; i++) {
    if(strcmp(argv[i++],"seed") == 0) {
      seed = atoi(argv[i]);
    }
  }

// Declare variables for the GA parameters and set them to some default values.

  int popsize  = 500;
  int ngen     = 100;
  float pmut   = 0.01;
  float pcross = 0.6;

// Create a phenotype for two variables.  The number of bits you can use to
// represent any number is limited by the type of computer you are using.  In
// this case, we use 16 bits to represent a floating point number whose value
// can range from -5 to 5, inclusive.  The bounds on x1 and x2 can be applied
// here and/or in the objective function.

  GABin2DecPhenotype map;

/* --------------------- Uncomment the required Function Limits --------------*/

  //map.add(16, -4.5, 4.5);  //Beale Function Limits
  //map.add(16, -4.5, 4.5);

  //map.add(16, -10, 10);    //Booth Function Limits
  //map.add(16, -10, 10);

  //map.add(16, -10, 10);    //Matyas Function Limits
  //map.add(16, -10, 10);

  //map.add(16, -1.5, 1.5);  //Rosenbrock Function Limits
  //map.add(16, -0.5, 2.5);

   map.add(16, -10, 0);	   //Mishra's Bird Function Limits
   map.add(16, -6.5, 0);
 
  //map.add(16,-2,2);      //Spherical Function Limits
  //map.add(16,-2,2);

  //map.add(16,-1.5,4);   //Mc-Cormick Function Limits
  //map.add(16,-3,4);

  map.add(16, -5, 5);
  map.add(16, -5, 5);   //Tangs Function Limits


// Create the template genome using the phenotype map we just made.

  GABin2DecGenome genome(map, objective);

// Now create the GA using the genome and run it.  We'll use sigma truncation
// scaling so that we can handle negative objective scores.

  GASimpleGA ga(genome);
  GASigmaTruncationScaling scaling;
  ga.populationSize(popsize);
  ga.nGenerations(ngen);
  ga.pMutation(pmut);
  ga.pCrossover(pcross);
  ga.scaling(scaling);
  ga.scoreFilename("bog.dat");
  ga.scoreFrequency(10);
  ga.flushFrequency(50);
  ga.evolve(seed);

// Dump the results of the GA to the screen.

  genome = ga.statistics().bestIndividual();
  cout << "the ga found an optimum at the point (";
  cout << genome.phenotype(0) << ", " << genome.phenotype(1) << ")\n\n";
  cout << "best of generation data are in '" << ga.scoreFilename() << "'\n";

  return 0;
}
 

// This objective function tries to maximize the value of the function
//
//                  y = -f(x1,x2)
//
float
objective(GAGenome & c)
{
  GABin2DecGenome & genome = (GABin2DecGenome &)c;

  float y;

/* --------------------- Uncomment the required Function  ---------------------*/

  //y = -(pow(1.5-genome.phenotype(0)+genome.phenotype(0)*genome.phenotype(1),2)+pow(2.25-genome.phenotype(0)+genome.phenotype(0)*pow(genome.phenotype(1),2),2)+pow(2.625-genome.phenotype(0)+genome.phenotype(0)*pow(genome.phenotype(1),3),2));	   // Beale Function

 // y = -(pow(genome.phenotype(0)+2*genome.phenotype(1)-7,2)+pow(2*genome.phenotype(0)+genome.phenotype(1)-5,2));  //Booth's Function

 // y=-(0.26*(pow(genome.phenotype(0),2)+pow(genome.phenotype(1),2))-0.48*genome.phenotype(0)*genome.phenotype(1));  // Matyas Function

 // y=-(pow(1-genome.phenotype(0),2) + 100*pow(genome.phenotype(1)-pow(genome.phenotype(0),2),2));  //Rosenbrock Function

  y=-(sin(genome.phenotype(1))*exp(pow(1-cos(genome.phenotype(0)),2))+cos(genome.phenotype(0))*exp(pow(1-sin(genome.phenotype(0)),2)) + pow(genome.phenotype(0)-genome.phenotype(1),2)); 		//Mishra's Bird Function


  //y=-(pow(genome.phenotype(0),2) + pow(genome.phenotype(1),2)); //Spherical Function
  
  //y=-(sin(genome.phenotype(0) + genome.phenotype(1)) + pow((genome.phenotype(0) - genome.phenotype(1)),2) - 1.5*genome.phenotype(0) +2.5*genome.phenotype(1) +1);

  //y = -(pow(genome.phenotype(0),4)-16*pow(genome.phenotype(0),2)+5*genome.phenotype(0)+ pow(genome.phenotype(1),4)-16*pow(genome.phenotype(1),2)+5*genome.phenotype(1))/2;       //Tangs Function
  return y;
}
