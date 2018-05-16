#include "mesodyn.h"

//Constructor
//TODO: Read D, noise seed from file.
Mesodyn::Mesodyn(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Newton*> New_, string name_)
    : name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sys{Sys_},
      New{New_},
      D{0.5},
      mean{1},
      stdev{1},
      timesteps{0},
      timebetweensaves{0},
      zNeighbor{1}, // Usage for e.g. layer z: foo[z]+foo[z+1] becomes foo[z] + foo[z+xNeighbor]
      yNeighbor{Lat[0]->MZ},
      xNeighbor{Lat[0]->MZ * Lat[0]->MY},
      cNeighbor{Lat[0]->M},                                     // Neighboring component
      componentNo{(int)In[0]->MolList.size()},                  //find how many compontents there are (e.g. head, tail, solvent)
      size{componentNo * Lat[0]->MX * Lat[0]->MY * Lat[0]->MZ}, //find out how large the density vector is (needed for sizing the flux vector)
                                                                //which will be 1 flux per lattice site per component per dimension
      initRho{0.5},                                             // default is a homogeneous system.
      dimensions{3}                                             //TODO: Get this number from somewhere
{
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  KEYS.push_back("diffusionconstant");
  KEYS.push_back("seed");
  KEYS.push_back("mean");
  KEYS.push_back("stdev");
}

Mesodyn::~Mesodyn() {
}

void Mesodyn::AllocateMemory() {
  //alocate memory for the fluxes of all components in all dimensions
  try {
    J.resize(dimensions * size);
    L.reserve(combinations(componentNo, 2) * Lat[0]->M);
    rho.reserve(size);
  } catch (...) {
    cout << "Failed to reserve enough memory. System too large for RAM?";
    abort();
  }
  if (debug)
    cout << "nothing to allocate in Mesodyn" << endl;
}

bool Mesodyn::CheckInput(int start) {
  if (debug)
    cout << "Check Mesodyn" << endl;
  bool success = true;

  success = In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES);
  if (success) {
    vector<string> options;
    if (GetValue("timesteps").size() > 0) {
      success = In[0]->Get_int(GetValue("timesteps"), timesteps, 1, 10000, "The number of timesteps should be between 1 and 10000");
    }
    cout << "Timesteps is " << timesteps << endl;

    if (GetValue("timebetweensaves").size() > 0) {
      success = In[0]->Get_int(GetValue("timebetweensaves"), timebetweensaves);
    }
    cout << "Timesteps is " << timesteps << endl;

    if (GetValue("diffusionconstant").size() > 0) {
      success = In[0]->Get_Real(GetValue("diffusionconstant"), D);
    }
    cout << "Diffusion const is " << D << endl;

    if (GetValue("seed").size() > 0) {
      success = In[0]->Get_Real(GetValue("seed"), seed);
    }
    cout << "Seed is " << seed << endl;

    if (GetValue("seed").size() > 0) {
      success = In[0]->Get_Real(GetValue("mean"), mean);
    }
    cout << "Mean is " << mean << endl;

    if (GetValue("seed").size() > 0) {
      success = In[0]->Get_Real(GetValue("stdev"), stdev);
    }
    cout << "Stdev is " << stdev << endl;
  }

  return success;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  AllocateMemory(); //this HAS to be done before fillRho
  fillRho(initRho);
  if (success) {
    cout << "Mesodyn is all set, starting calculations.." << endl;
    //TODO: GUESS FOR U?
    for (int t = 0; t < timesteps; t++) {      // get segment potentials by iteration so that they match the given rho.
      New[0]->Solve(&rho[0], &dummyVector[0]); // TODO: DUMMY VECTOR!
      onsagerCoefficient();
      langevinFlux();
      updateDensity();
    }
  }
  return success;
}

//defaults to homogeneous system for now
void Mesodyn::fillRho(Real givenDensity) {
  for (int i = 0; i < size; ++i) {
    rho[i] = givenDensity;
  }
  for (int i = 1; i <= componentNo; ++i) {
    ptrComponentStart.push_back(&rho[Lat[0]->M * i]);
  }
}

void Mesodyn::abort() {
  //Once false is returned, Mesodyn automatically quits to main.
  success = false;
}

/******** Calculations ********/
void Mesodyn::onsagerCoefficient() {

  //TODO: maybe do this inline / per J calculation to preserve memory
  vector<Real>::iterator lIterator;
  lIterator = L.begin();

  //TODO: Untested code
  //all combinations of ptrComponentStart multiplications
  for (int i = 0; i < componentNo - 1; ++i) {
    for (int j = i + 1; j < componentNo; ++j) {
      //at every coordinate (pointer arithmatic)
      for (int xyz = 0; xyz < (Lat[0]->M); ++xyz) {
        *lIterator = (*(ptrComponentStart[i] + xyz) * *(ptrComponentStart[j] + xyz));
        ++lIterator;
      }
    }
  }
}

void Mesodyn::langevinFlux() {
  //TODO: safer to change size calculation to xx.size()?
  vector<Real> u(size); //segment potential A

  //TODO: boundary condition in lattice?
  int z = 1;
  u[z] = New[0]->xx[z]; //which alphas? component a & b or one component at two sites?

  vector<Real>::iterator jIterator;
  jIterator = J.begin();

  for (int i = 0; i <= componentNo; ++i) {
    for (int z = 0; z <= Lat[0]->M; ++z) {
      gaussianNoise(dummyMean, dummyStdev, 1);
      //something like: for the number of onsager's coefficients, calculate flux according to x, y and z.
      *jIterator = -D * ((L[z] + L[z + xNeighbor]) * (u[z + xNeighbor] - u[z])) - ((L[z - xNeighbor] + L[z]) * (u[z] - u[z - xNeighbor])) + noise[0];
      ++jIterator;
    }
    for (int z = 0; z <= Lat[0]->M; ++z) {
      gaussianNoise(dummyMean, dummyStdev, 1);
      *jIterator = -D * ((L[z] + L[z + yNeighbor]) * (u[z + yNeighbor] - u[z])) - ((L[z - yNeighbor] + L[z]) * (u[z] - u[z - yNeighbor])) + noise[0];
      ++jIterator;
    }
    for (int z = 0; z <= Lat[0]->M; ++z) {
      gaussianNoise(dummyMean, dummyStdev, 1);
      *jIterator = -D * ((L[z] + L[z + zNeighbor]) * (u[z + zNeighbor] - u[z])) - ((L[z - zNeighbor] + L[z]) * (u[z] - u[z - zNeighbor])) + noise[0];
      ++jIterator;
    }
  }
}

inline Real Mesodyn::at(int x, int y, int z, int c) {
  return 1; //[ z*Lat[0]->MZ*Lat[0]->MY + y*Lat[0]->MX + x ];
}

void Mesodyn::updateDensity() {
  //old density + langevinFluxTwo
}

/* Generates a vector of length count, contianing gaussian noise of given mean, standard deviation.
	 Noise is stored in vector<Real> Mesodyn::noise
	 Possible errors: What if count > sizeof(unsinged long)?
	 Called by langevinFlux()
*/
void Mesodyn::gaussianNoise(Real mean, Real stdev, unsigned long count) {

  random_device generator;

  seed_seq seed("something", "something else");

  //Mersenne Twister 19937 bit state size PRNG
  mt19937 prng(seed);

  normal_distribution<> dist(mean, stdev);

  this->noise.resize(count);

  for (unsigned int i = 0; i < count; ++i) {
    this->noise[i] = dist(prng);
  }

  /* Debugging code (output value in all elements):
for (auto const &element: mesodyn.thisNoise)
				std::cout << element << ' ';
*/
}

/******* Tools ********/
int Mesodyn::factorial(int n) {
  if (n > 1) {
    return n * factorial(n - 1);
  } else
    return 1;
}

int Mesodyn::combinations(int n, int k) {
  return factorial(n) / (factorial(n - k) * factorial(k));
}

void Mesodyn::PutParameter(string new_param) {
  KEYS.push_back(new_param);
}
string Mesodyn::GetValue(string parameter) {
  int i = 0;
  int length = PARAMETERS.size();
  while (i < length) {
    if (parameter == PARAMETERS[i]) {
      return VALUES[i];
    }
    i++;
  }
  return "";
}
void Mesodyn::push(string s, Real X) {
  Reals.push_back(s);
  Reals_value.push_back(X);
}
void Mesodyn::push(string s, int X) {
  ints.push_back(s);
  ints_value.push_back(X);
}
void Mesodyn::push(string s, bool X) {
  bools.push_back(s);
  bools_value.push_back(X);
}
void Mesodyn::push(string s, string X) {
  strings.push_back(s);
  strings_value.push_back(X);
}
void Mesodyn::PushOutput() {
  strings.clear();
  strings_value.clear();
  bools.clear();
  bools_value.clear();
  Reals.clear();
  Reals_value.clear();
  ints.clear();
  ints_value.clear();
}
Real* Mesodyn::GetPointer(string s) {
  //vector<string> sub;
  //nothing yet
  return NULL;
}
int Mesodyn::GetValue(string prop, int& int_result, Real& Real_result, string& string_result) {
  int i = 0;
  int length = ints.size();
  while (i < length) {
    if (prop == ints[i]) {
      int_result = ints_value[i];
      return 1;
    }
    i++;
  }
  i = 0;
  length = Reals.size();
  while (i < length) {
    if (prop == Reals[i]) {
      Real_result = Reals_value[i];
      return 2;
    }
    i++;
  }
  i = 0;
  length = bools.size();
  while (i < length) {
    if (prop == bools[i]) {
      if (bools_value[i])
        string_result = "true";
      else
        string_result = "false";
      return 3;
    }
    i++;
  }
  i = 0;
  length = strings.size();
  while (i < length) {
    if (prop == strings[i]) {
      string_result = strings_value[i];
      return 3;
    }
    i++;
  }
  return 0;
}
