#ifndef MESODYNxH
#define MESODYNxH
#include "namics.h"
#include "input.h"
#include "system.h"
#include "newton.h"
#include <random>

class Mesodyn {
public:
	Mesodyn(vector<Input*>,vector<Lattice*>,vector<Segment*>,vector<Molecule*>,vector<System*>,vector<Newton*>,string);

~Mesodyn();
	void AllocateMemory();
	void gaussianNoise(Real, Real, unsigned long);
	void langevinFlux();

	vector<Real> noise;

	string name;
	vector<Input*> In;
	vector<Lattice*> Lat;
	vector<Segment*> Seg;
	vector<Molecule*> Mol;
	vector<System*> Sys;
	vector<Newton*> New;
	string brand;

	vector<string> ints;
	vector<string> Reals;
	vector<string> bools;
	vector<string> strings;
	vector<Real> Reals_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;
	void push(string,Real);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	Real* GetPointer(string);
	int GetValue(string,int&,Real&,string&);
	int timesteps;
	int timebetweensaves;


	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	bool CheckInput(int);
	void PutParameter(string);
	string GetValue(string);


};
#endif