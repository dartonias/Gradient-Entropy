/*  Main class for the simulation.
    Should be able to have spatial variation in temperature
    and spatial variation in coupling eventually.
*/

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <assert.h>
#include "MersenneTwister.h"
#include <math.h>

// Class defining the bonds
class BOND{
    public:
        // The first spin of the bond
        int a;
        // The second spin of the bond
        int b;
        // The coupling of the bond
        double J;
        // The temperature of the bond
        double beta;
        // Constructor
        BOND();
        // Assignment
        void assign(int a, int b, double J, double beta);
};

BOND::BOND(){
}

void BOND::assign(int _a, int _b, double _J, double _beta){
    a = _a;
    b = _b;
    J = _J;
    beta = _beta;
}

// Main class containing the full simulation
class SIM{
    private:
        // Size of the lattice in the x-direction and y-direction
        int Lx,Ly;
        // Number of total spins
        int nSpins;
        // Assuming a gradient from betaLow to betaHigh on the lattice
        double betaLow, betaHigh;
        // Assuming a gradient from JLow to JHigh on the lattice
        double JLow, JHigh;
        // Random number seed
        int seed;
        // Random number generator
        MTRand* rand;
        // All the connections, stored independantly for every spin
        // i.e. spin[i] is connected to BOND elements neighbors[i][j] for j=1,2,3,...
        std::vector<std::vector<int> > neighbors;
        // BOND elements of the lattice, defining the connections, couplings, and local temperature
        std::vector<BOND> bonds;
        // Array containing the spin state of the system
        std::vector<int> spins;
        // Internal energy of the system -- should be the same as calcE()
        double Energy;
        // Object used for cluster update
        std::vector<int> cluster;
        // Accumulated measurements
        double mE;
    public:
        // Default constructor, will read from a parameter file
        SIM();
        // Reads the data from input file
        void readInput();
        // Build up all the neighbors with a gradient in coupling and strength
        void buildGradient();
        // Build up all the neighbors with a jump in coupling and strength
        void buildJump();
        // Calculate full energy
        double calcE();
        // Polarize the spins of the lattice
        void polarizeSpins();
        // Single spin update
        void singleUp();
        // Number of equilibriation, Monte Carlo steps, and total bins
        int Eq, MCS, bins;
        // Cluster update
        void wolff();
        // Function called by cluster update
        void addToCluster(int z);
        // Print the configuration
        void print();
        // Reset all the measurements
        void resetMeasures();
        // Print all the measurements to file
        void printMeasures();
        // One monte carlo sweep, corresponding to "nSpins" single updates and one "wolff" update`
        void sweep();
};

// Default constructor
SIM::SIM(){
    readInput();
    rand = new MTRand(seed);
    //buildGradient();
    buildJump();
    polarizeSpins();
    Energy = calcE();
    resetMeasures();
}

// Reads the input parameters from a file
void SIM::readInput(){
    std::string filename = "param.txt";
    // Garbage string for collecting plaintext in parameter file
    std::string g;
    std::fstream inFile(filename.c_str());
    inFile >> g >> Lx;
    inFile >> g >> Ly;
    inFile >> g >> betaLow;
    inFile >> g >> betaHigh;
    inFile >> g >> JLow;
    inFile >> g >> JHigh;
    inFile >> g >> Eq;
    inFile >> g >> MCS;
    inFile >> g >> bins;
    inFile >> g >> seed;
    inFile.close();

    // Input checking
    assert (Lx >= 2);
    assert (Ly >= 2);
    assert (betaLow <= betaHigh);
    assert (JLow <= JHigh);

    /* Debug for input
    std::cout<< Lx          << std::endl;
    std::cout<< Ly          << std::endl;
    std::cout<< betaLow     << std::endl;
    std::cout<< betaHigh    << std::endl;
    std::cout<< JLow        << std::endl;
    std::cout<< JHigh       << std::endl;
    std::cout<< Eq          << std::endl;
    std::cout<< MCS         << std::endl;
    std::cout<< bins        << std::endl;
    */
}

// Builds a rectangular lattice
// We will us open boundary conditions, to match better what we expect
// for any microscopic experimental systems
void SIM::buildGradient(){
    nSpins = Lx*Ly;

    std::vector<int> tvec;
    tvec.clear();
    neighbors.resize(nSpins,tvec);

    bonds.clear();
    BOND tBond;
    int bond_counter = 0;
    int s1, s2;
    double tJ, tB;
    for (int x=0;x<Lx;x++){
        for (int y=0;y<Ly;y++){
            tJ = JLow + x*1.0/Lx * (JHigh - JLow);
            tB = betaLow + x*1.0/Lx * (betaHigh - betaLow);

            // Bond to the right
            s1 = x + y*Lx;
            s2 = s1 + 1;
            if (x < Lx-1){
                tBond.assign(s1,s2,tJ,tB);
                bonds.push_back(tBond);
                neighbors[s1].push_back(bond_counter);
                neighbors[s2].push_back(bond_counter);
                bond_counter++;
            }

            // Upward bond
            s2 = s1 + Lx;
            if (y < Ly-1){
                tBond.assign(s1,s2,tJ,tB);
                bonds.push_back(tBond);
                neighbors[s1].push_back(bond_counter);
                neighbors[s2].push_back(bond_counter);
                bond_counter++;
            }
        }
    }

    /* Debugging for bonds
    for (int i=0;i<neighbors.size();i++){
        std::cout << i << " --> ";
        for (int j=0;j<neighbors[i].size();j++){
            std::cout << bonds[neighbors[i][j]].a << " "<< bonds[neighbors[i][j]].b <<", ";
        }
        std::cout << std::endl;
    }
    */
}

// Builds a rectangular lattice
// We will us open boundary conditions, to match better what we expect
// for any microscopic experimental systems
void SIM::buildJump(){
    nSpins = Lx*Ly;

    std::vector<int> tvec;
    tvec.clear();
    neighbors.resize(nSpins,tvec);

    bonds.clear();
    BOND tBond;
    int bond_counter = 0;
    int s1, s2;
    double tJ, tB;
    for (int x=0;x<Lx;x++){
        for (int y=0;y<Ly;y++){
            if(x<Lx/2){
                tJ = JLow;
                tB = betaLow;
            }
            else{
                tJ = JHigh;
                tB = betaHigh;
            }

            // Bond to the right
            s1 = x + y*Lx;
            s2 = s1 + 1;
            if (x < Lx-1){
                tBond.assign(s1,s2,tJ,tB);
                bonds.push_back(tBond);
                neighbors[s1].push_back(bond_counter);
                neighbors[s2].push_back(bond_counter);
                bond_counter++;
            }

            // Upward bond
            s2 = s1 + Lx;
            if (y < Ly-1){
                tBond.assign(s1,s2,tJ,tB);
                bonds.push_back(tBond);
                neighbors[s1].push_back(bond_counter);
                neighbors[s2].push_back(bond_counter);
                bond_counter++;
            }
        }
    }

    /* Debugging for bonds
    for (int i=0;i<neighbors.size();i++){
        std::cout << i << " --> ";
        for (int j=0;j<neighbors[i].size();j++){
            std::cout << bonds[neighbors[i][j]].a << " "<< bonds[neighbors[i][j]].b <<", ";
        }
        std::cout << std::endl;
    }
    */
}

// Full calculation of the energy, by looping over all the bonds
// Assumption is that bonds are satisfies by FERROMAGNETIC interactions
// and we are simulating an Ising model where satisfied, E = -1, unsatisfied, E = 0
double SIM::calcE(){
    double tE = 0;
    for (int i=0;i<bonds.size();i++){
        if (spins[bonds[i].a] == spins[bonds[i].b]) tE -= bonds[i].J;
    }
    return tE;
}

// Polarize the spins of the lattice
void SIM::polarizeSpins(){
    spins.assign(nSpins,0);
}

// Single spin update
void SIM::singleUp(){
    // Spin we are attempting to flip
    int z = rand->randInt(nSpins-1);
    // Change in probability, which depends on the local beta*dE for eacn bond
    double dP = 0.0;
    double dE = 0.0;
    for (int i=0;i<neighbors[z].size();i++){
        if (spins[bonds[neighbors[z][i]].a] == spins[bonds[neighbors[z][i]].b]) {
        dP += bonds[neighbors[z][i]].J*bonds[neighbors[z][i]].beta;
        dE += bonds[neighbors[z][i]].J;
        }
        else {
        dP -= bonds[neighbors[z][i]].J*bonds[neighbors[z][i]].beta;
        dE -= bonds[neighbors[z][i]].J;
        }
    }
    if (rand->randDblExc() < exp(-1.0*dP)){
        spins[z] = (spins[z]+1)%2;
        Energy += dE;
    }
    
    /* Debugging check
    if(fabs(Energy - calcE()) > 1e-5){
        std::cout << "Energy = " << Energy << std::endl;
        std::cout << "calcE = " << calcE() << std::endl;
        throw -1;
    }
    */
}

// Wolff cluster update, important for typical systems near T_c
void SIM::wolff(){
    // Reset the cluster object
    cluster.assign(nSpins,0);
    // Starting spin of the cluster
    int z = rand->randInt(nSpins-1);
    // Add spin to the cluster
    addToCluster(z);
    // Flip cluster, change enegy
    for(int i=0;i<nSpins;i++){
        if(cluster[i] == 1){
            spins[i] = (spins[i] + 1)%2;
        }
    }
    Energy = calcE();
}

// Recursive algorithm for adding spins to the cluster
void SIM::addToCluster(int z){
    cluster[z] = 1;
    int s1,s2;
    for(int i=0;i<neighbors[z].size();i++){
        s1 = bonds[neighbors[z][i]].a;
        s2 = bonds[neighbors[z][i]].b;
        if(cluster[s1] == 0){
            if(spins[s1] == spins[z]){
                if(rand->randDblExc() < (1 - exp(-1.0*bonds[neighbors[z][i]].beta*bonds[neighbors[z][i]].J))){
                    addToCluster(s1);
                }
            }
        }
        else if(cluster[s2] == 0){
            if(spins[s2] == spins[z]){
                if(rand->randDblExc() < (1 - exp(-1.0*bonds[neighbors[z][i]].beta*bonds[neighbors[z][i]].J))){
                    addToCluster(s2);
                }
            }
        }
    }
}

// Print the lattice
void SIM::print(){
    for(int y=0;y<Ly;y++){
        for(int x=0;x<Lx;x++){
            std::cout << spins[x + Lx*y] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Energy = " << Energy << std::endl << std::endl;
}

void SIM::resetMeasures(){
    mE = 0.0;
}

void SIM::printMeasures(){
    std::string filename = "bins.txt";
    std::fstream outFile(filename.c_str(), std::fstream::out | std::fstream::app);
    outFile << mE / MCS << std::endl;
    outFile.close();
}

void SIM::sweep(){
    // Update the system
    for(int i=0;i<nSpins/2;i++){
        singleUp();
    }
    wolff();
    for(int i=0;i<nSpins/2;i++){
        singleUp();
    }
    // Measure the system
    mE += Energy;
}
