/*  Main class for the simulation.
    Should be able to have spatial variation in temperature
    and spatial variation in coupling eventually.
*/

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

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
        // Number of equilibriation, Monte Carlo steps, and total bins
        int Eq, MCS, bins;
        // All the connections, stored independantly for every spin
        // i.e. spin[i] is connected to BOND elements neighbors[i][j] for j=1,2,3,...
        std::vector<std::vector<int> > neighbors;
        // BOND elements of the lattice, defining the connections, couplings, and local temperature
        std::vector<BOND> bonds;
    public:
        // Default constructor, will read from a parameter file
        SIM();
        // Reads the data from input file
        void readInput();
        // Build up all the neighbors
        void buildLatt();
};

// Default constructor
SIM::SIM(){
    readInput();
    buildLatt();
}

// Reads the input parameters from a file
void SIM::readInput(){
    char* filename = "param.txt";
    // Garbage string for collecting plaintext in parameter file
    std::string g;
    std::fstream inFile(filename);
    inFile >> g >> Lx;
    inFile >> g >> Ly;
    inFile >> g >> betaLow;
    inFile >> g >> betaHigh;
    inFile >> g >> JLow;
    inFile >> g >> JHigh;
    inFile >> g >> Eq;
    inFile >> g >> MCS;
    inFile >> g >> bins;
    inFile.close();

    // Input checking
    assert (Lx >= 2);
    assert (Ly >= 2);
    assert (betaLow < betaHigh);
    assert (JLow < JHigh);

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
void SIM::buildLatt(){
    nSpins = Lx*Ly;

    vector<int> tvec();
    tvec.erase();
    neighbors.resize(nSpins,tvec);

    bonds.erase();
    BOND tBond;
    int bond_counter = 0;
    int s1, s2;
    double tJ, tB;
    for (int x=0;x<Lx-1;x++){
        for (int y=0;y<Ly-1;y++){
            tJ = JLow + x*1.0/Lx * (JHigh - JLow);
            tB = betaLow + x*1.0/Lx * (betaHigh - betaLow);

            // Bond to the right
            s1 = x + y*Lx;
            s2 = s1 + 1;
            tBond.assign(s1,s2,tJ,tB);
            bonds.push_back(tBond);
            neighbor[s1].push_back(bond_counter);
            neighbor[s2].push_back(bond_counter);
            bond_counter++;

            // Upward bond
            s2 = s1 + Lx;
            tBond.assign(s1,s2,tJ,tB);
            bonds.push_back(tBond);
            neighbor[s1].push_back(bond_counter);
            neighbor[s2].push_back(bond_counter);
            bond_counter++;
        }
    }

    /* Debugging for bonds */
    for (int i=0;i<neighbors.size(),i++){
        // todo
    }
    /* */
}
