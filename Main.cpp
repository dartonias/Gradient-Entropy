#include "Simulation.h"

int main(){
    SIM* sim = new SIM();
    for(int j=0;j<sim->Eq;j++){
            sim->sweep();
    }
    for(int i=0;i<sim->bins;i++){
        sim->resetMeasures();
        for(int j=0;j<sim->MCS;j++){
            sim->sweep();
        }
        sim->printMeasures();
    }
}
