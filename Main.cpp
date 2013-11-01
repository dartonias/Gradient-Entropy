#include "Simulation.h"

int main(){
    SIM* sim = new SIM();
    for(int i=0;i<sim->bins;i++){
        for(int j=0;j<sim->Eq;j++){
            sim->sweep();
        }
        sim->resetMeasures();
        for(int j=0;j<sim->MCS;j++){
            sim->sweep();
        }
        sim->printMeasures();
    }
}
