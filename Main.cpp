#include "Simulation.h"

int main(){
    SIM* sim = new SIM();
    for(int i=0;i<sim->bins;i++){
        sim->singleUp();
    }
}
