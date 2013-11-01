#include "Simulation.h"

int main(){
    SIM* sim = new SIM();
    for(int i=0;i<sim->bins;i++){
        sim->singleUp();
        if((i%10)==0){
            sim->wolff();
            sim->print();
        }
    }
}
