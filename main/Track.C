#include "Geometry.h"
#include "Tracking.h"

int main(){

    srand(time(0));
    
    Tracking Track;

    Track.Propagation();

    Track.Drawing();

    return 0;
}