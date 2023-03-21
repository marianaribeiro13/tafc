#include "Geometry.h"
#include "Tracking.h"

int main(){

    srand(time(0));
    
    Tracking Track;

    Track.Propagate();

    Track.Draw();

    return 0;
}