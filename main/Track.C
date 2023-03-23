#include "Geometry.h"
#include "Tracking.h"

int main(){

    srand(time(0));
    
    auto Track= new Tracking();

    Track->Propagate();

    Track->Draw();

    return 0;
}
