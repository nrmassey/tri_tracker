#include "set_cout_precision.h"

void set_cout_precision(std::ostream& out, int places)
{
    // set precision and format
    std::ios_base::fmtflags ff = out.flags();
    out.flags((ff & std::cout.fixed) | std::cout.scientific);
    out.precision(places);
}
