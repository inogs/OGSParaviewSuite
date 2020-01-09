/* <<< Release Notice for library >>> */

#include "proj.hpp"
#include "projects.hpp"

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

char const pj_release[] =
    "Rel. "
    STR(PROJ_VERSION_MAJOR)"."
    STR(PROJ_VERSION_MINOR)"."
    STR(PROJ_VERSION_PATCH)", "
    "September 15th, 2018";

const char *pj_get_release() {
    return pj_release;
}