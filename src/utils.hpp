/* the tool used in setup.cpp will be written here*/
#include "setup.hpp"

using namespace std;

string id2object(int id, float3 rotation, LBM & lbm);
void run_simulation(float si_u, int id, float3 rotation, float size, bool is_train, bool render );
void run_simulation_2d(float si_u, int id, float size, bool is_train, bool render);