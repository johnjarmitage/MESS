
void resetbasegrid(double *grid);
void markertorhogrid(int e,double value,double *grid,double *weight);
void movemarkers1();
void markermomentumparamstogrid();
void printstate(int name);
void stressstrainupdate();
void stressrotation();

double rhogridtomarker(int e,double *grid);
double secdevinv2(double tensorxx,double tensorxz);
double maxvelocity();

int get_s(double xpos);
int get_t(double zpos);