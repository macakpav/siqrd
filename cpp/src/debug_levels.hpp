#ifndef DEBUGLEVELS_HPP
#define DEBUGLEVELS_HPP
/*
    Defines debugging levels
*/

// NINFO - to supress all informative outputs
#ifdef NDEBUG
#undef DODESOLVER // debug ode solver 
#undef DMETHODS // debug methods for ode solving
#undef DLVL3 // maybe too much of debugging info
#undef DLVL2 // deeper debuging
#undef DLVL1 // some high level debug info
#undef DLVL0 // prints command line input
#endif

#ifdef DLVL3
#define DLVL2
#define DLVL1 
#define DLVL0
#endif
#ifdef DLVL2
#define DLVL1
#define DLVL0
#endif
#ifdef DLVL1
#define DLVL0
#endif
#endif