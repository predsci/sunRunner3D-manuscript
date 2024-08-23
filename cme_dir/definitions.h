#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LimO3
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            12

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   EIGHT_WAVES
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  GAMMA                          0
#define  R0                             1
#define  OMEGA_ROT                      2
#define  RHO_PERT                       3
#define  V_PERT                         4
#define  CME_START_TIME                 5
#define  CME_RAMP                       6
#define  CME_DURATION                   7
#define  THETA0                         8
#define  PHI0                           9
#define  CME_RAD                        10
#define  B_PERT                         11

/* [Beg] user-defined constants (do not change this line) */

#define  RING_AVERAGE                   8
#define  CHAR_LIMITING                  YES
#define  LIMITER                        MINMOD_LIM
#define  VTK_VECTOR_DUMP                YES
#define  VTK_TIME_INFO                  YES
#define  CME_V                          YES
#define  CME_RHO                        YES
#define  CME_B                          YES

/* [End] user-defined constants (do not change this line) */
