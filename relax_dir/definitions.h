#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LimO3
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            3

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
#define  ROTATING_FRAME                 YES

/* -- user-defined parameters (labels) -- */

#define  GAMMA                          0
#define  R0                             1
#define  OMEGA_ROT                      2

/* [Beg] user-defined constants (do not change this line) */

#define  RING_AVERAGE                   NO
#define  CHAR_LIMITING                  YES
#define  LIMITER                        MINMOD_LIM
#define  VTK_VECTOR_DUMP                NO
#define  VTK_TIME_INFO                  NO

/* [End] user-defined constants (do not change this line) */
