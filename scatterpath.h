
#define	RAD_PER_DEG	0.01745329
#define	PI	3.1415927
#define	ON	1
#define	OFF	0
#define	RSURF	6371.0
#define	ROUTC	3479.5
#define	RINC	1217.5
#define MAXCORDS	10000
#define	TINY	0.1
#define	UNDEF	-1

struct cords
{
	   float x;
	   float y;
};

void downgoing(double p, double r0, double v0, int lay0);
void upgoing(double p, double r0, double v0, int lay0);
void straight(double p, double r0, double v0, int lay0);
void write_ray(FILE *fp,int np,struct cords *raycords);
void trace(double p,double rtop,double rbot,double vtop,double vbot,double blay,double *delta,double *tt,double *path);
void read_vel(FILE *fp,int npt,double *rad,double *vel);
void write_time(FILE *fp,double dis_sc2,double p_1,double trat,float sclat,float sclon,double sc_dep,double delp);
void write_moreoutput(int lay);
void find_b(double *r,double *v,double *b);
void find_source(double *rad,double *vel,double *b,double source_dep,
        int *source_lay,double *source_rad,double *source_vel);

void find_ref(double *rad,float ref1_depth,int lay1,
        float ref2_depth,int lay2,
        float ref3_depth,int lay3);

void find_deltt(double p,double *del,double *tt);
void downgoing_fast(double p, double r0, double v0, int lay0);
void upgoing_fast(double p, double r0, double v0, int lay0);
void straight_fast(double p, double r0, double v0, int lay0);
void find_p(double *rad,double *vel,double *b,double p_in,int lay_st,int lay_end,double st_rad,double st_vel,double end_rad,double end_vel,float dis_sc1,double *p1,double *delp,double *trat,double *pathp);
void find_scatter(double *rad,double *vel,double *b,double sc_dep,int *sc_lay,double *sc_rad,double *sc_vel);
