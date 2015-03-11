/* this is to cal the real raypath in 1-D model, in order to cal. the tt perturbation 
based on tomography result. Thus, fix output=0, moreourput=1, projection=1
*/
/* raytrace.c:  a simple ray traceer. Xiaodong Song, Aug 26, 1993. 
 * find ray path for a given distance in the given ray parameter range.
 * modified from trs.c. Program itereates until it finds a ray paramter
 * which gives the desired delta within a certain error. It then find 
 * the ray path and everything for that ray parameter.
 *
 * generate ray path using piecewise approximation v=ar^b
 *
 * model can contain discontinuities (thichness defined by TINY=0.1km)
 *
 * Definition of "layer": it starts from surface as 0 layer and inceases
 * with depth.  Each layer contains inclusively the top depth of 
 * the layer but does not contain the bottom depth of the lay. The bottom
 * depth belongs to the next layer.
 *
 * the program can do top-side and under-side reflections at two layers 
 * (depths). Both top-side and under-side reflections occurs at the top
 * of the layer. See example for reflection specification
 *
 */
/*
  The original program is modified by Xiaolong Ma to be used to calculate
  the scattered ray path. I divided the ray path into two parts. One is the 
  direct wave from the source to the scatterer, and the other is the scattered
  wave from the scatterer to the receiver.        
  Warning :       
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "scatterpath.h"
#include "distazsub.h"

int  ref1_lay=UNDEF;
int  ref2_lay=UNDEF;
int  ref3_lay=UNDEF;
char ref1_type[20]="";
char ref2_type[20]="";
char ref3_type[20]="";
char *ptr_ref1;
char *ptr_ref2;
char *ptr_ref3;

float output_rad=ROUTC;		/*output deltasum,ttsum for this radius*/
float delcmb_in,delcmb_out;	/*deltas crossing the output_rad*/
float ttcmb_in,ttcmb_out;	/*times crossing the output_rad*/

static int saveray=0;
static int output=0;
static int moreoutput=0;
static int cmb_in,cmb_out;	/*flags for saving delta across CMB*/
static int iteration=50;

static double *rad;
static double *vel;
static double *b;
static double p;
static double bot_rad;
struct cords *raycords;
static int ncords=0;
static int nlay;
static double ttsum;
static double deltasum;
static double pathsum;
static FILE *outfp;

double source_rad,source_dep=0.,source_vel;
int source_lay;
int sc_lay;
int st_lay;
int fflag;
double sc_vel,sc_rad,sc_dep=0;
double st_vel,st_rad,st_dep=0;
float sx,sy,sz,rx,ry,rz;
float delta=-999;
static float error=0.001;
int projection=0;	/*if =1, project ray path into real Earth in terms of 
				latitude and longitude*/

main(ac,av)
int ac;
char **av;
{
 if(ac<14) {
   fprintf(stderr,"Usage: scatterpath velfile st_val end_val stla stlo sclat sclon sc_depth evla evlo evdp flag phase \n");
   fprintf(stderr,"Output: sclat sclon sc_depth tscatter t1 t2 delta1 delta2\n");
   fprintf(stderr,
               "velfile:      the velocity model you like!\n\
 st_val:        starting p parameter value\n\
end_val:        ending p parameter value\n\
   stla:	    station latitude\n\
   stlo:	    station longitude\n\
  sclat:        scatter latitude\n\
  sclon:        scatter lonitude\n\
 sc_dep:        scatter depth\n\
   evla:	event latitude\n\
   evlo:	event longitude\n\
   evdp:	event depth in km\n\
   flag:    source side(1) or receiver side(0)\n\
  phase:    P(PKP,PKIKP) or PKiKP\n");
   exit(-1);
}
FILE *velfp, *rayfp;
	char *velfile,rayfile[80],outfile[80];
        float source_depth,sc_depth;
	float st_val=0., end_val=0.;
	float ref1_depth=UNDEF;
	float ref2_depth=UNDEF;
	float ref3_depth=UNDEF;
	double inc_val,inc_val0,inc_ang0;
	int i;
	int npt;
    int lay_st,lay_end;
    double st_rad,st_vel,end_rad,end_vel;
	double del0,del1,del,tt;
	double p0,p1,p_1,p_in;
        double path_q;
	float stla=0,stlo=150,lat,lon,az;
        float sclat=0,sclon=0;
        float evla=0,evlo=0;
        float dis_sc1,az_sc1,dis_sc2,az_sc2;
        char *ph;
        double delp,trat,pathp;
	char shooting[10]="p";

     
	sprintf(outfile, "OUTPUT");
	velfile=av[1];
	sscanf(av[2], "%f", &st_val);
	sscanf(av[3], "%f", &end_val);
	sscanf(av[4], "%f", &stla);
	sscanf(av[5], "%f", &stlo);
        sscanf(av[6], "%f", &sclat);
        sscanf(av[7], "%f", &sclon);
        sscanf(av[8], "%f", &sc_depth);
	sscanf(av[9], "%f", &evla);
	sscanf(av[10], "%f", &evlo);
	sscanf(av[11], "%f", &source_depth);
	sscanf(av[12], "%d", &fflag);
	ph = av[13];
	sscanf(av[13], "%s", ph);
    if(strcmp(ph,"PKiKP")==0)
        {
            ref1_depth = 5153.5;
            ref1_type[0]  = 't';
        }
    
//	printf("%s\n", velfile);
//	printf("%f %f %f %f %f %f %f\n", sclat,sclon, stla, stlo, evla, evlo, source_depth);

	/* calculate ray paths and travel times */
	velfp=fopen(velfile, "r");
	if(saveray){
		rayfp=fopen(rayfile, "w");
		fprintf(stderr,"ray output: %s\n",rayfile);
	}
	if(output|| moreoutput){
		outfp=fopen(outfile, "w");
		fprintf(stderr,"output file is: %s\n",outfile);
	}

	fscanf(velfp,"%d", &npt);
	nlay=npt-1;

   	if(( vel=(double *)malloc (npt*sizeof(double)) )==NULL){
   		fprintf(stderr, "memory allocation error!");
   		exit(-1);
	}
   	if(( rad=(double *) malloc(npt*sizeof(double)) )==NULL){
   		fprintf(stderr, "memory allocation error!");
   		exit(-1);
	}
   	if(( b=(double *) malloc(npt*sizeof(double)) )==NULL){
   		fprintf(stderr, "memory allocation error!");
   		exit(-1);
	}
   	if(( raycords=(struct cords *)malloc(MAXCORDS*sizeof(struct cords)) )==NULL){
   		fprintf(stderr, "memory allocation error!");
   		exit(-1);
	}

	source_dep=source_depth;
        sc_dep = sc_depth;
	read_vel(velfp,npt,rad,vel);
	find_b(rad,vel,b);
	find_source(rad,vel,b,source_dep,
		&source_lay,&source_rad,&source_vel);
    find_scatter(rad,vel,b,sc_dep,
		&sc_lay,&sc_rad,&sc_vel);
	if(fflag==0) 
	{ 
	  st_lay = source_lay;	  
	}
	find_ref(rad,ref1_depth,ref1_lay,
		     ref2_depth,ref2_lay,
		     ref3_depth,ref3_lay);
        distaz(90.-evla,evlo,90.-sclat,sclon,&dis_sc1,&az_sc1); 
      
        
     if(fflag==1){

        distaz(90.-evla,evlo,90.-sclat,sclon,&dis_sc1,&az_sc1); 
        distaz(90.-sclat,sclon,90.-stla,stlo,&dis_sc2,&az_sc2);
        if(source_lay <= sc_lay){
          lay_st = source_lay;
          lay_end = sc_lay;
          st_rad = source_rad;
          st_vel = source_vel;
          end_rad = sc_rad;
          end_vel = sc_vel;
          p_in = sc_rad/sc_vel;
        }
        else if(source_lay > sc_lay){
          lay_st = sc_lay;
          lay_end = source_lay;
	      st_rad = sc_rad;
          st_vel = sc_vel;
          end_rad = source_rad;
          end_vel = source_vel;
          p_in = sc_rad/sc_vel;
        }
        find_p(rad,vel,b,p_in,lay_st,lay_end,st_rad,st_vel,end_rad,end_vel,dis_sc1,&p_1,&delp,&trat,&pathp);
        if(abs(source_dep-sc_depth)<TINY && abs(evla-sclat)<TINY && abs(evlo-sclon)<TINY) trat=0;

     }
     else if(fflag==0){
        lay_st = 0;
        lay_end = sc_lay;
        st_rad = rad[lay_st];
        st_vel = vel[lay_st];
        end_rad = sc_rad;
        end_vel = sc_vel;
        p_in =  sc_rad/sc_vel;;
        distaz(90.-stla,stlo,90.-sclat,sclon,&dis_sc1,&az_sc1); 
        distaz(90.-sclat,sclon,90.-evla,evlo,&dis_sc2,&az_sc2);
        find_p(rad,vel,b,p_in,lay_st,lay_end,st_rad,st_vel,end_rad,end_vel,dis_sc1,&p_1,&delp,&trat,&pathp);
        if(abs(source_dep-sc_depth)<TINY && abs(evla-sclat)<TINY && abs(evlo-sclon)<TINY) trat=0;

     }
	if(delta<0.){/*calculate new delta only when input delata <0*/
		distaz(90.-evla,evlo,90.-stla,stlo,&delta,&az);

	}
	if(projection){
		evla=geographic_to_geocentric(evla)*RAD_PER_DEG; 
		evlo*=RAD_PER_DEG;
                stla=geographic_to_geocentric(stla)*RAD_PER_DEG; 
		stlo*=RAD_PER_DEG;
		sx=cos(evla)*cos(evlo);
		sy=cos(evla)*sin(evlo);
		sz=sin(evla);
                rx=cos(stla)*cos(stlo);
                ry=cos(stla)*sin(stlo);
                rz=sin(stla);
	}

	/*find the ray parameter which gives desired distance*/
	if(shooting[0]=='a'){
		p0=sc_rad*sin(st_val*RAD_PER_DEG)/sc_vel;
		p1=sc_rad*sin(end_val*RAD_PER_DEG)/sc_vel;
	}
	else{
		p0=st_val;
		p1=end_val;
	}

	find_deltt(p0,&del0,&tt);
	find_deltt(p1,&del1,&tt);
	if(del0>del1){
		p=p0; p0=p1; p1=p;
	}
   
	for(i=1;i<iteration;i++){
		p=(p0+p1)/2.;
		find_deltt(p,&del,&tt);
//        printf("%f %f\n",p,del/RAD_PER_DEG);

		/*ray parameter found*/
		if(fabs(del-dis_sc2)<=error) break;	
		if(del>dis_sc2)  p1=p;
		else p0=p;
	}
//         printf(" %f %f %f %f %f %f\n",p_1*RAD_PER_DEG, p*RAD_PER_DEG,dis_sc1,dis_sc2,trat,tt);	
	if(i==iteration){
		fprintf(stderr,"%d iterations. solution not found.\n",iteration);
		exit(-1);
	}


	inc_ang0 =asin(p*sc_vel/sc_rad);
	inc_val0=inc_ang0/RAD_PER_DEG;
	

	/* the following finds the ray path for the ray parameter*/
   		deltasum = 0.0;
   		ttsum = 0.0;
                pathsum = 0.0;
		ptr_ref1=ref1_type;
		ptr_ref2=ref2_type;
		ptr_ref3=ref3_type;

   		raycords[0].x = 0.0;
   		raycords[0].y = sc_rad/RSURF;
   		ncords = 1;
	
       		if(output){
			fprintf(outfp,"takeoff angle=%10.5f\n",inc_val0);
   			fprintf(outfp,"point_index depth incident_ang ttsum delsum lat lon\n");
			cmb_in=cmb_out=OFF;
		}
		if(moreoutput){
			if(projection) gcpt(sx,sy,sz,rx,ry,rz,delta,deltasum/RAD_PER_DEG,&lat,&lon);
			fprintf(outfp,"%4d %6.1f %7.2f %9.4f %8.4f %7.2f %7.2f\n",
				sc_lay,sc_dep,inc_val0,
				ttsum,deltasum/RAD_PER_DEG,lat,lon);
		}

		if(sc_dep==0.){
			if(inc_ang0 >=PI/2.){
				fprintf(stderr,"incident angle must be smaller than 90 deg for surface source!\n");
				exit(-1);
			}
			else
				downgoing(p,sc_rad,sc_vel,sc_lay);
		}	
		else if(inc_ang0 >PI/2.){/* go up */
			upgoing(p,sc_rad,sc_vel,sc_lay);
		}		
		else{
			downgoing(p,sc_rad,sc_vel,sc_lay);
		}

   		if(saveray)
			write_ray(rayfp,ncords,raycords);
   		if(output){
			fprintf(outfp,"bottoming depth =%7.1f\n", RSURF-bot_rad);
			if(cmb_in==ON)
				fprintf(outfp,"Entering radius %7.1f km: delta = %10.4f  time=%10.4f\n",
					output_rad,delcmb_in/RAD_PER_DEG,ttcmb_in);
			if(cmb_out==ON)
				fprintf(outfp,"Exiting radius %7.1f km: delta = %10.4f  time=%10.4f\n",
					output_rad,delcmb_out/RAD_PER_DEG,ttcmb_out);
			fprintf(outfp,"total delta= %10.4f\n",deltasum/RAD_PER_DEG);
			fprintf(outfp,"total time= %10.4f\n",ttsum);
			fprintf(outfp,"ray found: ");

        		write_time(outfp,dis_sc2,p_1,trat,sclat,sclon,sc_dep,delp);
		}
   
   write_time(outfp,dis_sc2,p_1,trat,sclat,sclon,sc_dep,delp);
   free(rad);
   free(vel);
   free(b);
   free(raycords);
}


void find_deltt(double p,double *del,double *tt)
{
	double inc_ang0;

	/* rays go from larger incident valle to smaller one.*/
	inc_ang0 =asin(p*sc_vel/sc_rad);
 	deltasum = 0.0;
   	ttsum = 0.0;
        pathsum = 0.0;
	ptr_ref1=ref1_type;
	ptr_ref2=ref2_type;
	ptr_ref3=ref3_type;

	if(inc_ang0 <PI/2.){/* go down */
		downgoing_fast(p,sc_rad,sc_vel,sc_lay);
	}		
	else{
		upgoing_fast(p,sc_rad,sc_vel,sc_lay);
	}
	*del=deltasum/RAD_PER_DEG;
	*tt=ttsum;
}

/* 
	a downgoing ray starting at r0 with v0 in lay0
*/
void downgoing_fast(double p, double r0, double v0, int lay0)
{
	int lay;
	double tt,delta,path,rnorm,rtop,rbot,vtop,vbot,cosD,sinD,tmp;
	double blay;
            
	rtop=r0;
	vtop=v0;
	for(lay=lay0; lay<=nlay-1;lay++){
      		rbot = rad[lay+1];
      		vbot = vel[lay+1];
      		blay = b[lay];

		if(lay==ref1_lay && *ptr_ref1=='t'){
			*ptr_ref1++;
			break;
		}
		if(lay==ref2_lay && *ptr_ref2=='t'){
			*ptr_ref2++;
			break;
		}
		if(lay==ref3_lay && *ptr_ref3=='t'){
			*ptr_ref2++;
			break;
		}

		if(rtop-rbot<TINY){	/* discontinuity */
	      		if(rbot/vbot-p<0.) break;
			vtop=vbot;
			continue;
		}
			
		/* travel within layer (lay) */
	      	if(rbot/vbot-p<0. || lay==nlay-1)  
		{
			straight_fast(p,rtop,vtop,lay);
			break;
		}

		else{ /* go down to next layer */
			trace(p,rtop,rbot,vtop,vbot,blay,&delta,&tt,&path);
      			ttsum += tt;
      			deltasum += delta;
                 pathsum += path;
      		}
		rtop=rbot;
		vtop=vbot;
	}
	upgoing_fast(p,rad[lay],vel[lay],lay-1);
}

/*
	an upwards-going ray starting at r0 with v0 in lay0
*/
void upgoing_fast(double p, double r0, double v0, int lay0)
{
	int lay;
	double delta,path,tt,rnorm,rtop,rbot,cosD,sinD;
	double blay,vtop,vbot;
    int lay_st;
	rbot=r0;
	vbot=v0;

	if(fflag==1){
	  lay_st = 0;
	}
	else if(fflag==0){
	  lay_st = st_lay;
	}

   	for(lay=lay0;lay>=lay_st;lay--)
      	{
      		rtop = rad[lay];
      		vtop = vel[lay];
      		blay = b[lay];


		if(rtop-rbot<TINY){	/* discontinuity */
			goto checkbounce;
		}

		trace(p,rtop,rbot,vtop,vbot,blay,&delta,&tt,&path);
      		ttsum += tt;
      		deltasum += delta;
                pathsum += path;

checkbounce:	if(lay==ref1_lay && *ptr_ref1=='u'){
                        ptr_ref1++;
                        downgoing_fast(p,rad[lay],vel[lay],lay);
			break;
                }
                if(lay==ref2_lay && *ptr_ref2=='u'){
                        ptr_ref2++;
                        downgoing_fast(p,rad[lay],vel[lay],lay);
                        break;
                }
                if(lay==ref3_lay && *ptr_ref3=='u'){
                        ptr_ref3++;
                        downgoing_fast(p,rad[lay],vel[lay],lay);
                        break;
                }

		rbot=rtop;
		vbot=vtop;
      	}
}

/* 
	a ray traveling within one layer
	starting at r0 with v0 in lay0 to the top of lay0
*/
void straight_fast(double p, double r0, double v0, int lay0)
{
	double rtop,rbot,vtop,vbot,blay,rnorm;
	double cosD,sinD,delta,tt,path;
	
	
	rtop=rad[lay0];
	vtop=vel[lay0];
	blay=b[lay0];
	delta=(acos(p*v0/r0)+acos(p*vtop/rtop))/(1.-blay);
	tt=( sqrt(abs((r0*r0)/(v0*v0)-p*p)) + sqrt(abs((rtop*rtop)/(vtop*vtop)-p*p) ))/(1.-blay);
        path = sqrt(abs(r0*r0-v0*v0*p*p))+sqrt(abs(rtop*rtop-vtop*vtop*p*p))/(1.-blay);
       	ttsum += tt;
   	deltasum += delta;
        pathsum += path;
}


/* 
	a downgoing ray starting at r0 with v0 in lay0
*/
void downgoing(double p, double r0, double v0, int lay0)
{
	int lay;
	double tt,delta,path,rnorm,rtop,rbot,vtop,vbot,cosD,sinD,tmp;
	double blay;

	rtop=r0;
	vtop=v0;
	for(lay=lay0; lay<=nlay-1;lay++){
      		rbot = rad[lay+1];
      		vbot = vel[lay+1];
      		blay = b[lay];

		if(rbot <output_rad && cmb_in==OFF){	/*save delta entering CMB*/
			delcmb_in=deltasum;
			ttcmb_in=ttsum;
			cmb_in=ON;
		}

		if(lay==ref1_lay && *ptr_ref1=='t'){
			bot_rad=rad[lay];
	 		fprintf(stderr,"topside reflection on %lf km\n",RSURF-bot_rad);
			*ptr_ref1++;
			break;
		}
		if(lay==ref2_lay && *ptr_ref2=='t'){
			bot_rad=rad[lay];
	 		fprintf(stderr,"topside reflection on %lf km\n",RSURF-bot_rad);
			*ptr_ref2++;
			break;
		}
		if(lay==ref3_lay && *ptr_ref3=='t'){
			bot_rad=rad[lay];
	 		fprintf(stderr,"topside reflection on %lf km\n",RSURF-bot_rad);
			*ptr_ref2++;
			break;
		}

		if(rtop-rbot<TINY){	/* discontinuity */
	      		if(rbot/vbot-p<0.) break;
			if(moreoutput)
				write_moreoutput(lay+1);
			vtop=vbot;
			continue;
		}
			
		/* travel within layer (lay) */
	      	if(rbot/vbot-p<0. || lay==nlay-1)  
		{
			straight(p,rtop,vtop,lay);
			break;
		}

		else{ /* go down to next layer */
			trace(p,rtop,rbot,vtop,vbot,blay,&delta,&tt,&path);
      			ttsum += tt;
      			deltasum += delta;
                        pathsum += path;
			if(moreoutput)
				write_moreoutput(lay+1);
      			cosD = cos(deltasum);
      			sinD = sin(deltasum);
      			rnorm = rbot/RSURF;
      			raycords[ncords].x = rnorm*sinD;
      			raycords[ncords].y = rnorm*cosD;
      			ncords++;

      		}
		rtop=rbot;
		vtop=vbot;
	}
	upgoing(p,rad[lay],vel[lay],lay-1);
}

/*
	an upwards-going ray starting at r0 with v0 in lay0
*/
void upgoing(double p, double r0, double v0, int lay0)
{
	int lay;
	double delta,tt,path,rnorm,rtop,rbot,cosD,sinD;
	double blay,vtop,vbot;

	rbot=r0;
	vbot=v0;
   	for(lay=lay0;lay>=0;lay--)
      	{
      		rtop = rad[lay];
      		vtop = vel[lay];
      		blay = b[lay];

		if(rtop>output_rad && cmb_in==ON && cmb_out==OFF){
			delcmb_out=deltasum;
			ttcmb_out=ttsum;
			cmb_out=ON;
		}

		if(rtop-rbot<TINY){	/* discontinuity */
			if(moreoutput)
				write_moreoutput(lay);
			goto checkbounce;
		}

		trace(p,rtop,rbot,vtop,vbot,blay,&delta,&tt,&path);
      		ttsum += tt;
      		deltasum += delta;
                pathsum += path;
		if(moreoutput)
			write_moreoutput(lay);

      		cosD = cos(deltasum);
      		sinD = sin(deltasum);
      		rnorm = rtop/RSURF;
      		raycords[ncords].x = rnorm*sinD;
      		raycords[ncords].y = rnorm*cosD;
      		ncords++;

checkbounce:	if(lay==ref1_lay && *ptr_ref1=='u'){
                        fprintf(stderr,"underside reflection at %lf km\n",
                                RSURF-rad[lay]);
                        ptr_ref1++;
                        downgoing(p,rad[lay],vel[lay],lay);
                        break;
                }
                if(lay==ref2_lay && *ptr_ref2=='u'){
                        fprintf(stderr,"underside reflection at %lf km\n",
                                RSURF-rad[lay]);
                        ptr_ref2++;
                        downgoing(p,rad[lay],vel[lay],lay);
                        break;
                }
                if(lay==ref3_lay && *ptr_ref3=='u'){
                        fprintf(stderr,"underside reflection at %lf km\n",
                                RSURF-rad[lay]);
                        ptr_ref3++;
                        downgoing(p,rad[lay],vel[lay],lay);
                        break;
                }

		rbot=rtop;
		vbot=vtop;
      	}
}

/* 
	a ray traveling within one layer
	starting at r0 with v0 in lay0 to the top of lay0
*/
void straight(double p, double r0, double v0, int lay0)
{
	double rtop,rbot,vtop,vbot,blay,rnorm;
	double cosD,sinD,delta,tt,path;
	
	
	rtop=rad[lay0];
	vtop=vel[lay0];
	blay=b[lay0];
	delta=(acos(p*v0/r0)+acos(p*vtop/rtop))/(1.-blay);
	tt=( sqrt(abs((r0*r0)/(v0*v0)-p*p)) + sqrt(abs((rtop*rtop)/(vtop*vtop)-p*p)) )/(1.-blay);
        path = sqrt(abs(r0*r0-v0*v0*p*p))+sqrt(abs(rtop*rtop-vtop*vtop*p*p))/(1.-blay);
       	ttsum += tt;
   	deltasum += delta;
        pathsum += path;
	if(moreoutput)
		write_moreoutput(lay0);

       	cosD = cos(deltasum);
       	sinD = sin(deltasum);

        rnorm = rtop/RSURF;
        raycords[ncords].x = rnorm*sinD;
        raycords[ncords].y = rnorm*cosD;
       	ncords++;

	if(p==0.) {
		bot_rad=0.;
		return;
	}
	bot_rad=exp( (log(v0*p)-blay*log(r0)) /(1.-blay) );

}

/*
	find dt,ddelta within a lay with r=ar^b approximation
*/
void trace(p,rtop,rbot,vtop,vbot,blay,delta,tt,path)
double p,rtop,rbot,vtop,vbot,blay,*delta,*tt,*path;
{
	double a;

	if(blay!=1.0){
		*delta=acos(p*vtop/rtop)-acos(p*vbot/rbot);
		*delta /=1.-blay;
		*tt=sqrt(abs(rtop/vtop*rtop/vtop-p*p))
			-sqrt(abs(rbot/vbot*rbot/vbot-p*p));
                if(rtop==rbot) *tt=(0.05/vtop-0.05/vbot) ;
		*tt /=1.-blay;
                *path = sqrt(abs(rtop*rtop-vtop*vtop*p*p))-sqrt(abs(rbot*rbot-vbot*vbot*p*p));
               // *path /=1.-blay;
	}
	else{	/* v=ar */
		fprintf(stderr,"r=%lf b=1\n",rtop);
		a=vtop/rtop;
		*tt=log(rbot/rtop)/a/sqrt(abs(1.-p*p*a*a));
		*delta=p*a*a*(*tt);
                *path = fabs(rbot-rtop)/sqrt(abs(1.-p*p*a*a));
	}
}

/*
	read in model velocity
*/
void read_vel(FILE *fp,int npt,double *rad,double *vel)
{
	int i;
	double dep;

	for(i=0; i<=npt-1; i++){
   		if(fscanf(fp,"%lf%lf%*f%*f",&dep,&vel[i]) == EOF){
			fprintf(stderr,"Warning:  End  of model. %d points read.\n",i);
			break;
		}
		rad[i]=RSURF-dep;
        }
}

/*
	find b value of each layer
*/
void find_b(double *r,double *v,double *b)
{
	int i;

	for(i=0; i<=nlay-2; i++){
		if(r[i]-r[i+1]<TINY){
//			fprintf(stderr,"discontinuity at %7.1lf km\n",
//				RSURF-r[i]);
			continue;
		}
		b[i]=(log(v[i+1])-log(v[i])) /(log(r[i+1])-log(r[i]));
	}
	b[i]=0.;
        
}

/*
	find some source parameters
*/
void find_source(double *rad,double *vel,double *b,double source_dep,
	int *source_lay,double *source_rad,double *source_vel)
{
	int i;

        if(source_dep >6371.){
                fprintf(stderr,"source too deep\n");
                exit(-1);
        }
        *source_rad=RSURF-source_dep;
	for(i=0;i<=nlay-1;i++)
		if(rad[i]>= *source_rad && rad[i+1]< *source_rad) break;

	/* if source at a discontinuity */
	if(rad[i]-rad[i+1]<TINY) i++;
	*source_lay=i;
        *source_vel=vel[i]*exp(b[i]*log(*source_rad/rad[i]));
}

void find_scatter(double *rad,double *vel,double *b,double sc_dep,
	int *sc_lay,double *sc_rad,double *sc_vel)
{
        int i;

        if(sc_dep >2891.){
                fprintf(stderr,"scatter too deep\n");
                exit(-1);
        }
        *sc_rad=RSURF-sc_dep;
	for(i=0;i<=nlay-1;i++)
		if(rad[i]>= *sc_rad && rad[i+1]< *sc_rad) break;

	/* if scatter is at a discontinuity */
	if(rad[i]-rad[i+1]<TINY) i++;
	*sc_lay=i;
        *sc_vel=vel[i]*exp(b[i]*log(*sc_rad/rad[i]));
}
// Find the distance and delta between the source and scatters, Warning the scattering depth should be deeper than like 800!!!!!
void find_p(double *rad,double *vel,double *b,double p_in,int lay_st,int lay_end,double st_rad,double st_vel,double end_rad,double end_vel,float dis_sc1,double *p_1,double *delp,double *trat,double *pathp)
{     

       int k = 0,i=0; 
       double rtop,rbot,vtop,vbot,tt,delseg,path;
       double blay; 
       double p_sta,p_end,p_test;
       double sumdelp,sumtrat,sumpathp;
       p_sta = p_in;
       p_end = 0;

       for(i=1;i<iteration;i++){
       p_test = (p_sta+p_end)/2.0;
       sumdelp = 0;
       sumtrat = 0;
       sumpathp = 0;

         rtop= st_rad;
         vtop= st_vel;

       for(k=lay_st;k<=lay_end;k++){

      	rbot = rad[k+1];
      	vbot = vel[k+1];
      	blay = b[k];
        if(k==lay_end) {rbot = end_rad;vbot=end_vel;}
        if(rtop-rbot<TINY){	
	      		if(rbot/vbot-p<0.) break;
			vtop=vbot;
			continue;
		}
    if(rbot/vbot-p_test<0.0)  
	{
		float r1,v1,b1;
        r1=rad[k];
	    v1=vel[k];
	    b1=b[k];
	    delseg=(acos(p*vtop/rtop)+acos(p_test*v1/r1))/(1.-b1);
	    tt=( sqrt(abs((rtop*rtop)/(vtop*vtop)-p_test*p_test)) + sqrt(abs((r1*rtop)/(v1*v1)-p_test*p_test)) )/(1.-b1);
        path = sqrt(abs(rtop*rtop-vtop*vtop*p_test*p_test))+sqrt(abs(r1*r1-v1*v1*p_test*p_test))/(1.-b1);
       	sumtrat += tt;
   	    sumdelp += delseg;
        sumpathp += path;
		break;
	}


		else{ /* go down to next layer */
		trace(p_test,rtop,rbot,vtop,vbot,blay,&delseg,&tt,&path);   
      		sumtrat += tt;
      		sumdelp += delseg;
                sumpathp += path;
//                fprintf(fp,"%f %f %f %f \n",p_test,sumdelp/RAD_PER_DEG,dis_sc1,delseg);
      		}
         
		rtop=rbot;
		vtop=vbot;
      }/*ray parameter found*/

            
            if(fabs(sumdelp/RAD_PER_DEG-dis_sc1)<=error) { 
            *p_1 =p_test;
            *trat = sumtrat;
            *delp = sumdelp/RAD_PER_DEG;
            *pathp = sumpathp; 
//            printf("%f %f %f\n",sumtrat,p_test*RAD_PER_DEG,*delp);
             break;}
           	
		if(sumdelp/RAD_PER_DEG>dis_sc1) {
		  p_sta=p_test;
		}
		else if(sumdelp/RAD_PER_DEG<dis_sc1) {
		  p_end=p_test;
		}
		
	  }       
       
}


/*
	find reflection layers
*/ 
void find_ref(double *rad,
	float ref1_depth,int lay1, 
	float ref2_depth,int lay2,
	float ref3_depth,int lay3)
{
	int i;
	float depth;

	if(lay1 !=UNDEF) ref1_lay=lay1;
	if(lay2 !=UNDEF) ref2_lay=lay2;
	if(lay3 !=UNDEF) ref3_lay=lay3;

	if(ref1_depth !=UNDEF){
		depth=RSURF-ref1_depth;
	        for(i=0;i<=nlay-1;i++)
	        if(fabs(rad[i]-depth) <TINY){
       			ref1_lay=i;
       		        break;
               	}	
	}
        if(ref2_depth !=UNDEF){
                depth=RSURF-ref2_depth;
                for(i=0;i<=nlay-1;i++)
                if(fabs(rad[i]-depth) <TINY){
                        ref2_lay=i;
                        break;
                }
        }
        if(ref3_depth !=UNDEF){
                depth=RSURF-ref3_depth;
                for(i=0;i<=nlay-1;i++)
                if(fabs(rad[i]-depth) <TINY){
                        ref3_lay=i;
                        break;
                }
        }
}

/*
	output coordinats of a ray
*/
void write_ray(FILE *fp,int np,struct cords *raycords)
{
   int i;

   fprintf(fp,"MOVE %8.4f,%8.4f\n", raycords[0].x,raycords[0].y);
   for(i=1; i<=np-1; i++)
      fprintf(fp,"DRAW %8.4f,%8.4f\n", raycords[i].x,raycords[i].y);
}

/*
	output travel time of a ray
*/
void write_time(FILE *fp,double dis_sc2,double p_1,double trat,float sclat,float sclon,double sc_dep,double delp)
{
/*
	fprintf(fp, "DRAW %10.4f, %11.4f !%10.5f\n",deltasum/RAD_PER_DEG,ttsum,inc_val0);
*/
	float del_found,t_scatter;

	del_found=deltasum/RAD_PER_DEG;
//	t=ttsum+(dis_sc2-del_found)*p*RAD_PER_DEG;
        t_scatter = ttsum + trat;

        printf("%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f\n",
		sclat,sclon,sc_dep,t_scatter,trat,ttsum,del_found,delp);
}

/*
	write out:
	the layer the ray is traveling
	radius of the top of the layer
	incident angle at the top of  the layer
	accumulative time  from source to the radius
	accumulative delta  from source to the radius
*/
void write_moreoutput(int lay)  
{
	float lat,lon;

	if(projection) gcpt(sx,sy,sz,rx,ry,rz,delta,deltasum/RAD_PER_DEG,&lat,&lon);
	fprintf(outfp,"%4d %6.1f %7.2f %9.4f %8.4f %7.2f %7.2f\n",
		lay,RSURF-rad[lay],asin(p*vel[lay]/rad[lay])/RAD_PER_DEG,
		ttsum,deltasum/RAD_PER_DEG,lat,lon);
}
