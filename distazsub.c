#include <stdio.h>
#include <math.h>
#define RAD_PER_DEG	0.0174532925199432955

/*
#define B2A_SQ	0.99327733 
#define B2A_SQ	0.9931177
*/
#define B2A_SQ	0.993305521
/* 
  convert geographic latitude (degree) to geocentric latitude (degree)
*/
float geographic_to_geocentric(float lat)
{
	return(atan(B2A_SQ*tan(RAD_PER_DEG*lat))/RAD_PER_DEG);
}

/*
  convert geocentric latitude (degree) to geographic latitude (degree)
*/
float geocentric_to_geographic(float lat)
{
	return(atan(tan(RAD_PER_DEG*lat)/B2A_SQ)/RAD_PER_DEG);
}

/*
  calculate distance, azimuth between source (c0,l0) and receiver (c1,l1).
  colatitude:[0,180] (geographic), longitude [0,360].
  distance <=180.
  azimuth [0,360), measured at source location clockwise from north.
*/
void distaz(float c0,float l0,float c1,float l1,float *dist,float *az)
{
	double colat0,lon0,colat1,lon1,del,azimuth;
	double cosdel,cosaz,sinaz,tmp;
	double eps=1.0e-8;
	float geographic_to_geocentric(float lat);

	/*convert geographic to geocentric*/
	c0=90.-geographic_to_geocentric(90.-c0);
	c1=90.-geographic_to_geocentric(90.-c1);

	colat0 =(double)c0*RAD_PER_DEG;
	colat1 =(double)c1*RAD_PER_DEG;
	lon0 =(double)l0*RAD_PER_DEG;
	lon1 =(double)l1*RAD_PER_DEG;

	/*calculate distance*/
	cosdel=cos(colat0)*cos(colat1)+
		sin(colat0)*sin(colat1)*cos(lon1-lon0);
	if(cosdel >1.0) cosdel -=eps;
	if(cosdel <-1.0) cosdel +=eps;
	del=acos(cosdel);
	*dist =del/RAD_PER_DEG;

	/*calculate azimuth*/
	tmp=sin(colat0)*sin(del);
	if (tmp <=eps){ 
		/*special case: source at pole or del=0 or 180*/
		if(c0 <=eps) *az=180.;
		if(c0 >=(180.-eps)) *az=0.;
		if((*dist) <=eps) *az=-999;
		if((*dist) >=(180.-eps)) *az=-999;
		return;
	}
	cosaz=(cos(colat1)-cos(colat0)*cos(del))/tmp;
	sinaz=sin(colat1)*sin(lon1-lon0)/sin(del);
	azimuth=atan2(sinaz,cosaz)/RAD_PER_DEG;
	if(azimuth <0.) azimuth +=360.;
	*az =azimuth;
}

/* find a point dis1 degree away from the first point on the 
   great cirlce of two points.  Use OA=a1*OA1+a2*OA2 
   x1,y1,z1; x2,y2,z2:  cartesian coordinates of the points. 
			They should be in geocentric coordinates.
   dis:	  distance of the two point in degree, typically calculated by distaz above
   dis1:  distance from the first point. positive means towards the second point
   lat,lon:  latitude and longitude of the point dis1 degree away from the first point.
             it is in geographic coordinates.		
*/
void gcpt(float x1,float y1,float z1,float x2,float y2,float z2,
	float dis,float dis1,float *lat,float *lon)
{
	float sin2,a1,a2;
	float x,y,z;
	float eps=1.0e-7;

	dis *=RAD_PER_DEG;
	dis1 *=RAD_PER_DEG;
	sin2=sin(dis);
	sin2 *=sin2;
	a1=(cos(dis1)-cos(dis-dis1)*cos(dis))/sin2;
	a2=(cos(dis-dis1)-cos(dis)*cos(dis1))/sin2;
	x=a1*x1+a2*x2;
	y=a1*y1+a2*y2;
	z=a1*z1+a2*z2;
/*
	if(z>=1.0) z-=eps; if(z<=-1.0) z+=eps;
*/
	if(z>=1.0) z=1.0-eps; if(z<=-1.0) z=-1.0+eps;
	*lat=90.-acos( (double)z )/RAD_PER_DEG;
	*lat=geocentric_to_geographic((*lat));
	*lon=atan2((double)y,(double)x)/RAD_PER_DEG;
}

