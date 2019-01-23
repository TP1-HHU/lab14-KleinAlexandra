#include "EM.hxx"
#include <cmath>
#include <iostream>
//----------------------------------
using namespace std;
//----------------------------------
class particle{
public:
	particle(){x=0; y=0; z=0;  // location at t = 0
		         px=0; py=0; pz=0;  // momentum at t = -dt
						 gam=0;};  // gamma at t= -dt};
	// Print particle position, momentum u and gamma-factor
	void print_info(const double t) const;
	
  double get_x() const{return x;}
  void boris_push(double Ex, double Ey, double Ez, double Bx, double By, double Bz, double dt);
private:
	double x,y,z;
	double px,py,pz;
  double gam;
	// Calculate gamma-factor from current momentum

  void gamma();
  
  
};
//----------------------------------
//----------------------------------
int main(){
  const double FWHM = 10;
  const double a0 = 10;
  const double T0 = 200;
  const double Tend = 800;
  const double dt = 0.001;
  const int N = int(Tend/dt + 0.5);

  double Ex,Ey,Ez, Bx, By, Bz;
  double t=0;

  EM em(a0, FWHM,T0);
  particle p;

  for(int i=0; i<N; i++){
	  em.fields(Ex,Ey,Ez,Bx,By,Bz, t + dt*0.5, p.get_x());
	  p.boris_push(Ex, Ey, 0, Bx, 0, Bz, dt);

	  if (i%5 == 0 ) p.print_info(t);
	  t += dt;

  }


  return 0;
}
//----------------------------------
void particle::gamma()
{
    gam = sqrt(1+(px*px+py*py+pz*pz));
    
}
//----------------------------------
void particle::boris_push(double Ex, double Ey, double Ez, double Bx, double By, double Bz, double dt)
{
    double pmx = px + dt/2*Ex;
    double pmy = py + dt/2*Ey;
    double pmz = pz + dt/2*Ez;
    
    gamma();
    
    double tx = Bx*dt/2/gam;
    double ty = By*dt/2/gam;
    double tz = Bz*dt/2/gam;   
    
    double prix = pmx + pmy*tz - ty*pmz;  
    double priy = pmy + pmz*tx - tz*pmx;
    double priz = pmz + pmx*ty - tx*pmy;
    
    double sx = 2*tx/(1+(tx*tx+ty*ty+tz*tz));
    double sy = 2*ty/(1+(tx*tx+ty*ty+tz*tz));
    double sz = 2*tz/(1+(tx*tx+ty*ty+tz*tz));
    
    double ppx = pmx + priy*sz - priz*sy;
    double ppy = pmy + priz*sx - prix*sz;
    double ppz = pmz + prix*sy - priy*sx;
    
    px = ppx + 0.5*dt*Ex;
    py = ppy + 0.5*dt*Ey;
    pz = ppz + 0.5*dt*Ez;

    gamma();
    
    x = x + px/gam*dt;
    y = y + py/gam*dt;
    z = z + pz/gam*dt;
    
}


//--------------------------------------
void particle::print_info(const double t) const
{
	cout << t << "\t" << x << "\t" << y << "\t" << z
			 << "\t" << px << "\t" << py << "\t" << pz << "\t" << gam << endl;
}
