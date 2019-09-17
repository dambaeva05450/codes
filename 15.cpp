#include <iostream>
#include <vector>
#include<random>
#include <algorithm>
#include <ctime>
#include <fstream>
using namespace std;
const int n=3;
typedef vector<double> person;
double gen_chislo(double left, double right){
	std::random_device rd;   // non-deterministic generator
	std::mt19937 gen(rd());
	uniform_real_distribution<> dist(left,right);
	return dist(gen);
}
double searchd(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double a=y1-y2, b=x2-x1, c=x1*y2-x2*y1;
	return abs(a*x3+b*y3+c)/sqrt(a*a+b*b);
}

vector<double> pr1(double x)
{
	vector<double>v(3);
	 v[0]=4*x/3;
	 v[1]=3*x/4;
	 v[2]=(2*x+14)/5;
	sort(v.begin(), v.end());

	return v;
}
double gen_rad(double x, double y, vector<double>tr)
{
	vector<double> t(3);
	t[0]=searchd(tr[0],tr[1],tr[2],tr[3],x,y);
	t[1]=searchd(tr[2],tr[3],tr[4],tr[5],x,y);
	t[2]=searchd(tr[4],tr[5],tr[0],tr[1],x,y);
	return min(t[0],min(t[1],t[2]));
}
double ss2(person p, double a)
{
		double s2=0;
		person c;
	for (int i=0; i<n-1; i++)
	{
		for (int j=i+1; j<n; j++)
		{
			
			
			
					//расстояние между центрами двух окружностей
			
			double ro=sqrt((p[i*3]-p[j*3])*(p[i*3]-p[j*3])+(p[3*i+1]-p[3*j+1])*(p[3*i+1]-p[3*j+1]));
			
			//если окружности пересекаются
			if(p[3*i+2]+p[3*j+2]>ro && ro+p[3*j+2]>p[3*i+2] && p[3*i+2]+ro>p[3*j+2])
			{			
				double f1=2*acos((p[3*i+2]*p[3*i+2]-p[3*j+2]*p[3*j+2]+ro*ro)/(2*p[3*i+2]*ro));
				double f2=2*acos((p[3*j+2]*p[3*j+2]-p[3*i+2]*p[3*i+2]+ro*ro)/(2*p[3*j+2]*ro));
				//cout<<f1<<' '<<f2<<' '<<i<<' '<<j<<endl;
				double l=p[3*i+2]*p[3*i+2]*(f1-sin(f1))/2+p[3*j+2]*p[3*j+2]*(f2-sin(f2))/2;
				s2+=l;
			}
		}
	}
	
	
	
		
		return a*s2;	
}
double ss3( person p)
{
	double s3=0;
		for (int i=2; i<p.size(); i+=3)
	{
		s3+=acos(-1)*p[i]*p[i];
	}
	return s3;
}
double fin( person p, double a)
{
	
	return ss3(p)-ss2(p,a);
}

vector<double> exper(int np, int max_iter, double cr, double w)
{

vector<double>trien={0,0,3,4,8,6};

 double a=1000;
		
				vector<person>p;
				//сгенерировали 10 особей (каждая это есть система из трех окружностей)
				for(int i=0;i<np;++i)
				{
					person hh(3*n);
					for(int k=0; k<n; k++)
						{
							hh[k*3]=gen_chislo(0,8);
							vector<double> b=pr1(hh[k*3]);
							hh[k*3+1]=gen_chislo(b[0], b[1]);
							hh[k*3+2]=gen_rad(hh[k*3],hh[k*3+1],trien);
											
						}
								
					p.push_back(hh);
				}
						
									//диф эволюция
									
				for(int iter=0;iter<max_iter;++iter)
					{
										
						vector<person> x(3, person(n*3));
						for(int j=0; j<np; ++j)
							{
								for(int k=0; k<3; ++k)
								{
									int help_num = rand()%np;
									while(help_num==j)
									{
										help_num = rand()%np;
									}
									for(int l=0; l<n*3; ++l)
									x[k][l]=p[help_num][l];
								
								}
								person x_cl(n*3), x_s(n*3);
								for(int k=0; k<n; k++)
								{
									x_cl[k*3]=x[2][k*3]+w*(x[0][k*3]-x[1][k*3]);
												
									if(x_cl[k*3]>8 || x_cl[k*3]<0)
										x_cl[k*3]=gen_chislo(0,8);
									
									x_cl[k*3+1]=x[2][k*3+1]+w*(x[0][k*3+1]-x[1][k*3+1]);
									
									vector<double> b=pr1(x_cl[k*3]);
									
									if(x_cl[k*3+1]>b[1] || x_cl[k*3+1]<b[0])
										x_cl[k*3+1]=gen_chislo(b[0],b[1]);
									
									x_cl[k*3+2]=gen_rad(x_cl[3*k],x_cl[3*k+1],trien);
												
								}
											
								for(int i=0; i<n-1; ++i)
								{
									double u_i=gen_chislo(0,1);
									if(u_i<=cr)
									{
										x_s[3*i]=x_cl[3*i];
										x_s[3*i+1]=x_cl[3*i+1];
										x_s[3*i+2]=x_cl[3*i+2];
													
									}
													
									else
									{
										x_s[3*i]=p[j][3*i];
										x_s[3*i+1]=p[j][3*i+1];
										x_s[3*i+2]=p[j][3*i+2];
									}
													
								}
											
								x_s[n*3-3]=x_cl[n*3-3];
								x_s[n*3-2]=x_cl[n*3-2];
								x_s[n*3-1]=x_cl[n*3-1];
											
										
								if(fin(x_s,a)>fin(p[j],a))
									p[j]=x_s;
							}
					}
							
							person best_x=p[0],  worst_x=p[0], f_range(np), rez(4);
							double best_f=fin(p[0],a), worst_f=fin(p[0],a), fav=fin(p[0],a),s=0;
							f_range[0]=fin(p[0],a);
							for(int j=1;j<np;++j)
								{
									f_range[j]=fin(p[j],a);
									fav+=fin(p[j],a);
									if(fin(p[j],a)>fin(best_x,a))
									{
											
										best_x=p[j];
										best_f=fin(p[j],a);
											
									}
									if(fin(p[j],a)<fin(worst_x,a))
									{
											
										worst_x=p[j];
										worst_f=fin(p[j],a);
											
									}
									
										
								}
								fav/=np;
								for(auto i:f_range)
								{
									s+=(i-fav)*(i-fav);
								}
			rez[0]=best_f;
			rez[1]=worst_f;
			rez[2]=fav;
			rez[3]=s;
			
return rez;

}
						
			//4.51342 4.04039 0.524261 3.4295 3.43025 0.685449 2.57288 2.5882 0.505386 3.14192 0



int main(){
	ofstream file;
	file.open("result.txt");
	int max_iter=1000,np;
	for (int ind_np=5; ind_np<=10; ind_np++)
	{
		np=ind_np*6;
		for(double w=0.5; w<=0.9; w++)
		{
			for (double cr=0.9; cr>=0.1; cr-=0.2)
			{
				
				double fmx=0, fmn=64, fkv, fav=0;
				vector<double>f(4,0), c;
				for(int exp_i=0;exp_i<10;exp_i++)
				{
					
					c=exper(np,max_iter,cr,w);
					for(int i=0; i<4; i++)
					f[i]+=c[i]/10;
					
				}
				fav/=10;
				for(int i=0; i<4; i++)
					file<<f[i]<<"         ";
				file<<endl;
				
			}
		}
	}
		
	
				file.close();
}
//2.30931 2.30209 0.456088 4.0025 3.76737 0.588315 3.05926 3.04462 0.600143
//2.57184 2.55896 0.504069 4.43827 3.99919 0.534918 3.39181 3.40142 0.672594
//3.4865 3.45958 0.675767 4.58262 4.08551 0.508381 2.61413 2.61247 0.5215
//2.87938 2.68704 0.422008 4.06868 3.79063 0.591294 3.50174 3.07698 0.360543
//5.20128 4.42696 0.420796 3.41731 3.41258 0.679679 4.44995 3.93449 0.477626
//2.61993 2.62186 0.522828 3.44543 3.44796 0.677984 4.44395 4.00145 0.534786 3.17584
//4.49471 4.03259 0.52486 3.44916 3.44288 0.684086 2.55623 2.53389 0.493375 3.10034
//4.02623 3.52296 0.402623 5.35795 4.48082 0.369886 4.66667 4.08333 0.466667
//5.33657 4.46853 0.37288 4.65561 4.07366 0.465561 4.01352 3.51183 0.401352 1.61253
//3.29526 3.29477 0.658664 2.54435 2.64612 0.447805 5.15336 4.36318 0.39853 2.45501
//4.5101 4.0503 0.514132 3.42312 3.4385 0.675403 2.52347 2.52204 0.503553 3.06012
//2.52408 2.52575 0.503817 4.64637 4.1293 0.491394 3.42381 3.42579 0.683577 3.02403
//2.62748 2.64851 0.512875 3.49013 3.47004 0.674089 4.5244 4.05794 0.51235 3.07857
//3.45127 3.44433 0.683523 2.50395 2.52012 0.491084 4.4757 4.01991 0.529577 3.10647
//3.414 3.41944 0.679537 4.53427 4.05052 0.519856 2.59755 2.5856 0.509948 3.11667
//4.30262 3.93778 0.54155 3.4592 3.44766 0.682604 2.62976 2.62978 0.52594 3.12061
//2.64575 2.63623 0.521532 4.50665 4.01776 0.510218 3.472 3.45549 0.680865 3.1287
//3.37765 3.37977 0.674256 4.4297 3.99354 0.536974 2.40881 2.39228 0.468539 3.02375
//2.62716 2.63612 0.520056 4.32223 3.93259 0.55273 3.41539 3.39524 0.66696 3.11019 0.096762
//4.37376 3.97066 0.537441 2.619 2.61744 0.522556 3.4667 3.43838 0.670685 3.11449 0.0639402
//2.60447 2.61218 0.516268 4.54286 4.05184 0.515758 3.47252 3.46036 0.676531 3.11091 0
//4.50753 4.03687 0.524983 3.45279 3.44991 0.678906 2.58312 2.57181 0.507577 3.12323 0
//2.58296 2.58159 0.515497 3.42504 3.43643 0.678173 4.50107 4.04103 0.519386 3.1272 0
//2.56668 2.56877 0.51208 3.43056 3.43267 0.684845 4.49328 4.01724 0.517824 3.13965 0
//2.57419 2.57179 0.512916 3.44503 3.44108 0.684226 4.50594 4.03503 0.524458 3.1614 0
//3.40933 3.40551 0.678811 4.51372 4.04461 0.520763 2.56327 2.56892 0.509268 3.11436 0
//3.43572 3.43487 0.686465 4.48295 4.02824 0.524529 2.5847 2.58589 0.516225 3.18197 0
//3.42719 3.43077 0.683292 4.49253 4.029 0.527387 2.57002 2.56635 0.511067 3.16112 0