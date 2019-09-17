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
