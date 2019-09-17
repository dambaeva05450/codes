#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
using namespace std;
int const n=3;

typedef vector <double> person;
vector<double>coor(n*2);
vector<double> arr(n);
vector<double> quickSortR(int b, int e) {
vector<double>qq(n*3);
	int l=b;
	int r=e;
	double piv=arr[(l+r)/2];
	while (l <= r)
  {
    while (arr[l] < piv)
      l++;
    while (arr[r] > piv)
      r--;
    if (l <= r)
	{
		swap (arr[l], arr[r]);
		
		swap(coor[l*2],coor[r*2]);
		swap(coor[l*2+1],coor[r*2+1]);
		l++; r--;
	}
      
  }
  
  if (b < r)
    quickSortR (b, r);
  if (e > l)
    quickSortR (l, e);

for(int i=0; i<n; i++)
{
	qq[i*3]=coor[2*i];
	qq[i*3+1]=coor[2*i+1];
	qq[i*3+2]=arr[i];
}
return(qq);
}
double gen_chislo(double left, double right){
	random_device rd;   // non-deterministic generator
	mt19937 gen(rd());
	uniform_real_distribution<> dist(left,right);
	return dist(gen);
}
vector<double> gen_bounderies (vector<double> tr)
{
	vector<double>bound;
	double xmn=min(tr[0],min(tr[2],tr[4]));
	double xmx=max(tr[0],max(tr[2],tr[4]));	
	double ymn=min(tr[1],min(tr[3],tr[5]));
	double ymx=max(tr[1],max(tr[3],tr[5]));
	double mx=max(xmx,ymx);
	double mn=min(xmn,ymn);
	//границы для координат 

	bound.push_back(mn);
	bound.push_back(mx);
	
	// границы для радиуса
bound.push_back(min(xmx-xmn, ymx-ymn)/2);
	return bound;
}
//генератор окружностей
vector<double> gencir(vector<double> tr)
{
	vector<double> b=gen_bounderies(tr);
	coor.push_back(gen_chislo(b[0],b[1]));
	coor.push_back(gen_chislo(b[0],b[1]));
	
	return coor;
}

//генератор особи

person gen_person(vector<double> tr)
{
	
	for(int i=0; i<n; i++)
	{
		vector<double>c=gencir(tr);
		
		for(int j=0; j<2; j++)
		{
			coor[i]=c[j];
		}
		
		
	}
	for(int i=0; i<n-1; i++)
	arr[i]=sqrt((coor[i]-coor[i+2])*(coor[i]-coor[i+2])+(coor[i+1]-coor[i+3])*(coor[i+1]-coor[i+3]));
	 arr[n-1]=sqrt((coor[0]-coor[2*n-2])*(coor[0]-coor[2*n-2])+(coor[1]-coor[2*n-1])*(coor[1]-coor[2*n-1]));
	vector<double> s= quickSortR(0,n-1);
	s[2]=s[2]/2;
	for(int i=1; i<n; i++)
		s[3*i+2]-=s[3*i-1];
	
	return s;
	
}
//знаковая площадь
double triansign (double x1, double y1, double x2, double y2, double x3, double y3)
{
		return 0.5*(x2 - x1) * (y3 - y1) - 0.5*(y2 - y1) * (x3 - x1);
}
//площадь сектора образованного стороной треугольника и окружностью которая ее пересекает
double plosh (double d, double r)
{
	return r*r*acos(d/r)-r*d*sqrt(1-d*d/(r*r));
}
//находим расстояние от центра окружности до стороны
double searchd(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double a=y1-y2, b=x2-x1, c=x1*y2-x2*y1;
	return abs(a*x3+b*y3+c)/sqrt(a*a+b*b);
}
// находим суммарную площадь выхода за границы треугольника
double ss1(vector<double> tr, person p)
{
	double s1=0;
	
		for (int i=0; i<p.size(); i+=3)
	{
	
		//находим расстояние от центра данной окружности до сторон треугольника
		vector<double>t(3);
		 t[0]=searchd(tr[0],tr[1],tr[2],tr[3],p[i],p[i+1]);
		 t[1]=searchd(tr[2],tr[3],tr[4],tr[5],p[i],p[i+1]);
		 t[2]=searchd(tr[4],tr[5],tr[0],tr[1],p[i],p[i+1]);
		//если расстояние меньше радиуса, то сторона пересекает окружностью
		// так же выясним какую тройку образуют две вершины треугольника и центр окружности 
		// суммируем площади сегментов
		
		//если знаковая площадь меньше нуля, значит центр лежит внутри треугольника и к площади суммируем сам сегмент
		//в противном случае, центр лежит вне треугольника, а значит берем разность полной площади и найденного сегмента
	
		if(triansign(tr[0],tr[1],tr[2],tr[3],p[i],p[i+1])<0 && t[0]<=p[i+2])
		{
			s1+=plosh(t[0],p[i+2]);
			
		}
			
		if(triansign(tr[0],tr[1],tr[2],tr[3],p[i],p[i+1])>=0 && t[0]<=p[i+2])
		{
			s1+=acos(-1)*p[i+2]*p[i+2]-plosh(t[0],p[i+2]);
			
		}	
		
		if(triansign(tr[2],tr[3],tr[4],tr[5],p[i],p[i+1])<0 && t[1]<=p[i+2])
		{
			s1+=plosh(t[1],p[i+2]);
			
		}
			
		if(triansign(tr[2],tr[3],tr[4],tr[5],p[i],p[i+1])>=0 && t[1]<=p[i+2])
		{
			s1+=acos(-1)*p[i+2]*p[i+2]-plosh(t[1],p[i+2]);
			
		}
			
		if(triansign(tr[4],tr[5],tr[0],tr[1],p[i],p[i+1])<0 && t[2]<=p[i+2])
		{
			s1+=plosh(t[2],p[i+2]);
			
		}
			
		if(triansign(tr[4],tr[5],tr[0],tr[1],p[i],p[i+1])>=0 && t[2]<=p[i+2])
		{
			s1+=acos(-1)*p[i+2]*p[i+2]-plosh(t[2],p[i+2]);
			
		}	
		
	}
	return s1;
}
// находим суммарную площадь пересечений окружностей
double ss2(vector<double> tr, person p)
{
		double s2=0;
		
	for (int i=0; i<p.size()-1; i+=3)
	{
		for (int j=i+3; j<p.size(); j+=3)
		{
			//расстояние между центрами двух окружностей
			double ro=sqrt((p[i]-p[j])*(p[i]-p[j])+(p[i+1]-p[j+1])*(p[i+1]-p[j+1]));
			//если окружности пересекаются
			if(p[i+2]+p[j+2]>ro && ro+p[j+2]>p[i+2] && p[i+2]+ro>p[j+2])
			{			
				double f1=2*acos((p[i+2]*p[i+2]-p[j+2]*p[j+2]+ro*ro)/(2*p[i+2]*ro));
				double f2=2*acos((p[j+2]*p[j+2]-p[i+2]*p[i+2]+ro*ro)/(2*p[j+2]*ro));
				//cout<<f1<<' '<<f2<<' '<<i<<' '<<j<<endl;
				double l=p[i+2]*p[i+2]*(f1-sin(f1))/2+p[j+2]*p[j+2]*(f2-sin(f2))/2;
				s2+=l;
			}
		}
}

return s2;
}
//находим суммарную площадь всех окружностей
double ss3(vector<double> tr, person p)
{
	double s3=0;
		for (int i=2; i<p.size(); i+=3)
	{
		s3+=acos(-1)*p[i]*p[i];
	}
	return s3;
}
//s3-s2-s1
double fin(double alpha, double beta, vector<double>tr, person p)
{
	
	return ss3(tr, p)-alpha*ss1(tr,  p)-beta*ss2(tr, p);
}
int main()
{
	//xa,ya,xb,yb,xc,yc;
	//n
	//vector<double>xo(n),yo(n),ro(n)
	
	//s1 - outer from triangle
	//s2 - intersection
	//s3->max sum of squares
	//f = s3 - alpha*s1 - beta*s2	
	
	// f(xa,ya,xb,yb,xc,yc,xo,yo,ro,alpha,beta)
	
	

	fstream in("in.txt");	
	fstream out("out.txt");
	vector<double>trien(6);

	for(int i=0; i<6; i++)
	in>>trien[i];
	//for(int i=0; i<3*n; i++)
	//	a>>pq[i];


	//cout<<fin(xa,ya,xb,yb,xc,yc, xo, yo, ro);
	cout<<'!';
	double alpha, beta, max_iter=10;
	for(alpha=101;alpha<=1001;alpha+=100)
	{
		cout<<'.';
		for(beta=101;beta<=1001;beta+=100)
		{
			double av_f=0;
			for(int exp_i=0;exp_i<100;++exp_i)
			{
				vector<person>p;
				for(int i=0;i<10;++i)
				{
					person x=gen_person(trien);
					p.push_back(x);
				}
				//сгенерировали 10 особей (каждая это есть система из трех окружностей)
				
				//диф эволюция
				//параметры
				double w=0.4;//весовой коэф-т
				double cr=0.9;//параметр скрещивания
				
				
				for(int iter=0;iter<max_iter;++iter)
				{
					/**
					ШАГИ МЕТОДА diffevo
					оценочная функция f() вызывается с параметрами alpha и beta
					**/
					vector<person> x(3, person(n*3));
					for(int j=0; j<10; ++j)
					{
						for(int k=0; k<3; ++k)
						{
							int help_num = rand()%10;
							while(help_num==j)
							{
								help_num = rand()%10;
							}
							for(int l=0; l<n*3; ++l)
							x[k][l]=p[help_num][l];
						}
						person x_cl(n*3), x_s(n*3);
						for(int k=0; k<n*3; ++k)
						{
							x_cl[k]=x[2][k]+w*(x[0][k]-x[1][k]);
							//проверить на попадание в квадрат	
							vector<double>b=gen_bounderies(trien);
							if((k+1)%3==0)
							{
								if(x_cl[k]>b[2] || x_cl[k]<0)
									x_cl[k]=gen_chislo(0,b[2]);
								
							}	
							else
							{
								if(x_cl[k]>b[1] || x_cl[k]<b[0])
									x_cl[k]=gen_chislo(b[0],b[1]);
							}
								
							
							
						}
						for(int i=0; i<n*3-1; ++i)
						{
							double u_i=gen_chislo(0,1);
							if(u_i<=cr)
								x_s[i]=x_cl[i];
							else
								x_s[i]=p[j][i];
								
							
						}
						
						x_s[n*3-1]=x_cl[n*3-1];
						if(fin(alpha, beta, trien, x_s)>fin(alpha, beta, trien, p[j]))
							p[j]=x_s;
						
					}
					//w=w+0.6/max_iter;
					//cr=cr-0.8/max_iter;
				}
				person best_x=p[0];
				double best_f=fin(alpha,beta,trien, p[0]);
				for(int i=1;i<10;++i){
					if(fin(alpha,beta,trien, p[i])>best_f){
						best_f=fin(alpha,beta,trien, p[i]);
						best_x=p[i];
						
					}
				}
				
				av_f+=best_f;
				
			}
			av_f=av_f*1.0/100;
			out<<alpha<<";"<<beta<<";"<<av_f<<";"<<endl;
		}

	}

}