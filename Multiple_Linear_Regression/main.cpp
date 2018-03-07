//
//  main.cpp
//  Multiple_Linear_Regression
//
//  Created by 许玥 on 2017/10/17.
//  Copyright © 2017年 许玥. All rights reserved.
//

#include "stdio.h"
#include "math.h"
#include "stdlib.h"

int chlk(double a[],int n,int m,double d[])
{
    int i,j,k,u,v ;
    if((a[0]+1.0==1.0)||(a[0]<0.0))
    {
        printf("fail");
        return(-2);
    }
    a[0]=sqrt(a[0]);
    for(j=1;j<=n-1;j++)a[j]=a[j]/a[0];
    for(i=1;i<=n-1;i++)
    {
        u=i*n+i ;
        for(j=1;j<=i;j++)
        {
            v=(j-1)*n+i ;
            a[u]=a[u]-a[v]*a[v];
        }
        if((a[u]+1.0==1.0)||(a[u]<0.0))
        {
            printf("fail");
            return(-2);
        }
        a[u]=sqrt(a[u]);
        if(i!=(n-1))
        {
            for(j=i+1;j<=n-1;j++)
            {
                v=i*n+j ;
                for(k=1;k<=i;k++)
                    a[v]=a[v]-a[(k-1)*n+i]*a[(k-1)*n+j];
                    a[v]=a[v]/a[u];
            }
            
        }
        
    }
    for(j=0;j<=m-1;j++)
    {
        d[j]=d[j]/a[0];
        for(i=1;i<=n-1;i++)
        {
            u=i*n+i ;
            v=i*m+j ;
            for(k=1;k<=i;k++)
                d[v]=d[v]-a[(k-1)*n+i]*d[(k-1)*m+j];
            d[v]=d[v]/a[u];
        }
    }
    for(j=0;j<=m-1;j++)
    {
        u=(n-1)*m+j ;
        d[u]=d[u]/a[n*n-1];
        for(k=n-1;k>=1;k--)
        {
            u=(k-1)*m+j ;
            for(i=k;i<=n-1;i++)
            {
                v=(k-1)*n+i ;
                d[u]=d[u]-a[v]*d[i*m+j];
            }
            v=(k-1)*n+k-1 ;
            d[u]=d[u]/a[v];
        }
    }
    return(2);
}
                              
/*形参与函数类型    参数意义
double  x[m][n] 每一列存放m个自变量的观测值
double  y[n]    存放随机变量y的n个观测值
int  m  自变量个数
int  n  观测数据的组数
double  a[m+1]  返回回归系数a0,a1,…,am-1,am
double  dt[4]   dt(0)返回偏差平方和q，dt(1)返回平均标准偏差s，dt(2)返回复相关系数r，dt(3)返回回归平方和u
double  v[m]    返回m个自变量的偏相关系数
void  sqt1( )   过程
*/

void sqt2(double x[],double y[],int m,int n,double a[],double dt[],double v[])
{
    int i,j,k,mm ;
    double q,e,u,p,yy,s,r,pp,*b ;
    b=(double *)malloc((m+1)*(m+1)*sizeof(double));
    mm=m+1 ;
    b[mm*mm-1]=n ;
    for(j=0;j<=m-1;j++)
    {
        p=0.0 ;
        for(i=0;i<=n-1;i++)
            p=p+x[j*n+i];
        b[m*mm+j]=p ;
        b[j*mm+m]=p ;
    }
    for(i=0;i<=m-1;i++)
        for(j=i;j<=m-1;j++)
        {
            p=0.0 ;
            for(k=0;k<=n-1;k++)
                p=p+x[i*n+k]*x[j*n+k];
            b[j*mm+i]=p ;
            b[i*mm+j]=p ;
        }
    a[m]=0.0 ;
    for(i=0;i<=n-1;i++)
        a[m]=a[m]+y[i];
    for(i=0;i<=m-1;i++)
    {
        a[i]=0.0 ;
        for(j=0;j<=n-1;j++)
            a[i]=a[i]+x[i*n+j]*y[j];
    }
    chlk(b,mm,1,a);
    yy=0.0 ;
    for(i=0;i<=n-1;i++)
        yy=yy+y[i]/n ;
    q=0.0 ;
    e=0.0 ;
    u=0.0 ;
    for(i=0;i<=n-1;i++)
    {
        p=a[m];
        for(j=0;j<=m-1;j++)
            p=p+a[j]*x[j*n+i];
        q=q+(y[i]-p)*(y[i]-p);
        e=e+(y[i]-yy)*(y[i]-yy);
        u=u+(yy-p)*(yy-p);
    }
    s=sqrt(q/n);
    r=sqrt(1.0-q/e);
    for(j=0;j<=m-1;j++)
    {
        p=0.0 ;
        for(i=0;i<=n-1;i++)
        {
            pp=a[m];
            for(k=0;k<=m-1;k++)
                if(k!=j)pp=pp+a[k]*x[k*n+i];
            p=p+(y[i]-pp)*(y[i]-pp);
        }
        v[j]=sqrt(1.0-q/p);
    }
    dt[0]=q ;
    dt[1]=s ;
    dt[2]=r ;
    dt[3]=u ;
    free(b);
    return ;
}
                              
int main()
{
    int i ;
    double a[4],v[3],dt[4];
    double x[3][5]={ {1.1,1.0,1.2,1.1,0.9},{2.0,2.0,1.8,1.9,2.1},{3.2,3.2,3.0,2.9,2.9}};
    double y[5]={10.1,10.2,10.0,10.1,10.0};
    sqt2(&x[0][0],y,3,5,a,dt,v);
    for(i=0;i<=3;i++)
    printf("a(%d)=%13.5e\n",i,a[i]);
    printf("y = %13.5e * x1 + %13.5e * x2 + %13.5e * x3 + %13.5e",a[0],a[1],a[2],a[3]);
    printf("\n\n");
    printf("q=%13.5e\ns=%13.5e\nr=%13.5e\n",dt[0],dt[1],dt[2]);
    printf("\n");
    for(i=0;i<=2;i++)
    printf("v(%d)=%13.5e\n",i,v[i]);
    printf("\n");
    printf("u=%13.5e",dt[3]);
    printf("\n");
}

