#include "cuttingGlass.h"
#include "../solution_generator/make_solutionfile.h"
#include <vector>
#include <iostream>
#include <string>
extern Data* DATA;

namespace evalBatch2
{
  extern int plateItem[MAX_ITEMS],i0Item[MAX_ITEMS],xminItem[MAX_ITEMS];
}


namespace gensol2
{

typedef struct {
int plate;
int item;
int x0;
int y0;
int x1;
int y1;
int xnextcut1;
int ynextcut2;
int nextcut;
int rectangle;
} SOL ;

SOL config[MAX_ITEMS];
SOL bestCut1[MAX_ITEMS];
SOL bestPlate[MAX_ITEMS];

char nodefectX[PLATE_WIDTH + 1];

int nbdp;
int Dx0[15];
int Dy0[15];
int Dx1[15];
int Dy1[15];
void initDefectsPlate(int pId)
{
	int d,x;
	for (x=0;x<=PLATE_WIDTH;x++)  nodefectX[x]=1;
	nbdp=0;
	for (d=0;d<DATA->nb_defects_;d++)
		if (DATA->defect_list_[d].plateId==pId)
			{
			Dx0[nbdp]=DATA->defect_list_[d].x;
			Dx1[nbdp]=DATA->defect_list_[d].x+DATA->defect_list_[d].width;
			Dy0[nbdp]=DATA->defect_list_[d].y;
			Dy1[nbdp]=DATA->defect_list_[d].y+DATA->defect_list_[d].height;
			for (x=Dx0[nbdp]+1;x<Dx1[nbdp];x++) nodefectX[x]=0;
			nbdp++;
			}
}

int nextY_Up(int x,int y,int w,int h)
	{
	int d,y0;
	y0=y;
recom:
	for (d=0;d<nbdp;d++)
		{
		if (x>=Dx1[d]) continue;
		if (y>=Dy1[d]) continue;
		if (x+w<=Dx0[d]) continue;
		if (y+h<=Dy0[d]) continue;
		if (Dx1[d]>x+w) return 9999;
		y=Dy1[d];
		if (y-y0<TWENTY) y=y0+TWENTY;
		if (y+h<PLATE_HEIGHT) goto recom;
		}
	return y;
	}

int nextX_Right(int x,int y,int w,int h)
	{
	int d,x0;
	x0=x;
recom:
	for (d=0;d<nbdp;d++)
		{
		if (x>=Dx1[d]) continue;
		if (y>=Dy1[d]) continue;
		if (x+w<=Dx0[d]) continue;
		if (y+h<=Dy0[d]) continue;
		x=Dx1[d];
		if (x-x0<TWENTY) x=x0+TWENTY;
		if (x+w<PLATE_WIDTH) goto recom;
		}
	return x;
	}

int yDefectY(const int y,const int x0,const int x1)
	{
	int d;
	for (d=0;d<nbdp;d++)
		if (Dy1[d]>y && Dy0[d]<y)
			if (Dx1[d]>x0 && Dx0[d]<x1) 
				return Dy1[d];
	return 0;
	}

int noDefectInsideSquare(int x,int y,int x1,int y1)
	{
	int d;
	for (d=0;d<nbdp;d++)
		if (x<Dx1[d] && x1>Dx0[d]  && y<Dy1[d] && y1>Dy0[d]) return 0;
	return 1;
	}

int noDefectOnCut2(const int ycut2,const int x0,const int x1)
	{
	int d;
	for (d=0;d<nbdp;d++)
		if (Dy1[d]>ycut2 && Dy0[d]<ycut2 && Dx1[d]>x0 && Dx0[d]<x1) 
			return 0;			
	return 1;
	}


int noDefectOnVerticalCut(const int x,const int y0,const int y1)
	{
	int d;
	if (nodefectX[x]) return 1;
	for (d=0;d<nbdp;d++)
		if (Dx1[d]>x && Dx0[d]<x && Dy1[d]>y0 && Dy0[d]<y1) 
			return 0;			
	return 1;
	}

int yDefect(const int y,const int x0,const int x1)
{
	int d,ymin=PLATE_HEIGHT;
	for (d=0;d<nbdp;d++)
		if (Dx1[d]>x0 && Dx0[d]<x1 && Dy0[d]>y && Dy0[d]<ymin) ymin=Dy0[d];
	return ymin;
}

int yDefect(const int y,const int x)
{
	int d,ymin=PLATE_HEIGHT;
	for (d=0;d<nbdp;d++)
		if (Dx1[d]>x && Dx0[d]<x && Dy0[d]>y && Dy0[d]<ymin) ymin=Dy0[d];
	return ymin;
}
int xDefect(const int x,const int y)
{
	int d,xmin=PLATE_WIDTH;
	for (d=0;d<nbdp;d++)
		if (Dy1[d]>y && Dy0[d]<y && Dx0[d]>x && Dx0[d]<xmin) xmin=Dx0[d];
	return xmin;
}

#define MAX_CUT_SIZE 3500
int calcMaxXForNextCut1(int X0CUT1) 
	{
	int xmaxcut1;
	xmaxcut1=X0CUT1+MAX_CUT_SIZE;
	if (xmaxcut1>=PLATE_WIDTH) return PLATE_WIDTH;
	if (xmaxcut1>PLATE_WIDTH_20) 
		xmaxcut1=PLATE_WIDTH_20;	
	while (!nodefectX[xmaxcut1]) xmaxcut1--;
	return xmaxcut1;
	}

int chkHz(Solution* sol,int ymax_, int ymax,const int xnc1)
	{
	int y;
	if (ymax>PLATE_HEIGHT_20 && ymax<PLATE_HEIGHT) return 9999;
	if (ymax - ymax_<TWENTY)
		{
		ymax+=TWENTY;
		if (ymax>PLATE_HEIGHT) return ymax;
		while (ymax<PLATE_HEIGHT) 
			{
			y=yDefectY(ymax,sol->x0cut1,xnc1);
			if (y==0) break;
			ymax=y;
			}
		if (ymax>PLATE_HEIGHT_20) return PLATE_HEIGHT;
		return ymax;
		}
	if (ymax<PLATE_HEIGHT)
		{
		y=yDefectY(ymax,sol->x0cut1,xnc1);
		if (y==0) return ymax;
		ymax+=TWENTY;
		if (ymax>PLATE_HEIGHT) return ymax;
		
		while (ymax<PLATE_HEIGHT) 
			{
			y=yDefectY(ymax,sol->x0cut1,xnc1);
			if (y==0) break;
			ymax=y;
			}
		if (ymax>PLATE_HEIGHT_20) return PLATE_HEIGHT;
		}
	return ymax;
	}

int chkVt(Solution* sol,const int xmax_, int xmax)
	{
	if (xmax>PLATE_WIDTH_20 && xmax<PLATE_WIDTH) return 9999;
	if (xmax - xmax_<TWENTY)
		{
		xmax+=TWENTY;
		if (xmax>PLATE_WIDTH) return xmax;
		while (xmax<PLATE_WIDTH && !nodefectX[xmax]) xmax++;
		if (xmax>PLATE_WIDTH_20) return PLATE_WIDTH;
		return xmax;
		}
	if (xmax<PLATE_WIDTH && !nodefectX[xmax])
		{
		xmax+=TWENTY;
		if (xmax>PLATE_WIDTH) return xmax;
		while (xmax<PLATE_WIDTH && !nodefectX[xmax]) xmax++;
		if (xmax>PLATE_WIDTH_20) return PLATE_WIDTH;
		}
	return xmax;
	}

int evalCut3(Solution* sol,int &surf,int &nextx,int &xmax, int &y0,int &ymax)
{
	if(((xmax - nextx)*y0)+((nextx - sol->x0cut1) * ymax) - surf <= sol->wastemin) return 1;
	return 0;
}

int evalCut2(Solution* sol,int &surf,int &xmax,int &ymax)
{
	if (((xmax - sol->x0cut1) * ymax) - surf <= sol->wastemin) return 1;
	return 0;
}


void evalCut1(Solution* sol,int i,int xmax,int &surf)
{
	int waste = ((sol->xmaxcut1 - sol->x0cut1) * PLATE_HEIGHT) - surf;
	if (waste < sol->wastemin || (waste==sol->wastemin && surf>sol->surfCut1))
		{
		sol->wastemin = waste;
		sol->surfCut1 = surf;
		for (int i=sol->from_i0;i<=sol->to_i1;i++)
			bestCut1[i]=config[i];
	} 
}

#define EVALUATE_SHIFT_UP {\
	yd=yDefect(nexty,x0,nextx);\
	if (yd>=ymax && noDefectOnVerticalCut(nextx,y0,y))\
		{\
		YUB_=YUB;\
		if (yd<YUB) YUB=yd;\
		xmaxsave=xmax;xmax_save=xmax_;\
		if (nextx>xmax) {xmax_=xmax;xmax=nextx;}else if (nextx!=xmax && nextx>xmax_) xmax_=nextx;\
		xnc1=chkVt(sol, xmax_,xmax);\
		if (xnc1<=XUB)\
			{\
			ymaxsave=ymax;ymax_save=ymax_;\
			if (nexty>ymax) {ymax_=ymax;ymax=nexty;}else if (nexty!=ymax && nexty>ymax_) ymax_=nexty;\
			ync2=chkHz(sol, ymax_,ymax,xnc1);\
			if (ync2<=YUB)\
				{\
				storeItemConf(x0,y)\
				config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;\
				config[i].ynextcut2=ync2;\
				xd=xDefect(nextx,ync2);\
				if (xd>=xnc1)\
					{\
					if (evalCut2(sol, surf,xnc1,ync2)!=0)\
						{\
						if (xd<XUB) enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,xd,ymax_,ymax,PLATE_HEIGHT);\
						else enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,PLATE_HEIGHT);\
						}\
					}\
				if (evalCut3(sol, surf,nextx,xnc1,y0,ync2)!=0)\
					enumCut1(sol,cut2, i+1,nextx,y0,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				}\
			ymax=ymaxsave;ymax_=ymax_save;\
			}\
		YUB=YUB_;\
		xmax=xmaxsave;xmax_=xmax_save;\
		}\
	}

#define EVALUATE_SHIFT_RIGHT {\
	yd=yDefect(nexty,nextx);\
	if (yd>=ymax)\
		{\
		YUB_=YUB;\
		if (yd<YUB) YUB=yd;\
		xmaxsave=xmax;xmax_save=xmax_;\
		if (nextx>xmax) {xmax_=xmax;xmax=nextx;}else if (nextx!=xmax && nextx>xmax_) xmax_=nextx;\
		xnc1=chkVt(sol, xmax_,xmax);\
		if (xnc1<=XUB)\
			{\
			ymaxsave=ymax;ymax_save=ymax_;\
			if (nexty>ymax) {ymax_=ymax;ymax=nexty;}else if (nexty!=ymax && nexty>ymax_) ymax_=nexty;\
			ync2=chkHz(sol, ymax_,ymax,xnc1);\
			if (ync2<=YUB)\
				{\
				storeItemConf(x,y0)\
				config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;\
				config[i].ynextcut2=ync2;\
				tryCut4(sol,cut2, i+1,x,y0,nextx,nexty,xnc1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				xd=xDefect(nextx,ync2);\
				if (xd>=xnc1)\
					{\
					if (evalCut2(sol, surf,xnc1,ync2)!=0)\
						{\
						if (xd<XUB) enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,xd,ymax_,ymax,PLATE_HEIGHT);\
						else enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,PLATE_HEIGHT);\
						}\
					}\
				if (evalCut3(sol, surf,nextx,xnc1,y0,ync2)!=0)\
					enumCut1(sol,cut2,i+1,nextx,y0,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				}\
			ymax=ymaxsave;ymax_=ymax_save;\
			}\
		YUB=YUB_;\
		xmax=xmaxsave;xmax_=xmax_save;\
		}\
	}

#define EVALUATE {\
	yd=yDefect(nexty,nextx);\
	if (yd>=ymax)\
		{\
		YUB_=YUB;\
		if (yd<YUB) YUB=yd;\
		xmaxsave=xmax;xmax_save=xmax_;\
		if (nextx>xmax) {xmax_=xmax;xmax=nextx;}else if (nextx!=xmax && nextx>xmax_) xmax_=nextx;\
		xnc1=chkVt(sol, xmax_,xmax);\
		if (xnc1<=XUB)\
			{\
			ymaxsave=ymax;ymax_save=ymax_;\
			if (nexty>ymax) {ymax_=ymax;ymax=nexty;}else if (nexty!=ymax && nexty>ymax_) ymax_=nexty;\
			ync2=chkHz(sol, ymax_,ymax,xnc1);\
			if (ync2<=YUB)\
				{\
				storeItemConf(x0,y0)\
				config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;\
				config[i].ynextcut2=ync2;\
				tryCut4(sol,cut2,i+1,x0,y0,nextx,nexty,xnc1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				xd=xDefect(nextx,ync2);\
				if (xd>=xnc1)\
					{\
					if (evalCut2(sol, surf,xnc1,ync2)!=0)\
						{\
						if (xd<XUB)enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,xd,ymax_,ymax,PLATE_HEIGHT);\
						else enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,PLATE_HEIGHT);\
						}\
					}\
				if (evalCut3(sol, surf,nextx,xnc1,y0,ync2)!=0)\
					enumCut1(sol,cut2,i+1,nextx,y0,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				}\
			ymax=ymaxsave;ymax_=ymax_save;\
			}\
		YUB=YUB_;\
		xmax=xmaxsave;xmax_=xmax_save;\
		}\
		}




#define storeItemConf(x,y) \
	config[i].item=j;\
	config[i].x0=x;config[i].y0=y;\
	config[i].x1=nextx;config[i].y1=nexty;



#define nextxCORRECT (nextx<=XUB && (nextx<=PLATE_WIDTH_20 || nextx==PLATE_WIDTH))
#define nextyCORRECT (nexty<=YUB && (nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT))

void enumCut1(Solution* sol,int cut2,const int i,int x0,int y0,int surf,int xmax_,int xmax,int XUB,int ymax_,int ymax,int YUB);
void tryCut4(Solution* sol,int cut2, const int i,int x0,int y0,int nextx,int y1,int &xnc1,int &ync2,int surf,int& xmax_,int& xmax,int XUB,int& ymax_,int& ymax,int YUB)
	{
	int j,dx,dy,nexty,ok,xd;
	if (i<DATA->nb_items_)
		{
		j = sol->seq_[i];
		dy=DATA->item_list_[j].width;
		dx=DATA->item_list_[j].length;
		if (dx==nextx-x0)
			{
			nexty=y1+dy;
			if nextyCORRECT
				{
				ok=0;
				if (sol->nbcut4[cut2]!=0)
					{
					if (nexty==ymax && noDefectInsideSquare(x0,y1,nextx,nexty)) ok=1;
					}
				else
					{
					if (nexty>=ymax+TWENTY)
						{
						if (noDefectInsideSquare(x0,y1,nextx,nexty) && noDefectOnCut2(nexty,sol->x0cut1,xnc1)) ok=1;
						}
					else
					if (nexty==ymax && ymax>=ymax_+TWENTY && noDefectInsideSquare(x0,y1,nextx,nexty)) ok=1;
					}	
				if (ok)
					{
					xd=xDefect(nextx,nexty);
					if (xd>=xnc1)
						{
						sol->nbcut4[cut2]++;
						surf=surf+DATA->item_list_[j].surf;
						storeItemConf(x0,y1)
						config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;
						config[i].ynextcut2=nexty;
						if (xd<XUB) XUB=xd;
						if (evalCut2(sol, surf,xnc1,nexty)!=0)
							enumCut1(sol,cut2+1,i+1,sol->x0cut1,nexty,surf,xmax_,xmax,XUB,ymax_,nexty,PLATE_HEIGHT);
						if (evalCut3(sol,surf,nextx,xnc1,y0,nexty)!=0)
							enumCut1(sol,cut2,i+1,        nextx,   y0,surf,xmax_,xmax,XUB,ymax_,nexty,nexty);
						sol->nbcut4[cut2]--;
						}
					}
				}
			return;
			}
		if (dy==nextx-x0)
			{
			nexty=y1+dx;
			if nextyCORRECT
				{
				ok=0;
				if (sol->nbcut4[cut2]!=0)
					{
					if (nexty==ymax && noDefectInsideSquare(x0,y1,nextx,nexty)) ok=1;
					}
				else
					{
					if (nexty>=ymax+20)
						{
						if (noDefectInsideSquare(x0,y1,nextx,nexty) && noDefectOnCut2(nexty,sol->x0cut1,xnc1)) ok=1;
						}
					else
					if (nexty==ymax && ymax>=ymax_+TWENTY && noDefectInsideSquare(x0,y1,nextx,nexty)) ok=1;
					}
				if (ok)
					{
					xd=xDefect(nextx,nexty);
					if (xd>=xnc1)
						{
						sol->nbcut4[cut2]++;
						surf=surf+DATA->item_list_[j].surf;
						storeItemConf(x0,y1)
						config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;
						config[i].ynextcut2=nexty;
						if (xd<XUB) XUB=xd;
						if (evalCut2(sol, surf,xnc1,nexty)!=0)
							enumCut1(sol,cut2+1,i+1,sol->x0cut1,nexty,surf,xmax_,xmax,XUB,ymax_,nexty,PLATE_HEIGHT);
						if (evalCut3(sol,surf,nextx,xnc1,y0,nexty)!=0)
							enumCut1(sol,cut2,i+1,        nextx,   y0,surf,xmax_,xmax,XUB,ymax_,nexty,nexty);
						sol->nbcut4[cut2]--;
						}
					}
				}
			}
		}
}

void enumCut1(Solution* sol,int cut2,const int i,int x0,int y0,int surf,int xmax_,int xmax,int XUB,int ymax_,int ymax,int YUB)
	{
	int j,dx,dy,x,nextx,y,nexty,xd,yd,YUB_;
	int xmaxsave,xmax_save,ymaxsave,ymax_save,xnc1,ync2;
	if (i<=sol->to_i1)
		{
		j = sol->seq_[i];
		dy=DATA->item_list_[j].width;
		dx=DATA->item_list_[j].length;
		surf=surf+DATA->item_list_[j].surf;
		
		nextx=x0+dx;
		if nextxCORRECT
			{
			y=nextY_Up(x0,y0,dx,dy);
			nexty=y+dy;			
			if (y!=y0)
				{
				if nextyCORRECT
					EVALUATE_SHIFT_UP
				nexty=y0+dy;
				if nextyCORRECT
					{
					x=nextX_Right(x0,y0,dx,dy);
					nextx=x+dx;
					if nextxCORRECT
						EVALUATE_SHIFT_RIGHT
					}
				}
			else
			if nextyCORRECT
				EVALUATE
			}
		if (dx!=dy)
			{
			nextx=x0+dy;
			if nextxCORRECT
				{
				y=nextY_Up(x0,y0,dy,dx);
				nexty=y+dx;
				if (y!=y0)
					{
					if nextyCORRECT
						EVALUATE_SHIFT_UP
					nexty=y0+dx;
					if nextyCORRECT
						{
						x=nextX_Right(x0,y0,dy,dx);
						nextx=x+dy;
						if nextxCORRECT
							EVALUATE_SHIFT_RIGHT
						}
					}
				else
				if nextyCORRECT
					EVALUATE
				}
			}
		}
	else
		{
		evalCut1(sol,i,xmax,surf);
		}
	}


std::string batchPath;
std::string defectsPath;
std::string optParamsPath;

int to_i1[7000];

int evalBatch(Solution* sol)
	{
	int cut1,np,i,j,nbcut1,ymax,rectangle_;

	i=DATA->nb_items_-1;
	nbcut1=0;
nextcut1:
	j=evalBatch2::i0Item[i];
	to_i1[nbcut1++]=i;
	i=j-1;
	if (i>=0) goto nextcut1;

	cut1=nbcut1-1;
	sol->from_i0=0;
	np=0;
	sol->numplate=evalBatch2::plateItem[sol->from_i0];
	initDefectsPlate(np);
	sol->x0cut1=0;
	sol->numcut1=0;
	sol->to_i1=to_i1[cut1];
	sol->xmaxcut1=evalBatch2::xminItem[sol->to_i1];

enumcut1:

	sol->wastemin=INT_MAX;
	enumCut1(sol,0,sol->from_i0, sol->x0cut1, 0, 0,
	0, 0, sol->xmaxcut1, 
	0, 0, PLATE_HEIGHT);
	if (sol->wastemin<INT_MAX)
		{
		for (int i=sol->from_i0;i<=sol->to_i1;i++)
			{
			bestPlate[i]=bestCut1[i];
			bestPlate[i].plate=sol->numplate;
			bestPlate[i].xnextcut1=sol->xmaxcut1;
			}
		cut1--;
		if (cut1>=0)
			{
			sol->numcut1++;
			sol->from_i0=sol->to_i1+1;
			sol->to_i1=to_i1[cut1];

			if (evalBatch2::plateItem[sol->from_i0]!=np)
				{
				np=sol->numplate=evalBatch2::plateItem[sol->from_i0];
				initDefectsPlate(np);
				sol->numcut1=0;
				sol->x0cut1=0;
				}
			else
				{
				sol->x0cut1=sol->xmaxcut1;
				}
			sol->xmaxcut1=evalBatch2::xminItem[sol->to_i1];

			goto enumcut1;
			}
		}
	else 
		return -2;

		

					rectangle_=-1;
					for (i=0;i<DATA->nb_items_;i++)
						{

							rectangle_=bestPlate[i].rectangle;
							ymax=bestPlate[i].ynextcut2;
							for (j=0;j<DATA->nb_items_;j++)
								{
								if(bestPlate[j].rectangle==rectangle_ && ymax<bestPlate[j].ynextcut2)
									ymax=bestPlate[j].ynextcut2;
								}

	
							bestPlate[i].ynextcut2=ymax;


						}

					std::vector<std::vector<int> > data(DATA->nb_items_);
					for (i=0;i<DATA->nb_items_;i++)
					{
						data[i].push_back(bestPlate[i].plate);
						data[i].push_back(bestPlate[i].item);
						data[i].push_back(bestPlate[i].xnextcut1);
						data[i].push_back(bestPlate[i].ynextcut2);
						data[i].push_back(bestPlate[i].x0);
						data[i].push_back(bestPlate[i].y0);
						data[i].push_back(bestPlate[i].x1);
						data[i].push_back(bestPlate[i].y1);
					}

					for (i=0;i<DATA->nb_items_;i++)
						if (bestPlate[i].y1>bestPlate[i].ynextcut2) return -2;

				if (make_solution_file(data,sol->solfile))
				{
					return 1;
				}
				return -2;
	}

}
