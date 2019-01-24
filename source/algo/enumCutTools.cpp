#define PLATE_WIDTH 6000
#define PLATE_WIDTH_20 5980
#define PLATE_HEIGHT 3210
#define PLATE_HEIGHT_20 3190
#define PLATE_SURFACE 19260000
#define TWENTY 20
#define MIN_CUT_SIZE 100
#define MAX_PLATES 100
#define MAX_ITEMS 700
#define MAX_DEFECTS 1000
#define NMB_PLATES 100
#define MAXCUTSIZE 3500
#define MAX_CUT_SIZE 3500

#include "cuttingGlass.h"

extern Data* DATA;

char nodefectX[PLATE_WIDTH + 1];

int NBDP;
int DX0[15];
int DY0[15];
int DX1[15];
int DY1[15];
void initDefectsPlate(int pId)
{
	int d,x;
	for (x=0;x<=PLATE_WIDTH;x++)  nodefectX[x] = 1;
	NBDP=0;
	for (d=0;d<DATA->nb_defects_;d++)
		if (DATA->defect_list_[d].plateId==pId)
			{
			DX0[NBDP]=DATA->defect_list_[d].x;
			DX1[NBDP]=DATA->defect_list_[d].x+DATA->defect_list_[d].width;
			DY0[NBDP]=DATA->defect_list_[d].y;
			DY1[NBDP]=DATA->defect_list_[d].y+DATA->defect_list_[d].height;
			for (x=DX0[NBDP] + 1;x<DX1[NBDP];x++) nodefectX[x]=0;
			NBDP++;
			}
}


void initDefectsPlate(Solution* sol, int x0, int x1)
{
	sol->nbdp=0;
	for (int d=0;d<NBDP;d++)
	{
		if(DX0[d] <= x1 && DX1[d] >= x0)
		{
			sol->Dx0[sol->nbdp]=DX0[d];
			sol->Dx1[sol->nbdp]=DX1[d];
			sol->Dy0[sol->nbdp]=DY0[d];
			sol->Dy1[sol->nbdp]=DY1[d];
			sol->nbdp++;
		}
	}
}

void initDefectsPlateGenSol(Solution* sol, int pId)
{
	int d,x;
	for (x=0;x<=PLATE_WIDTH;x++)  nodefectX[x] = 1;
	sol->nbdp=0;
	for (d=0;d<DATA->nb_defects_;d++)
		if (DATA->defect_list_[d].plateId==pId)
			{
			sol->Dx0[sol->nbdp]=DATA->defect_list_[d].x;
			sol->Dx1[sol->nbdp]=DATA->defect_list_[d].x+DATA->defect_list_[d].width;
			sol->Dy0[sol->nbdp]=DATA->defect_list_[d].y;
			sol->Dy1[sol->nbdp]=DATA->defect_list_[d].y+DATA->defect_list_[d].height;
			for (x=sol->Dx0[sol->nbdp] + 1;x<sol->Dx1[sol->nbdp];x++) nodefectX[x]=0;
			sol->nbdp++;
			}
}

int nextY_Up(Solution* sol, int x,int y,int w,int h)
	{
	int d,y0;
	y0=y;
recom:
	for (d=0;d<sol->nbdp;d++)
		{
		if (x>=sol->Dx1[d]) continue;
		if (y>=sol->Dy1[d]) continue;
		if (x+w<=sol->Dx0[d]) continue;
		if (y+h<=sol->Dy0[d]) continue;

		if (sol->Dx1[d]>=x + w) return 9999;

		y=sol->Dy1[d];
		if (y-y0<TWENTY) y=y0+TWENTY;
		if (y+h<=PLATE_HEIGHT) goto recom;
		}
	return y;
	}

int nextX_Right(Solution* sol, int x,int y,int w,int h)
	{
	int d,x0;
	x0=x;
recom:
	for (d=0;d<sol->nbdp;d++)
		{
		if (x>=sol->Dx1[d]) continue;
		if (y>=sol->Dy1[d]) continue;
		if (x+w<=sol->Dx0[d]) continue;
		if (y+h<=sol->Dy0[d]) continue;
		x=sol->Dx1[d];
		if (x-x0<TWENTY) x=x0+TWENTY;
		if (x+w<=PLATE_WIDTH) goto recom;
		}
	return x;
	}

int shiftCut2OverDefectY(Solution* sol, const int y,const int x0,const int x1)
	{
	int d;
	for (d=0;d<sol->nbdp;d++)
		if (sol->Dy1[d]>y && sol->Dy0[d]<y)
			if (sol->Dx1[d]>x0 && sol->Dx0[d]<x1)
				return sol->Dy1[d];
	return 0;
	}

int noDefectInsideSquare(Solution* sol, int x,int y,int x1,int y1)
	{
	int d;
	for (d=0;d<sol->nbdp;d++)
		if (x<sol->Dx1[d] && x1>sol->Dx0[d]  && y<sol->Dy1[d] && y1>sol->Dy0[d]) return 0;
	return 1;
	}

int noDefectOnCut2(Solution* sol, const int ycut2,const int x0,const int x1)
	{
	int d;
	for (d=0;d<sol->nbdp;d++)
		if (sol->Dy1[d]>ycut2 && sol->Dy0[d]<ycut2 && sol->Dx1[d]>x0 && sol->Dx0[d]<x1)
			return 0;
	return 1;
	}


int noDefectOnVerticalCut(Solution* sol, const int x,const int y0,const int y1)
	{
	int d;
	for (d=0;d<sol->nbdp;d++)
		if (sol->Dx1[d]>x && sol->Dx0[d]<x && sol->Dy1[d]>y0 && sol->Dy0[d]<y1)
			return 0;			
	return 1;
	}

int preventTooManyCut4(Solution* sol, const int y,const int x0,const int x1)
{
	int d,ymin=PLATE_HEIGHT;
	for (d=0;d<sol->nbdp;d++)
		if (sol->Dx1[d]>x0 && sol->Dx0[d]<x1 && sol->Dy0[d]>y && sol->Dy0[d]<ymin) ymin=sol->Dy0[d];
	return ymin;
}

int preventCut3InDefect(Solution* sol, const int y,const int x)
{
	int d,ymin=PLATE_HEIGHT;
	for (d=0;d<sol->nbdp;d++)
		if (sol->Dx1[d]>x && sol->Dx0[d]<x && sol->Dy0[d]>=y && sol->Dy0[d]<ymin) ymin=sol->Dy0[d];
	return ymin;
}
int xDefect(Solution* sol, const int x,const int y)
{
	int d,xmin=PLATE_WIDTH;
	for (d=0;d<sol->nbdp;d++)
		if (sol->Dy1[d]>y && sol->Dy0[d]<y && sol->Dx0[d]>=x && sol->Dx0[d]<xmin) xmin=sol->Dx0[d];
	return xmin;
}



#define MAX_CUT_SIZE 3500
int calcMaxXForNextCut1(Solution* sol, int X0CUT1)
	{
	int xmaxcut1;
	xmaxcut1=X0CUT1+MAX_CUT_SIZE;
	if (xmaxcut1>=PLATE_WIDTH)
		xmaxcut1=PLATE_WIDTH;
	else
		{
		if (xmaxcut1>PLATE_WIDTH_20)
			xmaxcut1=PLATE_WIDTH_20;
		while(!nodefectX[xmaxcut1]) xmaxcut1--;
		}
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
			y=shiftCut2OverDefectY(sol, ymax,sol->x0cut1,xnc1);
			if (y==0) break;
			ymax=y;
			}
		if (ymax>PLATE_HEIGHT_20) return PLATE_HEIGHT;
		return ymax;
		}
	if (ymax<PLATE_HEIGHT)
		{
		y=shiftCut2OverDefectY(sol, ymax,sol->x0cut1,xnc1);
		if (y==0) return ymax;
		ymax+=TWENTY;
		if (ymax>PLATE_HEIGHT) return ymax;

		while (ymax<PLATE_HEIGHT)
			{
			y=shiftCut2OverDefectY(sol, ymax,sol->x0cut1,xnc1);
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

int chkHz0(Solution* sol,int ymax_, int ymax)
	{
	if (ymax>PLATE_HEIGHT_20 && ymax<PLATE_HEIGHT) return 9999;
	if (ymax - ymax_<TWENTY)
		{
		ymax+=TWENTY;
		if (ymax>PLATE_HEIGHT) return ymax;
		if (ymax>PLATE_HEIGHT_20) return PLATE_HEIGHT;
		}
	return ymax;
	}

int chkVt0(Solution* sol,const int xmax_, int xmax)
	{
	if (xmax>PLATE_WIDTH_20 && xmax<PLATE_WIDTH) return 9999;
	if (xmax - xmax_<TWENTY)
		{
		xmax+=TWENTY;
		if (xmax>PLATE_WIDTH) return xmax;
		if (xmax>PLATE_WIDTH_20) return PLATE_WIDTH;
		}
	return xmax;
	}


int evalCut3(Solution* sol,int &surf,int &nextx,int &xmax, int &y0,int &ymax)
{
	if( ((xmax - nextx)*y0) + ((nextx - sol->x0cut1) * ymax) - surf + DATA->min_item_surf_ <= sol->wastemin) return 1;
	return 0;
}

int evalCut2(Solution* sol,int &surf,int &xmax,int &ymax)
{
	if ( ((xmax - sol->x0cut1) * ymax) - surf + DATA->min_item_surf_  <= sol->wastemin) return 1;
	return 0;
}

void chkFirstItemOfPlate(Solution* sol, const int i)
	{
	int j,dx,dy,x,nextx,y,nexty;
	j = sol->seq_[i];
	dy=DATA->item_list_[j].length;
	dx=DATA->item_list_[j].width;
	y=nextY_Up(sol, 0,0,dx,dy);
	if (y==0) return;		
	nexty=y+dy;
	if (nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT) return;
	x=nextX_Right(sol, 0,0,dx,dy);
	if (x==0) return;
	nextx=x+dx;
	if (nextx<=sol->xmaxcut1) return;
	
	if (dx!=dy)
		{
		y=nextY_Up(sol, 0,0,dy,dx);
		if (y==0) return;
		nexty=y+dx;
		if (nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT) return;
		x=nextX_Right(sol, 0,0,dy,dx);
		if (x==0) return;
		nextx=x+dy;
		if (nextx<=sol->xmaxcut1) return;
		}
			
	sol->x0cut1=nextx-MAXCUTSIZE;
	while(sol->x0cut1<x && !nodefectX[sol->x0cut1]) sol->x0cut1++;
	sol->xmaxcut1 = calcMaxXForNextCut1(sol, sol->x0cut1);
	}
