#include "cuttingGlass.h"
#include "../solution_generator/make_solutionfile.h"
#include <vector>
#include <iostream>


namespace gensol
{

Data* DATA;

int i0prev;

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

#include "enumCutTools.cpp"

void evalCut1(Solution* sol, const int &i1,int &surf,const int &xEndCut1)
{
	int waste = ((xEndCut1 - sol->x0cut1) * PLATE_HEIGHT) - surf;
	if (waste < sol->wastemin || (waste==sol->wastemin && surf>sol->surfCut1))
	{
		sol->wastemin = waste;
		sol->i0nextcut = i1 + 1;
		sol->x0nextcut = xEndCut1;
		sol->surfCut1 = surf;
		for (int i=i0prev;i<=i1;i++)
			bestCut1[i]=config[i];
		}
}

#define EVALUATE_SHIFT_UP {\
	ymaxsave=ymax;ymax_save=ymax_;\
	if (nexty>ymax) {ymax_=ymax;ymax=nexty;}else if (nexty!=ymax && nexty>ymax_) ymax_=nexty;\
	yd=preventTooManyCut4(sol, nexty,x0,nextx);\
	if (yd>=ymax && noDefectOnVerticalCut(sol, nextx,y0,y))\
		{\
		YUB_=YUB;\
		if (yd<YUB) YUB=yd;\
		xmaxsave=xmax;xmax_save=xmax_;\
		if (nextx>xmax) {xmax_=xmax;xmax=nextx;}else if (nextx!=xmax && nextx>xmax_) xmax_=nextx;\
		xnc1=chkVt(sol, xmax_,xmax);\
		if (xnc1<=XUB)\
			{\
			ync2=chkHz(sol, ymax_,ymax,xnc1);\
			if (ync2<=YUB)\
				{\
				storeItemConf(x0,y)\
				config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;\
				config[i].ynextcut2=ync2;\
				xd=xDefect(sol, nextx,ync2);\
				if (xd>=xnc1)\
					{\
					evalCut1(sol, i,surf,xnc1);\
					if (evalCut2(sol, surf,xnc1,ync2)!=0)\
						{\
						if (xd<XUB) enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,xd,ymax_,ymax,PLATE_HEIGHT);\
						else enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,PLATE_HEIGHT);\
						}\
					}\
				if (evalCut3(sol, surf,nextx,xnc1,y0,ync2)!=0)\
					enumCut1(sol,cut2, i+1,nextx,y0,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				}\
			}\
		YUB=YUB_;\
		xmax=xmaxsave;xmax_=xmax_save;\
		}\
		ymax=ymaxsave;ymax_=ymax_save;\
	}

#define EVALUATE_SHIFT_RIGHT {\
	ymaxsave=ymax;ymax_save=ymax_;\
	if (nexty>ymax) {ymax_=ymax;ymax=nexty;}else if (nexty!=ymax && nexty>ymax_) ymax_=nexty;\
	yd=preventCut3InDefect(sol, nexty,nextx);\
	if (yd>=ymax)\
		{\
		YUB_=YUB;\
		if (yd<YUB) YUB=yd;\
		xmaxsave=xmax;xmax_save=xmax_;\
		if (nextx>xmax) {xmax_=xmax;xmax=nextx;}else if (nextx!=xmax && nextx>xmax_) xmax_=nextx;\
		xnc1=chkVt(sol, xmax_,xmax);\
		if (xnc1<=XUB)\
			{\
			ync2=chkHz(sol, ymax_,ymax,xnc1);\
			if (ync2<=YUB)\
				{\
				storeItemConf(x,y0)\
				config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;\
				config[i].ynextcut2=ync2;\
				xd=xDefect(sol, nextx,ync2);\
				if (xd>=xnc1)\
					{\
					evalCut1(sol, i,surf,xnc1);\
					if (evalCut2(sol, surf,xnc1,ync2)!=0)\
						{\
						if (xd<XUB) enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,xd,ymax_,ymax,PLATE_HEIGHT);\
						else enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,PLATE_HEIGHT);\
						}\
					}\
				if (evalCut3(sol, surf,nextx,xnc1,y0,ync2)!=0)\
					enumCut1(sol,cut2,i+1,nextx,y0,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				tryCut4(sol,cut2, i+1,x,y0,nextx,nexty,xnc1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				}\
			}\
		YUB=YUB_;\
		xmax=xmaxsave;xmax_=xmax_save;\
		}\
		ymax=ymaxsave;ymax_=ymax_save;\
	}

#define EVALUATE {\
	ymaxsave=ymax;ymax_save=ymax_;\
	if (nexty>ymax) {ymax_=ymax;ymax=nexty;}else if (nexty!=ymax && nexty>ymax_) ymax_=nexty;\
	yd=preventCut3InDefect(sol, nexty,nextx);\
	if (yd>=ymax)\
		{\
		YUB_=YUB;\
		if (yd<YUB) YUB=yd;\
		xmaxsave=xmax;xmax_save=xmax_;\
		if (nextx>xmax) {xmax_=xmax;xmax=nextx;}else if (nextx!=xmax && nextx>xmax_) xmax_=nextx;\
		xnc1=chkVt(sol, xmax_,xmax);\
		if (xnc1<=XUB)\
			{\
			ync2=chkHz(sol, ymax_,ymax,xnc1);\
			if (ync2<=YUB)\
				{\
				storeItemConf(x0,y0)\
				config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;\
				config[i].ynextcut2=ync2;\
                xd=xDefect(sol, nextx,ync2);\
				if (xd>=xnc1)\
					{\
					evalCut1(sol, i,surf,xnc1);\
					if (evalCut2(sol, surf,xnc1,ync2)!=0)\
						{\
						if (xd<XUB)enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,xd,ymax_,ymax,PLATE_HEIGHT);\
						else enumCut1(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,PLATE_HEIGHT);\
						}\
					}\
				if (evalCut3(sol, surf,nextx,xnc1,y0,ync2)!=0)\
					enumCut1(sol,cut2,i+1,nextx,y0,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				tryCut4(sol,cut2,i+1,x0,y0,nextx,nexty,xnc1,ync2,surf,xmax_,xmax,XUB,ymax_,ymax,YUB);\
				}\
			}\
		YUB=YUB_;\
		xmax=xmaxsave;xmax_=xmax_save;\
		}\
		ymax=ymaxsave;ymax_=ymax_save;\
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
					if (nexty==ymax && noDefectInsideSquare(sol, x0,y1,nextx,nexty)) ok=1;
					}
				else
					{
					if (nexty>=ymax+20)
						{
						if (noDefectInsideSquare(sol, x0,y1,nextx,nexty) && noDefectOnCut2(sol, nexty,sol->x0cut1,xnc1)) ok=1;
						}
					else
					if (nexty==ymax && ymax>=ymax_+TWENTY && noDefectInsideSquare(sol, x0,y1,nextx,nexty)) ok=1;
					}
				if (ok)
					{
					xd=xDefect(sol, nextx,nexty);
					if (xd>=xnc1)
						{
					sol->nbcut4[cut2]++;
					surf=surf+DATA->item_list_[j].surf;
					storeItemConf(x0,y1)
					config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;
					config[i].ynextcut2=nexty;
					if (xd<XUB) XUB=xd;
					evalCut1(sol, i,surf,xnc1);
					if (evalCut2(sol, surf,xnc1,nexty)!=0)
						enumCut1(sol,cut2+1,i+1,sol->x0cut1,nexty,surf,xmax_,xmax,XUB,ymax_,nexty,PLATE_HEIGHT);
					if (evalCut3(sol,surf,nextx,xnc1,y0,nexty)!=0)
						enumCut1(sol,cut2,i+1,        nextx,   y0,surf,xmax_,xmax,XUB,ymax_,nexty,nexty);
					sol->nbcut4[cut2]--;
					}}
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
					if (nexty==ymax && noDefectInsideSquare(sol, x0,y1,nextx,nexty)) ok=1;
					}
				else
					{
					if (nexty>=ymax+20)
						{
						if (noDefectInsideSquare(sol, x0,y1,nextx,nexty) && noDefectOnCut2(sol, nexty,sol->x0cut1,xnc1)) ok=1;
						}
					else
					if (nexty==ymax && ymax>=ymax_+TWENTY && noDefectInsideSquare(sol, x0,y1,nextx,nexty)) ok=1;
					}
				if (ok)
					{
					xd=xDefect(sol, nextx,nexty);
					if (xd>=xnc1)
						{
					sol->nbcut4[cut2]++;
					surf=surf+DATA->item_list_[j].surf;
					storeItemConf(x0,y1)
					config[i].rectangle=((int) 10000*sol->numplate)+((int) 100*sol->numcut1)+(int) cut2;
					config[i].ynextcut2=nexty;
					evalCut1(sol, i,surf,xnc1);
					if (xd<XUB) XUB=xd;
					if (evalCut2(sol, surf,xnc1,nexty)!=0)
						enumCut1(sol,cut2+1,i+1,sol->x0cut1,nexty,surf,xmax_,xmax,XUB,ymax_,nexty,PLATE_HEIGHT);
					if (evalCut3(sol,surf,nextx,xnc1,y0,nexty)!=0)
						enumCut1(sol,cut2,i+1,        nextx,   y0,surf,xmax_,xmax,XUB,ymax_,nexty,nexty);
					sol->nbcut4[cut2]--;
					}}
				}
			}
		}
}

void enumCut1(Solution* sol,int cut2,const int i,int x0,int y0,int surf,int xmax_,int xmax,int XUB,int ymax_,int ymax,int YUB)
	{
	int j,dx,dy,x,nextx,y,nexty,xd,yd,YUB_;
	int xmaxsave,xmax_save,ymaxsave,ymax_save,xnc1,ync2;
	if (i<DATA->nb_items_)
		{
		j = sol->seq_[i];
		dy=DATA->item_list_[j].width;
		dx=DATA->item_list_[j].length;
		surf=surf+DATA->item_list_[j].surf;
		
		nexty=PLATE_HEIGHT;
		nextx=x0+dx;
		if nextxCORRECT
			{
			y=nextY_Up(sol, x0,y0,dx,dy);
			nexty=y+dy;			
			if (y!=y0)
				{
				if nextyCORRECT
					EVALUATE_SHIFT_UP
				nexty=y0+dy;
				if nextyCORRECT
					{
					x=nextX_Right(sol, x0,y0,dx,dy);
					nextx=x+dx;
					if nextxCORRECT
						EVALUATE_SHIFT_RIGHT
					}
				}
			else
			if nextyCORRECT
			{
			  EVALUATE
			}
			}
		if (dx!=dy)
			{
			nextx=x0+dy;
			if nextxCORRECT
				{
				y=nextY_Up(sol, x0,y0,dy,dx);
				nexty=y+dx;
				if (y!=y0)
					{
					if nextyCORRECT
						EVALUATE_SHIFT_UP
					nexty=y0+dx;
					if nextyCORRECT
						{
						x=nextX_Right(sol, x0,y0,dy,dx);
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
	}

int evalPlate(Solution* sol,int i0)
{
	int i,x0prev;
	sol->x0cut1=sol->numcut1=0;
	sol->xmaxcut1=3500;
	while (!nodefectX[sol->xmaxcut1]) sol->xmaxcut1--;
	
	chkFirstItemOfPlate(sol, i0);
	
	sol->i0nextcut=i0;
	sol->x0nextcut=0;
	do	{

		i0prev=i0;
		sol->wastemin = INT_MAX;
		enumCut1(sol,0,i0, sol->x0cut1, 0, 0,
		0, 0, sol->xmaxcut1, 
		0, 0, PLATE_HEIGHT);

		if (sol->wastemin<INT_MAX)
			{
			for (i=i0prev;i<sol->i0nextcut;i++)
				{
				bestPlate[i]=bestCut1[i];
				bestPlate[i].plate=sol->numplate;
				bestPlate[i].xnextcut1=sol->x0nextcut;
				}
			x0prev=sol->x0cut1;
			sol->x0cut1=sol->x0nextcut;
			sol->xmaxcut1 = calcMaxXForNextCut1(sol, sol->x0cut1);
			i0prev=i0;
			i0=sol->i0nextcut;
			sol->numcut1++;
			}
		}
	while (sol->wastemin<INT_MAX && i0<DATA->nb_items_);
	return sol->i0nextcut;
}

std::string batchPath;
std::string defectsPath;
std::string optParamsPath;

int evalBatch(Solution* sol)
{
	int i,i0;
	i0=0;
	int np = 0;
	int xlast = 0;
	do {
		sol->numplate=np;
    	initDefectsPlateGenSol(sol, np);
		i=evalPlate(sol,i0);
		if (i<DATA->nb_items_)
		{
			i0=i;
			np++;
		}

		if(np > NMB_PLATES) {

			for(int i = 0; i < sol->seq_.size(); i++)
			{
				std::cout << sol->seq_[i] << " ";
			}
			std::cout << "\n";
			std::cout << "gensol problem\n";
			return -2;
		}

	} while(i<DATA->nb_items_);

	xlast=sol->x0nextcut;

	{
				int j,ymax;
				int rectangle_;

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
				{
					if (bestPlate[i].y1>bestPlate[i].ynextcut2) return -2;
				}

				if (make_solution_file(data, sol->solfile))
				{
					return 1;
				}

				return -2;
		}
}

}

