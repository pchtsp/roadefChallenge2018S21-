#include "cuttingGlass.h"


namespace evalBatch0
{

Data* DATA;

#include "enumCutTools.cpp"

void evalCut1(Solution* sol, const int &i,int &surf,const int &xEndCut1)
{
	int waste = ((xEndCut1 - sol->x0cut1) * PLATE_HEIGHT) - surf;
	if (waste < sol->wastemin || (waste==sol->wastemin && surf>sol->surfCut1))
	{
		sol->wastemin = waste;
		sol->i0nextcut = i + 1;
		sol->x0nextcut = xEndCut1;
		sol->surfCut1 = surf;
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
				if(nexty <= ync2 - 20) rotat = 0;\
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

#define nextxCORRECT (nextx<=XUB && (nextx<=PLATE_WIDTH_20 || nextx==PLATE_WIDTH))
#define nextyCORRECT (nexty<=YUB && (nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT))

void enumCut1(Solution* sol,int cut2,const int i,int x0,int y0,int surf,int xmax_,int xmax,int XUB,int ymax_,int ymax,int YUB);
void tryCut4(Solution* sol,int cut2, const int i,int x0,int y0,int x1,int y1,int &xnc1,int &ync2,int surf,int& xmax_,int& xmax,int XUB,int& ymax_,int& ymax,int YUB)
	{
	int j,dx,dy,nexty,ok,xd;
	if (i<DATA->nb_items_)
		{
		j = sol->seq_[i];
		dy=DATA->item_list_[j].width;
		dx=DATA->item_list_[j].length;
		if (dx==x1-x0)
			{
			nexty=y1+dy;
			if nextyCORRECT
				{
				ok=0;
				if (sol->nbcut4[cut2]!=0)
					{
					if (nexty==ymax && noDefectInsideSquare(sol, x0,y1,x1,nexty)) ok=1;
					}
				else
				{
				if (nexty>=ymax+20)
					{
					if (noDefectInsideSquare(sol, x0,y1,x1,nexty) && noDefectOnCut2(sol, nexty,sol->x0cut1,xnc1)) ok=1;
					}
				else
				if (nexty==ymax && ymax>=ymax_+TWENTY && noDefectInsideSquare(sol, x0,y1,x1,nexty)) ok=1;
				}
				if (ok)
					{
					xd=xDefect(sol, x1,nexty);
					if (xd>=xnc1)
						{
						sol->nbcut4[cut2]++;
						surf=surf+DATA->item_list_[j].surf;
						evalCut1(sol, i,surf,xnc1);

						if (xd<XUB) XUB=xd;

						if (evalCut2(sol, surf,xnc1,nexty)!=0)
							enumCut1(sol,cut2+1,i+1,sol->x0cut1,nexty,surf,xmax_,xmax,XUB,ymax_,nexty,PLATE_HEIGHT);
						if (evalCut3(sol,surf,x1,xnc1,y0,nexty)!=0)
							enumCut1(sol,cut2,i+1,         x1,   y0,surf,xmax_,xmax,XUB,ymax_,nexty,nexty);
						sol->nbcut4[cut2]--;
						}
					}
				}
			return;
			}
		if (dy==x1-x0)
			{
			nexty=y1+dx;
			if nextyCORRECT
				{
				ok=0;
				if (sol->nbcut4[cut2]!=0)
					{
					if (nexty==ymax && noDefectInsideSquare(sol, x0,y1,x1,nexty)) ok=1;
					}
				else
				{
				if (nexty>=ymax+20)
					{
					if (noDefectInsideSquare(sol, x0,y1,x1,nexty) && noDefectOnCut2(sol, nexty,sol->x0cut1,xnc1)) ok=1;
					}
				else
				if (nexty==ymax && ymax>=ymax_+TWENTY && noDefectInsideSquare(sol, x0,y1,x1,nexty)) ok=1;
				}
				if (ok)
					{
					xd=xDefect(sol, x1,nexty);
					if (xd>=xnc1)
						{
						sol->nbcut4[cut2]++;
						surf=surf+DATA->item_list_[j].surf;
						evalCut1(sol, i,surf,xnc1);

						if (xd<XUB) XUB=xd;
						if (evalCut2(sol, surf,xnc1,nexty)!=0)
							enumCut1(sol,cut2+1,i+1,sol->x0cut1,nexty,surf,xmax_,xmax,XUB,ymax_,nexty,PLATE_HEIGHT);
						if (evalCut3(sol,surf,x1,xnc1,y0,nexty)!=0)
							enumCut1(sol,cut2,i+1,         x1,   y0,surf,xmax_,xmax,XUB,ymax_,nexty,nexty);
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
	int xmaxsave,xmax_save,ymaxsave,ymax_save,xnc1,ync2,rotat;
	if (i<DATA->nb_items_)
		{
		j = sol->seq_[i];
		dy=DATA->item_list_[j].width;
		dx=DATA->item_list_[j].length;
		surf=surf+DATA->item_list_[j].surf;

		nexty=PLATE_HEIGHT;
		nextx=x0+dx;
		rotat=1;
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
				EVALUATE
			}
		if (/*rotat && */dx!=dy)
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

#define EVALUATE0 {\
		xmaxsave=xmax;xmax_save=xmax_;\
		if (nextx>xmax) {xmax_=xmax;xmax=nextx;}else if (nextx!=xmax && nextx>xmax_) xmax_=nextx;\
		xnc1=chkVt0(sol, xmax_,xmax);\
		if (xnc1<=sol->xmaxcut1)\
			{\
			ymaxsave=ymax;ymax_save=ymax_;\
			if (nexty>ymax) {ymax_=ymax;ymax=nexty;}else if (nexty!=ymax && nexty>ymax_) ymax_=nexty;\
			ync2=chkHz0(sol, ymax_,ymax);\
			if (ync2<=YUB)\
				{\
				if(nexty <= ync2 - 20) rotat = 0;\
				evalCut1(sol, i,surf,xnc1);\
				if (evalCut2(sol, surf,xnc1,ync2)!=0)\
					enumCut1_0(sol,cut2+1,i+1,sol->x0cut1,ync2,surf,xmax_,xmax,ymax_,ymax,PLATE_HEIGHT);\
				if (evalCut3(sol, surf,nextx,xnc1,y0,ync2)!=0)\
					enumCut1_0(sol,cut2,i+1,nextx,y0,surf,xmax_,xmax,ymax_,ymax,YUB);\
				tryCut4_0(sol,cut2,i+1,x0,y0,nextx,nexty,xnc1,ync2,surf,xmax_,xmax,ymax_,ymax,YUB);\
				}\
			ymax=ymaxsave;ymax_=ymax_save;\
			}\
		xmax=xmaxsave;xmax_=xmax_save;\
		}


void enumCut1_0(Solution* sol,int cut2,const int i,int x0,int y0,int surf,int xmax_,int xmax,int ymax_,int ymax,int YUB);

void tryCut4_0(Solution* sol,int cut2, const int i,int x0,int y0,int x1,int y1,int &xnc1,int &ync2,int surf,int& xmax_,int& xmax,int& ymax_,int& ymax,int YUB)
	{
	int j,dx,dy,nexty,ok,xd;
	if (i<DATA->nb_items_)
		{
		j = sol->seq_[i];
		dy=DATA->item_list_[j].width;
		dx=DATA->item_list_[j].length;
		if (dx==x1-x0)
			{
			nexty=y1+dy;
			if nextyCORRECT
				{
				ok=0;
				if (sol->nbcut4[cut2]!=0)
				{
					if (nexty==ymax) ok=1;
				}
				else
				if (nexty>=ymax+20) ok=1;
				else
				if (nexty==ymax && ymax>=ymax_+TWENTY) ok=1;
				if (ok)
					{
					sol->nbcut4[cut2]++;
					surf=surf+DATA->item_list_[j].surf;
					evalCut1(sol, i,surf,xnc1);
					if (evalCut2(sol, surf,xnc1,nexty)!=0)
						enumCut1_0(sol,cut2+1,i+1,sol->x0cut1,nexty,surf,xmax_,xmax,ymax_,nexty,PLATE_HEIGHT);
					if (evalCut3(sol,surf,x1,xnc1,y0,nexty)!=0)
						enumCut1_0(sol,cut2,i+1,         x1,   y0,surf,xmax_,xmax,ymax_,nexty,nexty);
					sol->nbcut4[cut2]--;
					}
				}
			return;
			}
		if (dy==x1-x0)
			{
			nexty=y1+dx;
			if nextyCORRECT
				{
				ok=0;
				if (sol->nbcut4[cut2]!=0)
				{
					if (nexty==ymax) ok=1;
				}
				else
				if (nexty>=ymax+20) ok=1;
				else
				if (nexty==ymax && ymax>=ymax_+TWENTY) ok=1;
				if (ok)
					{
					sol->nbcut4[cut2]++;
					surf=surf+DATA->item_list_[j].surf;
					evalCut1(sol, i,surf,xnc1);
					if (evalCut2(sol, surf,xnc1,nexty)!=0)
						enumCut1_0(sol,cut2+1,i+1,sol->x0cut1,nexty,surf,xmax_,xmax,ymax_,nexty,PLATE_HEIGHT);
					if (evalCut3(sol,surf,x1,xnc1,y0,nexty)!=0)
						enumCut1_0(sol,cut2,i+1,         x1,   y0,surf,xmax_,xmax,ymax_,nexty,nexty);
					sol->nbcut4[cut2]--;
					}
				}
			}
		}
}

void enumCut1_0(Solution* sol,int cut2,const int i,int x0,int y0,int surf,int xmax_,int xmax,int ymax_,int ymax,int YUB)
	{
	int j,dx,dy,x,nextx,y,nexty,xd,yd,YUB_;
	int xmaxsave,xmax_save,ymaxsave,ymax_save,xnc1,ync2, rotat;
	if (i<DATA->nb_items_)
		{
		j = sol->seq_[i];
		dy=DATA->item_list_[j].width;
		dx=DATA->item_list_[j].length;
		surf=surf+DATA->item_list_[j].surf;
		nexty=PLATE_HEIGHT;
		nextx=x0+dx;
		rotat=1;
		if (nextx<=sol->xmaxcut1 && (nextx<=PLATE_WIDTH_20 || nextx==PLATE_WIDTH))
			{
			nexty=y0+dy;
			if nextyCORRECT
			{
				EVALUATE0
			}
			}
		if (/*rotat && */dx!=dy)
			{
			nextx=x0+dy;
			if (nextx<=sol->xmaxcut1 && (nextx<=PLATE_WIDTH_20 || nextx==PLATE_WIDTH))
				{
				nexty=y0+dx;
				if nextyCORRECT
					EVALUATE0
				}
			}
		}
	}



// returns total surface of items packed, index in the sequence for next plate start and last cut_1 position
std::vector<int> evalPlate(Solution* sol, int i0)
{
	std::vector<int> res;
	sol->x0cut1 =  0;
	sol->xmaxcut1 = 3500;
	while (!nodefectX[sol->xmaxcut1]) sol->xmaxcut1--;

	chkFirstItemOfPlate(sol, i0);


	int totalsurf = 0;
	sol->i0nextcut  = i0;
	sol->x0nextcut = 0;
	do	{
		sol->wastemin = INT_MAX;

		initDefectsPlate(sol, sol->x0cut1, sol->xmaxcut1);


		if(sol->nbdp)
		{
			enumCut1(sol,0,i0, sol->x0cut1, 0, 0,
			0, 0, sol->xmaxcut1,
			0, 0, PLATE_HEIGHT);
		}
		else
		{
			enumCut1_0(sol,0,i0, sol->x0cut1, 0, 0,
			0, 0,
			0, 0, PLATE_HEIGHT);
		}


		if(sol->wastemin < INT_MAX)
		{
			sol->x0cut1 = sol->x0nextcut;
			sol->xmaxcut1 = calcMaxXForNextCut1(sol, sol->x0cut1);
			i0 = sol->i0nextcut;
			totalsurf += sol->surfCut1;
			sol->numcut1++;
		}
	}
	while (sol->wastemin < INT_MAX && i0 < DATA->nb_items_);

	res.push_back(totalsurf);
	res.push_back(sol->i0nextcut);
	res.push_back(sol->x0nextcut);
	return res;
}

int evalBatch(Solution* sol, bool print)
{
	int i, i0;
	i0 = 0;
	int np = 0;
	double totalW = 0;
	do
	{
		initDefectsPlate(np);

		std::vector<int> eval_plate_result = evalPlate(sol, i0);
		i = eval_plate_result[1];


		if(print)
		{
			  double W = (100.0 * (PLATE_SURFACE -  eval_plate_result[0])) / PLATE_SURFACE;
			  totalW += W;
			  std::cout << std::setw(5) << np
			  << std::setw(15) << eval_plate_result[0]
			  << std::setw(15) << eval_plate_result[1]
			  << std::setw(15) << eval_plate_result[2]
			  << std::setw(15) << W
			  << std::setw(15) << totalW / (np + 1)
			  << std::endl;
		}

		if(i < DATA->nb_items_)
		{
			i0 = i;
			np++;
		}
	}
	while(i < DATA->nb_items_);

	return (PLATE_WIDTH * np) + sol->x0nextcut;
}

std::vector<int> evalPlate_to_Mark(Solution* sol, int i0)
{
	std::vector<int> res;
	sol->numcut1 = 0;
	sol->x0cut1  = 0;
	sol->xmaxcut1 = 3500;
	while (!nodefectX[sol->xmaxcut1]) sol->xmaxcut1--;
	int totalsurf = 0;
	sol->i0nextcut  = i0;
	sol->x0nextcut = 0;
	do	{
		sol->wastemin = INT_MAX;


		initDefectsPlate(sol, sol->x0cut1, sol->xmaxcut1);

		if(sol->nbdp)
		{
			enumCut1(sol,0,i0, sol->x0cut1, 0, 0,
			0, 0, sol->xmaxcut1,
			0, 0, PLATE_HEIGHT);
		}
		else
		{
			enumCut1_0(sol,0,i0, sol->x0cut1, 0, 0,
			0, 0,
			0, 0, PLATE_HEIGHT);
		}

		if(sol->wastemin < INT_MAX)
			{
			for (int i=i0;i<sol->i0nextcut;i++)
				sol->cut1OfIdx[i]=sol->numcut1;
			sol->tcut1[sol->numcut1].i0cut1=i0;
			sol->tcut1[sol->numcut1].surf=sol->surfCut1;
			sol->tcut1[sol->numcut1].x0cut1=sol->x0cut1;

			sol->x0cut1 = sol->x0nextcut;
			sol->xmaxcut1 = calcMaxXForNextCut1(sol, sol->x0cut1);
			i0 = sol->i0nextcut;
			totalsurf += sol->surfCut1;
			sol->numcut1++;
			}
	}
	while (sol->wastemin < INT_MAX && i0 < DATA->nb_items_);

	sol->tcut1[sol->numcut1].i0cut1=i0;

	res.push_back(totalsurf);
	res.push_back(sol->i0nextcut);
	res.push_back(sol->x0nextcut);
	return res;
}

std::vector<int> evalPlateFromCut1(Solution* sol, int i0)
{
	std::vector<int> res;
	sol->numcut1=sol->cut1OfIdx[i0];
	if (sol->numcut1)
	{
		i0 = sol->tcut1[sol->numcut1].i0cut1 - 1;
		sol->numcut1=sol->cut1OfIdx[i0];
	}

	i0=sol->tcut1[sol->numcut1].i0cut1;

	sol->x0cut1 = sol->tcut1[sol->numcut1].x0cut1;
	sol->xmaxcut1 = calcMaxXForNextCut1(sol, sol->x0cut1);

	int totalsurf = 0;

	for (int nc=0;nc<sol->numcut1;nc++) totalsurf+=sol->tcut1[nc].surf;

	sol->i0nextcut  = i0;
	do	{
		sol->wastemin = INT_MAX;


		initDefectsPlate(sol, sol->x0cut1, sol->xmaxcut1);


		if(sol->nbdp)
		{
			enumCut1(sol,0,i0, sol->x0cut1, 0, 0,
			0, 0, sol->xmaxcut1,
			0, 0, PLATE_HEIGHT);
		}
		else
		{
			enumCut1_0(sol,0,i0, sol->x0cut1, 0, 0,
			0, 0,
			0, 0, PLATE_HEIGHT);
		}



		if(sol->wastemin < INT_MAX)
			{
			sol->x0cut1 = sol->x0nextcut;
			sol->xmaxcut1 = calcMaxXForNextCut1(sol, sol->x0cut1);
			i0 = sol->i0nextcut;
			totalsurf += sol->surfCut1;
			}
	}
	while (sol->wastemin < INT_MAX && i0 < DATA->nb_items_);

	res.push_back(totalsurf);
	res.push_back(sol->i0nextcut);
	res.push_back(sol->x0nextcut);
	return res;
}




}
