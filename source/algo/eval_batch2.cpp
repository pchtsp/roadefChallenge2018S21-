#include "cuttingGlass.h"

extern Data* DATA;
Solution* SOL_TO_EVALUATE1;


namespace evalBatch2
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

int MAX_CUT_SIZE = 3500;

int nbdp;
int Dx0[15];
int Dy0[15];
int Dx1[15];
int Dy1[15];
char nodefectX[PLATE_WIDTH + 1];

void initDefectsPlate(int pId)
	{
	int x,d;
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
int nextY_Up(int x,int y,int w,int h);
	
int isDefectInsideSquare(int x,int y,int x1,int y1)
	{
	int d;
	for (d=0;d<nbdp;d++)
		if (x<Dx1[d] && x1>Dx0[d] && y<Dy1[d] && y1>Dy0[d]) 
			{
			return 1;
			}
	return 0;
	}

int isDefectOnCUT2(const int y,const int x1)
	{
	int d;
	for (d=0;d<nbdp;d++)
		if (Dy1[d]>y && Dy0[d]<y && Dx1[d]>SOL_TO_EVALUATE1->x0cut1 && Dx0[d]<x1) 
			return 1;
	return 0;
	}

int shiftUp_CUT2_AfterDefect(int y,const int x1)
	{
	int d;
recom:
	for (d=0;d<nbdp;d++)
		if (Dy1[d]>y && Dy0[d]<y && Dx1[d]>SOL_TO_EVALUATE1->x0cut1 && Dx0[d]<x1) 
			{
			y=Dy1[d];
			goto recom;
			}
	return y;
	}

int shiftRight_CUT1_AfterDefect(int x)
	{
	int d;
recom:
	for (d=0;d<nbdp;d++)
		if (Dx1[d]>x && Dx0[d]<x) 
			{
			x=Dx1[d];
			goto recom;
			}
	return x;
	}

int nextY_Up(int x,int y,int w,int h)
	{
	int d,y0,x1;
	x1=x+w;
	y0=y;
recom:
	for (d=0;d<nbdp;d++)
		{
		if (x>=Dx1[d]) continue;
		if (y>=Dy1[d]) continue;
		if (x1<=Dx0[d]) continue;
		if (y+h<=Dy0[d]) continue;
		y=Dy1[d];
		if (y-y0<TWENTY) y=y0+TWENTY;
		if (y+h<PLATE_HEIGHT) goto recom;
		}
	return y;
	}

int nextX_Right(int x,int y,int w,int h)
	{
	int d,x0,y1;
	x0=x;
	y1=y+h;
recom:
	for (d=0;d<nbdp;d++)
		{
		if (x>=Dx1[d]) continue;
		if (y>=Dy1[d]) continue;
		if (x+w<=Dx0[d]) continue;
		if (y1<=Dy0[d]) continue;
		x=Dx1[d];
		if (x-x0<TWENTY) x=x0+TWENTY;
		if (x+w<PLATE_WIDTH) goto recom;
		}
	return x;
	}

/*	int x0cut2[33][35];
	int x1cut2[33][35];
	int y0cut2[33][35];
	int y1cut2[33][35];
	int nbicut2[35];
	*/

int nbcut3[33];
int xcut[33][35];
int x0cut[33][35];
int ycut[33][35];
int y0cut[33][35];
int ynextcut2[33];

void xmaxCUT1(int X0) 
	{
	SOL_TO_EVALUATE1->x0cut1=X0;
	SOL_TO_EVALUATE1->xmaxcut1=X0+MAX_CUT_SIZE;
	if (SOL_TO_EVALUATE1->xmaxcut1>=PLATE_WIDTH)
		SOL_TO_EVALUATE1->xmaxcut1=PLATE_WIDTH;
	else
		{
		if (SOL_TO_EVALUATE1->xmaxcut1>PLATE_WIDTH_20) 
			SOL_TO_EVALUATE1->xmaxcut1=PLATE_WIDTH_20;
		while (!nodefectX[SOL_TO_EVALUATE1->xmaxcut1]) SOL_TO_EVALUATE1->xmaxcut1--;
		}
	}


int isDefectOnCUT3(const int x,const int y0,const int y1)
	{
	int d;
	for (d=0;d<nbdp;d++)
		if (Dx1[d]>x && Dx0[d]<x && Dy1[d]>y0 && Dy0[d]<y1) 
			return 1;			
	return 0;
	}

int chkCUT2(int cut2,int xmax)
	{
	int i,y,ymax_,ymax,ync2_;			
	ymax_=0;
	ymax=ycut[cut2][0];
	for (i=1;i<nbcut3[cut2];i++) 
		{
		y=ycut[cut2][i];
		if (y>ymax) {ymax_=ymax;ymax=y;}
		else
		if (y>ymax_ && y!=ymax) ymax_=y;			
		}	
	if (ymax>PLATE_HEIGHT_20 && ymax<PLATE_HEIGHT) return 9999;
	if (ymax - ymax_<TWENTY) 
		{
		ymax+=TWENTY;
		if (ymax>=PLATE_HEIGHT) return ymax;
		ymax=shiftUp_CUT2_AfterDefect(ymax,xmax);
		if (ymax>PLATE_HEIGHT_20) ymax=PLATE_HEIGHT;
		
		ync2_=cut2!=0?ynextcut2[cut2-1]:0;
		for (i=0;i<nbcut3[cut2];i++) 
			{
			if (isDefectOnCUT3(xcut[cut2][i],ymax_,ymax)) return 9999;
			if (x0cut[cut2][i]!=xcut[cut2][i])
				{
				if (y0cut[cut2][i]!=ync2_ && nextY_Up(x0cut[cut2][i],ycut[cut2][i],
					xcut[cut2][i]-x0cut[cut2][i],
					ymax-ycut[cut2][i])!=ycut[cut2][i]) return 9999;
				
				//isDefectInsideSquare(x0cut[cut2][j],ycut[cut2][j],xcut[cut2][j],ymax)) return 9999;
				}
			}
		return ymax;
		}
	y=shiftUp_CUT2_AfterDefect(ymax,xmax);
	if (y!=ymax) 
		{
		if (y<ymax+TWENTY)
			{
			ymax+=TWENTY;
			if (ymax>PLATE_HEIGHT) return ymax;
			ymax=shiftUp_CUT2_AfterDefect(ymax,xmax);
			}
		else
			ymax=y;
		}
	if (ymax>PLATE_HEIGHT_20) ymax=PLATE_HEIGHT;
	ync2_=cut2!=0?ynextcut2[cut2-1]:0;
	for (i=0;i<nbcut3[cut2];i++) 
		{
		if (isDefectOnCUT3(xcut[cut2][i],ymax_,ymax)) return 9999;
		if (x0cut[cut2][i]!=xcut[cut2][i])
			{
			//if (y0cut[cut2][i]!=ync2_ && isDefectInsideSquare(x0cut[cut2][i],ycut[cut2][i],xcut[cut2][i],ymax)) return 9999;
			if (y0cut[cut2][i]!=ync2_ && nextY_Up(x0cut[cut2][i],ycut[cut2][i],xcut[cut2][i]-x0cut[cut2][i],ymax-ycut[cut2][i])!=ycut[cut2][i]) return 9999;
			}
		}
	return ymax;
	}

int chkAllCUT2s(int cut2,const int x1cut1)
	{
	int i;			
	for (i=0;i<=cut2;i++) 
		if (isDefectOnCUT2(ynextcut2[i],x1cut1)) return 0;
	return 1;
	}
	
int chkCUT1(int cut2,int xmax)
	{
	int i,x,xmax_;	
	if (xmax>PLATE_WIDTH_20 && xmax<PLATE_WIDTH) return 9999;
	xmax_=0;
	for (i=0;i<=cut2;i++) 
		{
		x=xcut[i][nbcut3[i]-1];
		if (x!=xmax && x>xmax_) xmax_=x;			
		}
	if (xmax<xmax_+TWENTY) 
		{
		xmax+=TWENTY;//printf("(%d,%d)",cut2,xmax);
		if (xmax>=PLATE_WIDTH) return xmax;
		xmax=shiftRight_CUT1_AfterDefect(xmax);
		if (xmax>PLATE_WIDTH_20) return PLATE_WIDTH;
		return xmax;
		}
	x=shiftRight_CUT1_AfterDefect(xmax);
	if (x!=xmax)
		{
		if (x - xmax<TWENTY) 
			{
			xmax+=TWENTY;
			if (xmax>=PLATE_WIDTH) return xmax;
			xmax=shiftRight_CUT1_AfterDefect(xmax);
			if (xmax>PLATE_WIDTH_20) return PLATE_WIDTH;
			return xmax;
			}
		else 
			{
			if (x>PLATE_WIDTH_20) return PLATE_WIDTH;
			return x;
			}
		}
	return xmax;
	}

#define EVALUATE {\
	xmaxsave=xmax;if (nextx>xmax) xmax=nextx;\
	x0cut[cut2][nbcut3[cut2]]=x0;\
	xcut[cut2][nbcut3[cut2]]=nextx;\
	ycut[cut2][nbcut3[cut2]]=nexty;\
	y0cut[cut2][nbcut3[cut2]]=y;\
	nbcut3[cut2]++;\
	xnc1=chkCUT1(cut2,xmax);\
	if (xnc1<xmin)\
		{\
		ync2=chkCUT2(cut2,xnc1);\
		if (ync2<=PLATE_HEIGHT_20 || ync2==PLATE_HEIGHT)\
			{\
			ync2_=ynextcut2[cut2];ynextcut2[cut2]=ync2;\
			if (i<SOL_TO_EVALUATE1->to_i1)\
				{\
				szmin=DATA->item_list_[SOL_TO_EVALUATE1->seq_[i+1]].length;\
				if (nexty+szmin<=PLATE_HEIGHT) tryCut4(cut2,i+1,x0,y0,nextx,xnc1,nexty,xmax,xmin);\
				if (ync2+szmin<=PLATE_HEIGHT) enumCut1(cut2+1,i+1,SOL_TO_EVALUATE1->x0cut1,ync2,xmax,xmin);\
				if (nextx+szmin<xmin && nbcut3[cut2]<SOL_TO_EVALUATE1->params_->MAXCUT3) enumCut1(cut2,i+1,nextx,y0,xmax,xmin);\
				}\
			else\
			if (chkAllCUT2s(cut2,xnc1)) xmin=xnc1;\
			ynextcut2[cut2]=ync2_;\
			}\
		}\
	nbcut3[cut2]--;\
	xmax=xmaxsave;\
	}

#define EVALUATE_SHIFT_UP {\
	xmaxsave=xmax;if (nextx>xmax) xmax=nextx;\
	x0cut[cut2][nbcut3[cut2]]=x0;\
	xcut[cut2][nbcut3[cut2]]=nextx;\
	ycut[cut2][nbcut3[cut2]]=nexty;\
	y0cut[cut2][nbcut3[cut2]]=y;\
	nbcut3[cut2]++;\
	xnc1=chkCUT1(cut2,xmax);\
	if (xnc1<xmin)\
		{\
		ync2=chkCUT2(cut2,xnc1);\
		if (ync2<=PLATE_HEIGHT_20 || ync2==PLATE_HEIGHT)\
			{\
			ync2_=ynextcut2[cut2];ynextcut2[cut2]=ync2;\
			if (i<SOL_TO_EVALUATE1->to_i1)\
				{\
				szmin=DATA->item_list_[SOL_TO_EVALUATE1->seq_[i+1]].length;\
				if (ync2+szmin<=PLATE_HEIGHT) enumCut1(cut2+1,i+1,SOL_TO_EVALUATE1->x0cut1,ync2,xmax,xmin);\
				if (nextx+szmin<xmin &&  nbcut3[cut2]<SOL_TO_EVALUATE1->params_->MAXCUT3) enumCut1(cut2,i+1,nextx,y0,xmax,xmin);\
				}\
			else\
			if (chkAllCUT2s(cut2,xnc1)) xmin=xnc1;\
			ynextcut2[cut2]=ync2_;\
			}\
		}\
	nbcut3[cut2]--;\
	xmax=xmaxsave;\
	}



#define EVALUATE_SHIFT_RIGHT {\
	xmaxsave=xmax;if (nextx>xmax) xmax=nextx;\
	x0cut[cut2][nbcut3[cut2]]=x;\
	xcut[cut2][nbcut3[cut2]]=x;\
	ycut[cut2][nbcut3[cut2]]=nexty;\
	y0cut[cut2][nbcut3[cut2]]=y;\
	nbcut3[cut2]++;\
	x0cut[cut2][nbcut3[cut2]]=x;\
	xcut[cut2][nbcut3[cut2]]=nextx;\
	ycut[cut2][nbcut3[cut2]]=nexty;\
	y0cut[cut2][nbcut3[cut2]]=y;\
	nbcut3[cut2]++;\
	xnc1=chkCUT1(cut2,xmax);\
	if (xnc1<xmin)\
		{\
		ync2=chkCUT2(cut2,xnc1);\
		if (ync2<=PLATE_HEIGHT_20 || ync2==PLATE_HEIGHT)\
			{\
			ync2_=ynextcut2[cut2];ynextcut2[cut2]=ync2;\
			if (i<SOL_TO_EVALUATE1->to_i1)\
				{\
				szmin=DATA->item_list_[SOL_TO_EVALUATE1->seq_[i+1]].length;\
				if (nexty+szmin<=PLATE_HEIGHT) tryCut4(cut2,i+1,x,y0,nextx,xnc1,nexty,xmax,xmin);\
				if (ync2+szmin<=PLATE_HEIGHT) enumCut1(cut2+1,i+1,SOL_TO_EVALUATE1->x0cut1,ync2,xmax,xmin);\
				if (nextx+szmin<xmin && nbcut3[cut2]<=SOL_TO_EVALUATE1->params_->MAXCUT3) enumCut1(cut2,i+1,nextx,y0,xmax,xmin);\
				}\
			else\
			if (chkAllCUT2s(cut2,xnc1)) xmin=xnc1;\
			ynextcut2[cut2]=ync2_;\
			}\
		}\
	nbcut3[cut2]--;\
	nbcut3[cut2]--;\
	xmax=xmaxsave;\
	}

#define EVALUATE_CUT4 {\
	x0cut[cut2][nbcut3[cut2]]=x0;\
	xcut[cut2][nbcut3[cut2]]=nextx;\
	ycut[cut2][nbcut3[cut2]]=nexty;\
	y0cut[cut2][nbcut3[cut2]]=y1;\
	nbcut3[cut2]++;\
	ync2=chkCUT2(cut2,xnc1);\
	if (nexty==ync2 && (ync2<=PLATE_HEIGHT_20 || ync2==PLATE_HEIGHT) &&\
			(SOL_TO_EVALUATE1->nbcut4[cut2]==0 || nexty==SOL_TO_EVALUATE1->ycut4[cut2]))\
		{\
		ync2_=ynextcut2[cut2];ynextcut2[cut2]=ync2;\
		nexty_=SOL_TO_EVALUATE1->ycut4[cut2];SOL_TO_EVALUATE1->ycut4[cut2]=nexty;\
		SOL_TO_EVALUATE1->nbcut4[cut2]++;\
		if (i<SOL_TO_EVALUATE1->to_i1)\
			{\
			if (ync2+DATA->item_list_[SOL_TO_EVALUATE1->seq_[i+1]].length<=PLATE_HEIGHT)\
				enumCut1(cut2+1,i+1,SOL_TO_EVALUATE1->x0cut1,ync2,xmax,xmin);\
			if (nextx+DATA->item_list_[SOL_TO_EVALUATE1->seq_[i+1]].length<xmin)\
				enumCut1(cut2,i+1,nextx,y0,xmax,xmin);\
			}\
		else\
		if (chkAllCUT2s(cut2,xnc1)) xmin=xnc1;\
		SOL_TO_EVALUATE1->nbcut4[cut2]--;\
		SOL_TO_EVALUATE1->ycut4[cut2]=nexty_;\
		ynextcut2[cut2]=ync2_;\
		}\
	nbcut3[cut2]--;\
	}
		
void enumCut1(const int cut2,const int i,int x0,int y0,int xmax,int &xmin);
void tryCut4(const int cut2, const int i,int x0,int y0,const int nextx,const int xnc1,int y1,int xmax,int &xmin)
	{
	int j,dx,dy,nexty,nexty_,ync2,ync2_;
	j=SOL_TO_EVALUATE1->seq_[i];
	dy=DATA->item_list_[j].width;
	dx=DATA->item_list_[j].length;	
	if (dx==nextx-x0)
		{
		nexty=y1+dy;
		if ((nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT) &&
			(SOL_TO_EVALUATE1->nbcut4[cut2]==0 || nexty==SOL_TO_EVALUATE1->ycut4[cut2]) && 
			isDefectInsideSquare(x0,y1,nextx,nexty)==0) 
				EVALUATE_CUT4
		return;
		}
	if (dy==nextx-x0)
		{
		nexty=y1+dx;
		if ((nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT) &&
			(SOL_TO_EVALUATE1->nbcut4[cut2]==0 || nexty==SOL_TO_EVALUATE1->ycut4[cut2]) && 
			isDefectInsideSquare(x0,y1,nextx,nexty)==0)
				EVALUATE_CUT4
		}
	}

void enumCut1(const int cut2,const int i,int x0,int y0, int xmax,int &xmin)
	{
	int j,dx,dy,x,nextx,y,nexty,xmaxsave,szmin;
	int xnc1,ync2_,ync2;
	
	
	j=SOL_TO_EVALUATE1->seq_[i];
	dy=DATA->item_list_[j].width;
	dx=DATA->item_list_[j].length;	
	nextx=x0+dx;
	if (nextx<xmin)
		{
		y=nextY_Up(x0,y0,dx,dy);
		nexty=y+dy;
		if ((nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT) && (SOL_TO_EVALUATE1->nbcut4[cut2]==0 || 
			nexty<=SOL_TO_EVALUATE1->ycut4[cut2]-TWENTY || nexty==SOL_TO_EVALUATE1->ycut4[cut2]))
			{
			if (y!=y0)
				EVALUATE_SHIFT_UP
			else
				EVALUATE
			}
		if (y!=y0)
			{
			x=nextX_Right(x0,y0,dx,dy);
			nextx=x+dx;
			if (nextx<xmin)
				{
				nexty=y0+dy;
				if ((nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT) && (SOL_TO_EVALUATE1->nbcut4[cut2]==0 || 
					nexty<=SOL_TO_EVALUATE1->ycut4[cut2]-TWENTY || nexty==SOL_TO_EVALUATE1->ycut4[cut2]))
					EVALUATE_SHIFT_RIGHT
				}
			}
		}		
	if (dx!=dy)
		{
		nextx=x0+dy;
		if (nextx<xmin)
			{
			y=nextY_Up(x0,y0,dy,dx);
			nexty=y+dx;
			if ((nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT) && (SOL_TO_EVALUATE1->nbcut4[cut2]==0 || 
				nexty<=SOL_TO_EVALUATE1->ycut4[cut2]-TWENTY || nexty==SOL_TO_EVALUATE1->ycut4[cut2]))
				{
				if (y!=y0)
					EVALUATE_SHIFT_UP
				else
					EVALUATE
				}
			if (y!=y0)
				{
				x=nextX_Right(x0,y0,dy,dx);
				nextx=x+dy;
				if (nextx<xmin) 
					{
					nexty=y0+dx;
					if ((nexty<=PLATE_HEIGHT_20 || nexty==PLATE_HEIGHT) && (SOL_TO_EVALUATE1->nbcut4[cut2]==0 || 
						nexty<=SOL_TO_EVALUATE1->ycut4[cut2]-TWENTY || nexty==SOL_TO_EVALUATE1->ycut4[cut2]))
						EVALUATE_SHIFT_RIGHT
					}
				}					
			}
		}
	}

//int totalsurf;
int x0nextcut;

int plateItem[MAX_ITEMS],i0Item[MAX_ITEMS],xminItem[MAX_ITEMS],maxItemPlate;
int surfItem[MAX_ITEMS];
#define OUT_OF_PLATE 6001

int evalCUT1(int maxItemPlate)
	{
	int i1,xmin;	
	i1=SOL_TO_EVALUATE1->to_i1=SOL_TO_EVALUATE1->from_i0;
	do {	
		SOL_TO_EVALUATE1->to_i1=i1;
		if (xminItem[i1]<OUT_OF_PLATE) 
			{
			xmin=xminItem[i1];
			enumCut1(0,SOL_TO_EVALUATE1->from_i0,SOL_TO_EVALUATE1->x0cut1,0,0,xmin);	
			if (xmin<xminItem[i1])
				{
				plateItem[i1]=SOL_TO_EVALUATE1->numplate;
				xminItem[i1]=xmin;
				i0Item[i1]=SOL_TO_EVALUATE1->from_i0;
				}
			i1++;
			}
		else
			{
			xmin=SOL_TO_EVALUATE1->xmaxcut1+1;
			enumCut1(0,SOL_TO_EVALUATE1->from_i0,SOL_TO_EVALUATE1->x0cut1,0,0,xmin);	
			if (xmin<=SOL_TO_EVALUATE1->xmaxcut1)
				{
				if (xmin<xminItem[i1])
					{
					plateItem[i1]=SOL_TO_EVALUATE1->numplate;
					xminItem[i1]=xmin;
					i0Item[i1]=SOL_TO_EVALUATE1->from_i0;
					}
				i1++;
				}
			}







/*

		enumCut1(0,from_i0,i1,SOL_TO_EVALUATE1->x0cut1,0,0,xmin);	
		if (xmin<=SOL_TO_EVALUATE1->xmaxcut1)
			{
			if (xmin<xminItem[i1])
				{
				xminItem[i1]=xmin;
				i0Item[i1]=from_i0;
				}
				xminItem[i1]=xmin;
			i1++;
			}*/
		} 
	while (xmin<=SOL_TO_EVALUATE1->xmaxcut1 && i1<maxItemPlate);
	return i1;
	}



std::vector<int> evalPlate(Solution* sol, int i0)
{
	SOL_TO_EVALUATE1 = sol;
	std::vector<int> res;
	int i0_,i1,imax;
	surfItem[i0]=DATA->item_list_[SOL_TO_EVALUATE1->seq_[i0]].surf;
	xminItem[i0]=OUT_OF_PLATE; 	
	for (maxItemPlate=i0+1;maxItemPlate<DATA->nb_items_;maxItemPlate++)
		{
		surfItem[maxItemPlate]=surfItem[maxItemPlate-1]+DATA->item_list_[SOL_TO_EVALUATE1->seq_[maxItemPlate]].surf;
		xminItem[maxItemPlate]=OUT_OF_PLATE; 	
		if (surfItem[maxItemPlate]>PLATE_SURFACE) break;
		}

	imax=i0_=i0;
	xmaxCUT1(0);

tryCut1From_i0:
	SOL_TO_EVALUATE1->from_i0=i0;
	i1=evalCUT1(maxItemPlate);
	if (i1>imax) imax=i1;
	if (i1<=i0_)
		{
		x0nextcut=0;
		res.push_back(0);
		res.push_back(i0_);
		res.push_back(0);
		return res;
		}
nexti0:
		i0++;
		if (i0<=imax && i0<maxItemPlate)
			{
			int x=xminItem[i0-1];
			if (x>=OUT_OF_PLATE) 
				{
				res.push_back(-1);
				res.push_back(-1);
				res.push_back(-1);
				return res;			
				}
			if (x+DATA->item_list_[SOL_TO_EVALUATE1->seq_[i0]].length<=PLATE_WIDTH)
				{
				xmaxCUT1(x);
				goto tryCut1From_i0;
				}
			goto nexti0;
			}
	x0nextcut=xminItem[imax-1];
	res.push_back(surfItem[imax-1]);
	res.push_back(imax);
	res.push_back(xminItem[imax-1]);
	return res;
	}	

int evalBatch(Solution* sol, int timeLimit)
{
	SOL_TO_EVALUATE1 = sol;
	std::vector<int> eval_plate_result;
	int np,i, i0;
	i0 = np = 0;
	
	int t0 = time(0);
	
	do
	{
		sol->numplate=np;
		initDefectsPlate(np);				

		eval_plate_result = evalPlate(sol, i0);
		if (eval_plate_result[0]<=0) return -1;
		
		i = eval_plate_result[1];
		
		if(i < DATA->nb_items_)
		{
			i0 = i;
			np++;
			if (np>=NMB_PLATES) return -1;
		}

		if(time(0) - t0 >= timeLimit) return -2;

	}
	while(i < DATA->nb_items_);


	return (PLATE_WIDTH * np) + eval_plate_result[2];


}

	
}
