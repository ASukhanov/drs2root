
//Copy filter coeffs into array
double fcx[100],fcy[100];
double *xxxxxx,*yyyyyy; int shift;
int fcn=d4r->fgraph->GetN();
xxxxxx=d4r->fgraph->GetX(); yyyyyy=d4r->fgraph->GetY();
for(ii=0;ii<fcn;ii++) {fcx[ii]=xxxxxx[ii];fcy[ii]=yyyyyy[ii];}

//Shift TGrapf, double *xx,*yy; int shift; 
shift=26.; for(ii=0;ii<fcn;ii++) fcx[ii]+=shift;TGraph tg(d4r->fgraph->GetN(),fcx,fcy);
