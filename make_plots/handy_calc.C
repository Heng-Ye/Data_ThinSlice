#include "../headers/BasicParameters.h"
#include "../headers/BasicAnaFunc.h"

void handy_calc() {
  double p1=1011.8;
  double e1=1000.*p2ke(p1/1000.);

  double p1_mc=1007.9;
  double e1_mc=1000.*p2ke(p1_mc/1000.);


  double p2_stop_range=955.1;
  double e2_stop_range=1000.*p2ke(p2_stop_range/1000.);

  double p2_stop_range_mc=960.2;
  double e2_stop_range_mc=1000.*p2ke(p2_stop_range_mc/1000.);

  cout<<"e1:"<<e1<<endl;	
  cout<<"e1_mc:"<<e1_mc<<endl;	

  cout<<"e2_stop_range:"<<e2_stop_range<<endl;	
  cout<<"e2_stop_range_mc:"<<e2_stop_range_mc<<endl;	

 //
 double pp1=1.0; 
 double pp1_1=1.30;
 double dp_ov_p=100.*(pp1_1-pp1)/pp1;
 cout<<"dp_ov_p:"<<dp_ov_p<<endl;


 double ee1=p2ke(pp1);
 double ee1_1=p2ke(pp1_1);
 double de_ov_e=100.*(ee1_1-ee1)/ee1;
 cout<<"de_ov_e:"<<de_ov_e<<endl;

 //



}
